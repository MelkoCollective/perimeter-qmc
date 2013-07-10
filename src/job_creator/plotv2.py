#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    14.06.2013 12:35:45 EDT
# File:    plotv2.py

import sys
import copy
import matplotlib
matplotlib.use('Agg') #only save /no show /must happen before pylab
import pylab
import numpy as np
import collections as col #deque
import glob
from string import maketrans #maketrans

from addon import *

def style_string(s):
    r = s.translate(maketrans("_", " "))
    r2 = r.split("[")
    return r2[0].title() + ("[" + r2[1] if len(r2) > 1 else "");

#------------------- global variables ------------------- 
N = 0      #amount of lines
L = 0      #amount of labels
labels = 0
data = 0
#------------------- read the file and generate the global variables ------------------- 
def read_file(name):
    global N, L, labels, data
    ifs = open(name, "r")
    raw_lines = ifs.readlines()
    ifs.close()
    labels = raw_lines[0][:-1].split(' ') #get labels
    labels.append("index[-]") #add extra label
    lines = [x[:-1] for x in raw_lines[1:] if x[0] != '-'] #remove comments and labels
    
    N = len(lines)
    L = len(labels)
    data = np.zeros((L, N))
    
    for i in range(len(labels)):
        labels[i] = style_string(labels[i])
    
    for i in range(N):
        lines[i] = lines[i].split(' ') + [str(i)]
        for j in range(L):
            data[j][i] = lines[i][j]
#------------------- check if user input is sane and convert string labels to int  ------------------- 
def validate_plot_labels(opts):
    global N, L, labels, data
    p = parameter
    for opt in opts:
        if p.contains(opt):
            p[opt + "_"] = p[opt]
            if is_int(p[opt]):
                if p[opt] < 0 or p[opt] >= L:
                    RED("option -{0} = {1} is too big or too small".format(opt, p[opt]))
                    return False
            else:
                if p[opt] not in labels:
                    RED("option -{0} = {1} is not valid".format(opt, p[opt]))
                    return False
                else:
                    p[opt + "_"] = labels.index(p[opt])
    return True
#------------------- check if begin/end/spacing is sane ------------------- 
def validate_range_labels():
    global N, L, labels, data
    p = parameter
    dx = data[p["x_"]]
    for opt in ["b", "e"]:
        if p.contains(opt):
            if not is_number(p[opt]):
                RED("option -{} = {} must be an integer".format(opt, p[opt]))
                return False
            else:
                if p[opt] < min(dx) or p[opt] > max(dx):
                    RED("option -{} = {} is too big or too small".format(opt, p[opt]))
                    return False
                else:
                    if opt == "b":
                        p[opt+"_"] = np.where(dx >= p[opt])[0][0]
                    else:
                        p[opt+"_"] = np.where(dx <= p[opt])[0][-1] + 1
    
    if not p.contains("b"):
        p["b_"] = 0
    if not p.contains("sp"):
        p["sp_"] = 1
    else:
        min_spacing = float(dx[1] - dx[0])
        p["sp_"] = max(1, int(p["sp"] / min_spacing + 0.5))
    
    return True
#------------------- get the range we want to plot ------------------- 
def get_range(data):
    p = parameter
    
    if not p.contains("e_"):
        return data[p["b_"]:][0::p["sp_"]]
    else:
        return data[p["b_"]:p["e_"]][0::p["sp_"]]
#------------------- the main function that handles everything ------------------- 
def main():
    global N, L, labels, data
    p = parameter
    p["o"] = "plot.svg"
    p["legend"] = "best"
    parameter.read(sys.argv)
    files = p["arg"]
    #------------------- check if filenames are given ------------------- 
    if len(files) < 1:
        RED("no file given... abort")
        return False
    #------------------- check if files are readable ------------------- 
    for f in files:
        if not exists(f):
            if "*" in f:
                files += glob.glob(f).sort()
            else:
                RED("file {0} not readable... abort".format(f))
                return False
    #------------------- read first file ------------------- 
    read_file(files[0])
    validate_plot_labels(["x", "y", "err"])
    #------------------- print available labels ------------------- 
    if p.contains("l"):
        for l in labels:
            GREEN("{0:0>2}: {1}".format(labels.index(l), l))
        return True
    #------------------- plot output if -p is given ------------------- 
    if p.contains("p") and p.contains("x") and p.contains("y"):
        if p.contains("lin"):
            colors = col.deque(["r^", "r-", "b^", "b-", "y^", "y-", "g^", "g-", "m^", "m-", "c^", "c-", "k^", "k-"])
        else:
            colors = col.deque(["r^-", "b^-", "y^-", "g^-", "m^-", "c^-", "k^-"])
        
        ax = pylab.subplot(111)
        
        #------------------- prepare plot infos ------------------- 
        pylab.xlabel(labels[p["x_"]])
        pylab.ylabel(labels[p["y_"]])
        if p.contains("ylog"):
            ax.set_yscale("log")
            pylab.ylabel("log(" + labels[p["x_"]] + ")")
            
        if p.contains("xlog"):
            ax.set_xscale("log")
            pylab.xlabel("log(" + labels[p["x_"]] + ")")
        
        if len(files) == 1:
            textlabel = style_string(f[:-4])
            pylab.title(textlabel)
        else:
            pylab.title("Multi-Plot")
            
        pylab.grid(True)
        
        #------------------- the actual plot loop ------------------- 
        for f in files:
            GREEN("plotting file " + f)
            read_file(f)
            textlabel = style_string(f.split("/")[0])
            #------------------- validate x, y, err, b, e and sp if given ------------------- 
            if not validate_plot_labels(["x", "y", "err"]) or not validate_range_labels():
                return False
            
            pd = [];
            
            
            #------------------- accumulate y and err data ------------------- 
            if p.contains("acc"):
                data[p["y_"]] = np.add.accumulate(data[p["y_"]])
                if p.contains("err_"):
                    errmean = np.mean(data[p["err_"]])
                    #~ data[p["err_"]] = np.add.accumulate(data[p["err_"]])
                    data[p["err_"]] = np.sqrt(np.add.accumulate(np.square(data[p["err_"]])))
            if p.contains("dif"):
                data[p["y_"]] -= data[p["dif"]]
                data[p["err_"]] = np.sqrt(np.square(data[p["err_"]] + np.square(data[p["dif"]+1])))
            
            #~ for i in range(len(data[p["y_"]])):
                #~ print(str(data[p["y_"]][i]) + " " + str(data[p["err_"]][i]))
            
            for l in range(L):
                pd.append(get_range(data[l])) #get the datarange that needs to be plotted
                
            if p.contains("acc"):
                s = int(len(pd[p["y_"]])/2)
                ofs = open("cross_res.txt", "a")
                ofs.write(" {0} {1} {2} {3} {4} {5} {6}\n".format(errmean
                                                 , pd[p["y_"]][s-1]
                                                 , pd[p["err_"]][s-1]
                                                 , pd[p["y_"]][s]
                                                 , pd[p["err_"]][s]
                                                 , pd[p["y_"]][-1]
                                                 , pd[p["err_"]][-1]))
                ofs.close()
            
            if not p.contains("err"):
                ax.plot(pd[p["x_"]]
                      , pd[p["y_"]]
                      , colors[0]
                      , label = textlabel
                      )
            else:
                ax.errorbar(pd[p["x_"]]
                          , pd[p["y_"]]
                          , pd[p["err_"]]
                          , 0 #no error in x
                          , colors[0]
                          , label = textlabel
                          #, ecolor = "b"
                          )
            colors.rotate(-1)
            
            if p.contains("lin"):
                m,b = np.polyfit(pd[p["x_"]], pd[p["y_"]], 1) 
                x = pd[p["x_"]][1:]
                x[0]=0
                ax.plot(x, m*x+b, colors[0]) 
                colors.rotate(-1)
        #------------------- legend ------------------- 
        if len(files) > 1:
            leg = pylab.legend(loc = p["legend"]
                             , ncol = 3
                             , fancybox = True
                             , prop={'size':12}
                             )
            leg.get_frame().set_alpha(0.5)
        #------------------- parameter box ------------------- 
        else:
            if p.contains("acc"):
                text  = ""
                text += "{0}: {1:.0f}\n".format("H", data[:,-2][2])
                text += "{0}: {1:.0f}\n".format("L", data[:,-2][3])
                text += "{0}: {1:.0f}%\n".format("Accept", 100*data[:,-2][7])
                #~ text += "{0}: {1:.0f}\n".format("Seed", data[:,-2][labels.index("Seed")])
                text += "{0}: {1:.0f}us".format("Update Time", data[:,-2][labels.index("Loop Time[us]")])
                ax.text(0.37
                      , 0.22
                      , text
                      , transform=ax.transAxes
                      , fontsize=12
                      , verticalalignment='top'
                      , bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)
                      )
            if p.contains("lin"):
                text  = ""
                text += "{0}: {1:.3f}\n".format("m", m)
                text += "{0}: {1:.3f}".format("b", b)
                ax.text(0.37
                      , 0.22
                      , text
                      , transform=ax.transAxes
                      , fontsize=12
                      , verticalalignment='top'
                      , bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)
                      )
        if p.contains("acc"):
            pylab.xlim([pd[p["x_"]][0], pd[p["x_"]][-1]])
        
        #~ pylab.xlim(0, 1.2)
        #~ pylab.ylim(-1.2, .9)
            
        pylab.savefig(p["o"])
        GREENB(p["o"] + " is done")
    
    return True

if __name__ == "__main__":
    GREENB("plotv2.py")
    main()
    
