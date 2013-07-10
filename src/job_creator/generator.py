#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    06.06.2013 14:01:12 EDT
# File:    min_plot.py

import sys

from addon import *
import shift

import os
import shutil
import glob
import copy

loc = ord('A') - 1
idx_list = []

#------------------- style the numer with leading 0s ------------------- 
def str_elegant(nr):
    res = ""
    if nr < 10:
        res += "0"
    if nr < 100:
        res += "0"
    
    res += str(int(nr))
    
    if loc > ord('A'):
        dec = str(nr - int(nr)) + "0000";
        res += "_";
        res += dec[2:4]
        
    return res;
#----- visits all folders recursivly and calls fct there (and inbetween if intermed) -----
def recursion(path, dirs, fct, intermed = False):
    global loc
    global idx_list
    loc += 1
    if is_list(dirs[0]) and len(dirs[0]) == 2 and type(dirs[0][0]) != type(dirs[0][1]):
        for i in range(len(dirs)):
            idx_list.append(dirs[i][0]);
            if intermed:
                fct(path + "/" + chr(loc) +"-"+ str_elegant(dirs[i][0]))
            recursion(path + "/" + chr(loc) +"-"+ str_elegant(dirs[i][0]), dirs[i][1], fct, intermed)
            idx_list.pop()
        
        loc -= 1
        return
    
    
    for dir_ in dirs:
        idx_list.append(dir_);
        fct(path + "/"+ chr(loc)+"-" + str_elegant(dir_))
        idx_list.pop()
    loc -= 1
#------------------- make folders -------------------
def mkfolder(path_dir):
    if not exists(path_dir):
        os.mkdir(path_dir)
        GREEN(path_dir + " created")
    else:
        YELLOW(path_dir + " exists already")
#------------------- recursion forward ------------------- 
def mkfolders(path, dirs):
    recursion(path, dirs, mkfolder, True)
#------------------- remove folders -------------------
def rmfolder(path_dir):
    if exists(path_dir):
        shutil.rmtree(path_dir)
        RED(path_dir + " deleted")
    else:
        YELLOW(path_dir + " already deleted")
#------------------- recursion forward -------------------     
def rmfolders(path, dirs):
    recursion(path, dirs, rmfolder, True)
#------------------- copy to folders ------------------- 
def cp_to_folder(data, path_dir):
    if exists(path_dir):
        if exists(data):
            shutil.copy(data, path_dir)
            GREEN(data + " copied to " + path_dir)
        else:
            RED(data + " not found")
    else:
        YELLOW(path_dir + " not found")
#------------------- recursion forward ------------------- 
def cp_to_folders(data, path, dirs):
    if not is_list(data):
        data = [data]
    
    for d in data:
        fct = lambda path_dir: cp_to_folder(d, path_dir)
        recursion(path, dirs, fct)
#------------------- progress ------------------- 
def fancy_progress_bar(p):
    size = 50;
    
    bars = int(p * size);
    bar = "<"

    for i in range(bars):
        bar += "|"
    for i in range(bars, size):
        bar += " "
    bar += "> "
    end = str(int(p*100)) + "%"
    
    if(p < .33):
        RED(bar + NONE_ + end)
    elif(p < .66):
        YELLOW(bar + NONE_ + end)
    elif(p == 1):
        CYAN(bar + NONE_ + end)
    else:
        GREEN(bar + NONE_ + end)
#------------------- state files ------------------- 
def read_state_file_impl(path_dir):
    ofs = open(path_dir + "/state.txt", "r")
    ll = ofs.readlines()
    ofs.close()
    tel = {}
    for l in ll:
        l = l.split();
        tel[l[0]] = to_number(l[1])
    return tel
    
queue = 0
done = 0
run = 0
N = 10
total_p = 0
slowest_job = "00:00:00";
def read_state_files(path, dirs):
    global queue
    global done
    global run
    global total_p
    global slowest_job
    queue = 0
    done = 0
    run = 0
    total_p = 0
    slowest_job = "00:00:00";
    
    def read_state_file(path_dir):
        global queue
        global done
        global run
        global total_p
        global slowest_job
        if exists(path_dir + "/state.txt"):
            fail = True
            while(fail):
                tel = read_state_file_impl(path_dir);
                
                try: #it can happen that the file is read during c++ is writing in it and then p is not found
                    total_p += tel["p"];
                    if(tel["p"] == 1):
                        done += 1
                    else:
                        run += 1
                        if "a" in parameter.keys():
                            CYAN("job in " + path_dir + " is done in " + REDB_ + tel["t"], "")
                            fancy_progress_bar(tel["p"])
                        
                        if tel["t"] > slowest_job:
                            slowest_job = tel["t"];
                    fail = False
                except:
                    pass
        else:
            queue += 1 
    
    recursion(path, dirs, read_state_file)
    
    N = queue + done + run
    
    if queue:
        REDB(["", "all "][N == queue] + str(queue) + " job" + ["s are"," is"][queue == 1] + " still in queue or failed")
    if run:
        YELLOWB(["", "all "][N == run] + str(run) + " job" + ["s are"," is"][run == 1] + " running")
    if done:
        GREENB(["", "all "][N == done] + str(done) +  " job" + ["s are"," is"][done == 1] + " done")
    
    if N == done:
        CYAN("all done in " + path + "                              ", "")
    else:
        CYAN("last running job in " + path + "/* is done in " + REDB_ + slowest_job, "")
    
    fancy_progress_bar(float(total_p) / (queue + run + done))
    
#------------------- collect results -------------------
def collect_results(path, dirs):
    bash("cp addon/header.txt " + path + "/colres.txt")
    #~ fct = lambda path_dir: bash("tail -2 " + path_dir + "/fixres.txt >> " + path + "/colres.txt")
    fct = lambda path_dir: bash("tail -2 " + path_dir + "/results.txt >> " + path + "/colres.txt")
    
    recursion(path, dirs, fct)
#------------------- delete files ------------------- 
def delete_files(path, dirs, files):
    
    for f in files:
        fct = lambda path_dir: bash("rm " + path_dir + "/" + f)
        recursion(path, dirs, fct)
#------------------- bash files ------------------- 
def write_bash_file(path_dir):
    args = parameter["args"];
    ofs = open(path_dir + "/bash_in.txt", "w")
    
    ofs.write(parameter["bash"] + " ")
    
    for i in range(len(args)):
        arg = args[i]
        if is_list(arg):
            for a in arg:
                ofs.write("-" + a + " " + str(idx_list[i]) + " ")
        else:
            ofs.write("-" + arg + " " + str(idx_list[i]) + " ")
    
    ofs.write("\n")
    ofs.close()
    
    
    path_dir + "/shift.txt"
    shift.read_base_file(parameter["shift"])
    shift.max_steps(idx_list[0], idx_list[0])
    #~ shift.multiply(idx_list[0]/shift.H, idx_list[0]/shift.L)
    shift.grow(idx_list[1])
    
    if idx_list[1] != 0:
        bash("cp temp.txt " + path_dir + "/shift.txt")
    else:
        shift.write_shift_reagion(path_dir + "/shift.txt")
        
    shift.write_shift_reagion(path_dir + "/shift.txt", "a")
    shift.write_shift_reagion("temp.txt", "w")
    
    
    CYAN("bash file nr " + str(idx_list) + " written")
#------------------- recursion forward ------------------- 
def write_bash_files(path, dirs):
    recursion(path, dirs, write_bash_file)
    bash("rm temp.txt")
#------------------- prog start -------------------
def launch_program(path_dir, prog):
    if not exists(path_dir + "/state.txt"):
        #~ bash("sqsub -q serial " + parameter["sq"] + " -o run.log " + path_dir + "/" + prog) #sharcnet
        bash(path_dir + "/" + prog + " & ") #macBook
    else:
        tel = read_state_file_impl(path_dir);
        if "p" in tel.keys():
            if tel["p"] == 1:
                GREEN(path_dir + "/" + prog + " is done")
                return
                
        YELLOW(path_dir + "/" + prog + " is restarted")
        #~ bash("sqsub -q serial " + parameter["sq"] + " -o run.log " + path_dir + "/" + prog) #sharcnet
        bash(path_dir + "/" + prog + " & ") #macBook
#------------------- recursion forward ------------------- 
def launch_programs(prog, path, dirs):
    fct = lambda path_dir: launch_program(path_dir, prog)
    
    recursion(path, dirs, fct)
#------------------- if something goes wrong on sharcnet ------------------- 
def kill_programs(id0, end):
    for i in range(end-id0+1):
        bash("sqkill " + str(id0 + i))


def deviders(num):
    res = []
    for i in range(1, num):
        if num % i == 0:
            res.append(i)
    return res

spacing = lambda size: 1./min(deviders(size), key=lambda x: abs(size/x - 8))

def shift_root(ll, root, dirs):
    for idx in ll:
        for s in dirs:
            if s[0] == idx:
                dirs = s[1]
                break;
    temp = dirs;
    parameter["sel"] = root
    for idx in reversed(ll):
        global loc
        loc += 1
        parameter["sel"] += "/" + chr(loc) + "-" + str_elegant(idx)
        temp = [[idx, temp]]
    
    root = parameter["sel"]
    return root, dirs
#=================== main ===================
def main():
    global loc
    CLRSCR("")
    
    parameter.read(sys.argv)
    
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == "-a"):
        
        modules = glob.glob("./*.py")
        modules.remove("./sync.py")
        modules.remove("./plotv2.py")
        for i in range(len(modules)):
            modules[i] = modules[i][2:]
        
        modules.sort()
        BLUEB("progress for all modules")
        for mod in modules:
            if exists(mod[:-3] +  "/"): #only show modules that have a folder 
                execfile(mod)
                root  = mod[:-3]
                dirs  = parameter["dirs"]
                read_state_files(root, dirs)
            
        return
    
    if not exists(sys.argv[1]):
        REDB("module doesn't exist!")
        return
    
    execfile(sys.argv[1])
    
    root  = sys.argv[1][:-3]
    dirs  = parameter["dirs"]
    files = parameter["files"]
    
    #~ sqsub -q serial -r 2d --mpp 1g -o run.log 
    
    if "sel" in parameter.keys():
        ll = eval(parameter["sel"])
        root, dirs = shift_root(ll, root, dirs)
    
    if "split" in parameter.keys():
        split = [[x[0]] for x in dirs]
    else:
        split = [[]]
    
    old_root = copy.deepcopy(root);
    old_dirs = copy.deepcopy(dirs);
    old_loc = loc
    
    for sub in split:
        loc = old_loc
        root, dirs = shift_root(sub, old_root, old_dirs)
        print(root)
    
        bash_if("j", "sqjobs")
        bash_if("l", "cat run.log")
        bash_if("s", lambda:read_state_files(root, dirs))
        
        bash_if("c", lambda:collect_results(root, dirs))
        bash_if("p", "./plotv2.py " + root + "/colres.txt -p -x X -y Entropy -err Error -acc -sp 1 -o " + root + "/plot.svg")
        
        bash_if("clean", lambda:( rmfolder(root)
                                , bash("rm *.log")
                                , bash("cd ../../build;  make clean")
                                ))
        bash_if("comp", "cd ../../build; cmake ../ " + parameter["cmake"] + "; make sim")
        bash_if("make", lambda: ( mkfolder(root)
                                , mkfolders(root, dirs)
                                , cp_to_folders(files, root, dirs)
                                , write_bash_files(root, dirs)
                                ))
        bash_if("del", lambda: ( delete_files(root, dirs, ["sim"]) ))
        bash_if("run", lambda:launch_programs("sim", root, dirs))
        
        bash_if("kill", lambda:kill_programs(1583456, 1583458))
        

if __name__ == "__main__":
    print("min_plot.py")
    main();
