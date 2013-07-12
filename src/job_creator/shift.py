#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    24.06.2013 11:17:25 EDT
# File:    shift.py

import sys
from addon import *
import copy
import numpy as np

class periodic(object):
    def __init__(self, iterable):
        self.seq = iterable
        
    def __iter__(self):
        return self.seq.__iter__()
    
    def pr(self):
        print(self.seq)
    
    def __getitem__(self, i):
        return self.seq[(i + len(self.seq)) % len(self.seq)]
    
    def __setitem__(self, i, val):
        self.seq[(i + len(self.seq)) % len(self.seq)] = val
    
    
shift_region = 0;
H = 0
L = 0
maxS = 0
found2 = False

def empty():
    global shift_region, H, L
    for i in range(H):
        for j in range(L):
            if shift_region[i][j] == 1:
                return False
    return True

old_i = 0
old_j = 0

def multiply(mH, mL):
    global shift_region, H, L
    
    newH = mH * H
    newL = mL * L
    
    li = [[0 for i in range(newL)] for j in range(newH)]
    
    for i in range(H):
        for j in range(L):
            for k1 in range(mH):
                if shift_region[i][j] == 2:
                    for k2 in range(mL-1):
                        li[mH*i+k1][mL*j+k2] = min(1, shift_region[i][j])
                    li[mH*i+k1][mL*j+mL-1] = shift_region[i][j]
                else:
                    for k2 in range(mL):
                        li[mH*i+k1][mL*j+k2] = shift_region[i][j]
    
    shift_region = li
    H = newH
    L = newL
    make_periodic()
    
def max_steps(eH = 1, eL = 1):
    global shift_region, H, L, maxS, found2, old_i, old_j, parameter
    
    read_base_file(parameter["shift"])
    multiply(eH/L, eL/L)
    
    temp = copy.deepcopy(shift_region);
    
    while not empty():
        maxS += 1
        grow_step(-1)
    
    shift_region = temp
    old_i = 0
    old_j = 0
    
    return maxS

def make_periodic():
    global shift_region, H, L
    for i in range(H):
        shift_region[i] = periodic(shift_region[i])
    shift_region = periodic(shift_region)


def read_base_file(dir_name):
    global shift_region, H, L
    if not exists(dir_name):
        REDB("the file {0} is not readable. abort...".format(dir_name))
        return
    ifs = open(dir_name)
    shift_region = ifs.readlines()
    
    shift_region = [l[:-1].split() for l in shift_region]
    
    H = len(shift_region)
    L = len(shift_region[0])
    
    for i in range(H):
        for j in range(L):
            shift_region[i][j] = to_number(shift_region[i][j])
    
    make_periodic()
    
def print_shift_reagion():
    for i in range(H):
        for j in range(L):
            if shift_region[i][j]:
                GREENB(shift_region[i][j], "")
            else:
                RED(shift_region[i][j], "")
        print("")
    print("")

def write_shift_reagion(dir_name, mod = "w"):
    ofs = open(dir_name, mod)
    for i in range(H):
        for j in range(L):
             ofs.write(str([0, 1][shift_region[i][j]>0]))
             ofs.write(" ")
        ofs.write("\n")
    if mod == "w":
        ofs.write("\n")
        

def add_existent(num, cond):
    global shift_region, H, L
    for i in range(H):
        for j in range(L):
            if shift_region[i][j] > cond:
                shift_region[i][j] += num;


def grow_step(mod = 1, maxcount = 10):
    global shift_region, H, L, old_i, old_j
    #~ directions = ["right", "left", "up", "down"]
    #~ directions = ["left", "right"]
    #~ directions = ["up"]
    directions = ["right"]
    
    add = 100;

    count = 0
    subcount = 1
    trigger = False
    while count <= maxcount and subcount != 0:
        subcount = 0
        add_existent(add, 0);
        for i in range(old_i, H):
            for j in range(old_j, L):
                if shift_region[i][j] == add + 1:
                    if  shift_region[i  ][j+1] == 0 and "right" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i  ][j+1] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    elif  shift_region[i  ][j-1] == 0 and "left" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i  ][j-1] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    elif  shift_region[i-1][j  ] == 0 and "up" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i-1][j  ] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    elif  shift_region[i+1][j  ] == 0 and "down" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i+1][j  ] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    elif  shift_region[i-1][j-1] == 0 and "diag_up" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i-1][j-1] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    elif  shift_region[i+1][j+1] == 0 and "diag_down" in directions:
                        if count + subcount >= maxcount:
                            old_i = i
                            old_j = j
                            add_existent(-add, add-1);
                            return
                        if mod == 1:
                            shift_region[i+1][j+1] = 1
                        elif mod == -1:
                            shift_region[i][j] = add
                    else:
                        subcount -= 1
                    subcount+=1
                elif shift_region[i][j] == add + 2:
                    if count + subcount >= maxcount:
                        old_i = i
                        old_j = j
                        add_existent(-add, add-1);
                        return
                    shift_region[i][j] = 0
                    subcount+=1
            if old_j:
                old_j = 0
                trigger = True
        if old_i:
            old_i = 0
            trigger = True
        add_existent(-add, add-1);
        count += subcount
        if trigger:
            trigger = False
            subcount = 1
        
        

def sitelist(L, sp):
    sitelist = []
    
    for i in range(int(1/sp)):
        if i >= int(1/sp)/2.:
            sitelist.append(round((i+1)*sp*L) - round(i*sp*L))
        else:
            sitelist.append(round((i+1)*sp*L-0.00001) - round(i*sp*L-0.00001))
    
    add = 0
    for i in sitelist:
        add += i;
    assert(add == L)
    assert(sitelist[0] == sitelist[-1])
    return sitelist

def main():
    p = parameter
    read_base_file(p["arg"][0])
    multiply(p["mH"], p["mL"])
    max_steps()
    
    if p.contains("max"):
        print(str(maxS))
        return
    
def grow(g):
    global parameter, old_i, old_j
    p = parameter
    if p.contains("spacing"):
        s = sitelist(H, p["spacing"])
        
        if int(g) == g:
            for i in range(int(g), L):
                grow_step(-1, H)
            return
        
        for i in range(min(int(g + 1.0000001), L), L):
            grow_step(-1, H)
        g -= int(g)
        g /= p["spacing"]
        g = int(round(g))
        old_i = 0
        old_j = 0
        for i in range(g, int(round(1/p["spacing"]))):
            grow_step(-1, s[i])
    else:
        for i in range(g, maxS):
            grow_step(-1)
    

    
    
    
if __name__ == "__main__":
    print("shift.py")
    p = parameter
    p["mH"] = 1
    p["mL"] = 1
    p["g"] = 0
    
    parameter.read(sys.argv)
    p["arg"].append("shift1.txt")
    main()
