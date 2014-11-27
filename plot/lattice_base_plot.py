#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    27.08.2014 09:46:09 CEST
# File:    svg_text.py

from addon import *
from addon.svg import *
import math

def shift(mat, dx, dy, mult = 1):
    dx = -dx
    dy = -dy
    res = mat[mult*dy:] + mat[:mult*dy]
    
    for i in range(len(res)):
        res[i] = res[i][mult*dx:] + res[i][:mult*dx]
    
    return res

def generate_grid(**kwargs):
    grid = kwargs["grid"]
    mod = kwargs["mod"]
    site_mark = kwargs["site_mark"]
    tile_mark = kwargs["tile_mark"]
    svg = kwargs.get("svg", group())
    space = kwargs.get("space", 1)
    H = int(len(grid) / 2)
    L = int(len(grid[0]) / 2)
    
    sp = 50
    offH = sp*space
    offL = sp*space
    #=================== sqr ===================
    if mod == "sqr":
        dH = sp
        dL = sp
        size = ((L-1) * dL + 2 * offL, (H-1) * dH + 2 * offH)
    #=================== tri ===================
    elif mod == "tri":
        dH = math.sqrt(3) / 2 * sp
        dL = sp
        size = ((L-1 + (H-1)/2) * dL + 2 * offL, (H-1) * dH + 2 * offH)
    elif mod == "hex":
        dH = math.sqrt(3) / 2 * sp
        dL = sp
        size = (2*(L-1) * dL + 2 * offL, (H-1) * dH + 2 * offH)
    
    pos = []
    for h in range(2*(H+1)):
        pos.append([])
        ln = [0, .25, .5, .25][(h+1)%4]
        prop = [[1, 1, .5, .5], [.75, .75, .75, .75], [.5, .5, 1, 1], [.75, .75, .75, .75]][(h+1)%4]
        idx = 0
        for l in range(2*(L+1)):
            if mod == "sqr":
                pos[-1].append((offL +  (l-1) * dL/2.0, offH + (h-1) * dH/2.0))
            elif mod == "tri":
                pos[-1].append((offL +  ((l-1) + (2*H - 1 - h) / 2.0) * dL/2.0, offH + (h-1) * dH/2.0))
            elif mod == "hex":
                ln += prop[idx]
                
                idx += 1
                idx %= 4
                
                pos[-1].append((offL +  (ln-1) *dL, offH + (h-1) * dH/2.0))
    
    #=================== mark tile/site ===================
    sites = group(stroke = "none", fill_opacity = .5)
    svg.add(sites)
    
    for h in range(H):
        for l in range(L):
            if site_mark[h][l]:
                sites.add(polygon(points = "{},{} {},{} {},{} {},{}".format(pos[2*h][2*l][0], pos[2*h][2*l][1]
                                                                     , pos[2*h+2][2*l][0], pos[2*h+2][2*l][1]
                                                                     , pos[2*h+2][2*l+2][0], pos[2*h+2][2*l+2][1]
                                                                     , pos[2*h][2*l+2][0], pos[2*h][2*l+2][1]
                                                                     ), fill = site_mark[h][l]))
    
    tiles = group(stroke = "none", stroke_width = 0 , fill_opacity = .4) #width = 3
    svg.add(tiles)
    
    for h in range(H):
        for l in range(L):
            if mod == "sqr":
                if tile_mark[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[2*h+1][2*l+1][0], pos[2*h+1][2*l+1][1]
                                                                     , pos[2*h+3][2*l+1][0], pos[2*h+3][2*l+1][1]
                                                                     , pos[2*h+3][2*l+3][0], pos[2*h+3][2*l+3][1]
                                                                     , pos[2*h+1][2*l+3][0], pos[2*h+1][2*l+3][1]
                                                                     ), fill = tile_mark[h][l], stroke = tile_mark[h][l]))
            elif mod == "tri":
                if tile_mark[h][l][0]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[2*h+1][2*l+1][0], pos[2*h+1][2*l+1][1]
                                                                     , pos[2*h-1][2*l+1][0], pos[2*h-1][2*l+1][1]
                                                                     , pos[2*h+1][2*l+3][0], pos[2*h+1][2*l+3][1]
                                                                     , pos[2*h+3][2*l+3][0], pos[2*h+3][2*l+3][1]
                                                                     ), fill = tile_mark[h][l][0], stroke = tile_mark[h][l][0]))
                if tile_mark[h][l][1]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[2*h+1][2*l+1][0], pos[2*h+1][2*l+1][1]
                                                                     , pos[2*h-1][2*l+1][0], pos[2*h-1][2*l+1][1]
                                                                     , pos[2*h-1][2*l-1][0], pos[2*h-1][2*l-1][1]
                                                                     , pos[2*h+1][2*l-1][0], pos[2*h+1][2*l-1][1]
                                                                     ), fill = tile_mark[h][l][1], stroke = tile_mark[h][l][1]))
                if tile_mark[h][l][2]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[2*h+1][2*l+1][0], pos[2*h+1][2*l+1][1]
                                                                     , pos[2*h+1][2*l-1][0], pos[2*h+1][2*l-1][1]
                                                                     , pos[2*h+3][2*l+1][0], pos[2*h+3][2*l+1][1]
                                                                     , pos[2*h+3][2*l+3][0], pos[2*h+3][2*l+3][1]
                                                                     ), fill = tile_mark[h][l][2], stroke = tile_mark[h][l][2]))
            elif mod == "hex":
                if tile_mark[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{} {},{} {},{}".format(
                                                                       pos[2*h+1][2*l+1][0], pos[2*h+1][2*l+1][1]
                                                                     , pos[2*h-1][2*l+1][0], pos[2*h-1][2*l+1][1]
                                                                     , pos[2*h-1][2*l+3][0], pos[2*h-1][2*l+3][1]
                                                                     , pos[2*h+1][2*l+3][0], pos[2*h+1][2*l+3][1]
                                                                     , pos[2*h+3][2*l+3][0], pos[2*h+3][2*l+3][1]
                                                                     , pos[2*h+3][2*l+1][0], pos[2*h+3][2*l+1][1]
                                                                     ), fill = tile_mark[h][l], stroke = tile_mark[h][l]))
            

    
    #=================== draw grid ===================
    grid_line = group(stroke = "black", stroke_width = 1, stroke_dasharray = "2, 2")
    svg.add(grid_line)
    
    start = []
    end = []
    
    if mod in ["sqr", "tri"]:
        for h in range(H):
            grid_line.add(line(xy1 = pos[2*h+1][1], xy2 = pos[2*h+1][-3]))
        for l in range(L):
            grid_line.add(line(xy1 = pos[1][2*l+1], xy2 = pos[-3][2*l+1]))
    
    if mod == "tri":
        for h in range(H):
            if H-h < L: 
                grid_line.add(line(xy1 = pos[2*h+1][1], xy2 = pos[-3][2*(H-h)-1]))
            else:
                grid_line.add(line(xy1 = pos[2*h+1][1], xy2 = pos[2*(h+L)-1][-3]))
        for l in range(1,L):
            if L-l < H:
                grid_line.add(line(xy1 = pos[1][2*l+1], xy2 = pos[2*(L-l)-1][-3]))
            else:
                grid_line.add(line(xy1 = pos[1][2*l+1], xy2 = pos[-3][2*(l+H)-1]))
            
    if mod == "hex":
        for h in range(H):
            for l in range(L):
                if (l + h) % 2:
                    if(l < L-1):
                        grid_line.add(line(xy1 = pos[2*h+1][2*l+1], xy2 = pos[2*h+1][2*l+3]))
                else:
                    if(h < H-1):
                        grid_line.add(line(xy1 = pos[2*h+1][2*l+1], xy2 = pos[2*h+3][2*l+1]))
                    if(h > 0):
                        grid_line.add(line(xy1 = pos[2*h+1][2*l+1], xy2 = pos[2*h-1][2*l+1]))
        
    #=================== draw bonds ===================
    bonds = group(stroke_width = 8)  #, stroke_linecap = "round"
    svg.add(bonds)
    pbc_style = "{}, {}, {}, {}, {}".format(sp/3, 2, 2, 2, 2)
    
    for h in range(H):
        for l in range(L):
            p = pos[2*h+1][2*l+1]
            p_right = pos[2*h + 1][2*l + 3] #no overflow since pos larger than H, L
            p_down = pos[2*h + 3][2*l + 1]
            
            #------------------- horizontal ------------------- 
            if grid[2*h][2*l + 1]:
                col = grid[2*h][2*l + 1]
                if l == L - 1: # pbc
                    bonds.add(line(xy1 = p, xy2 = pos[2*h+1][2*l+2], stroke = col, stroke_dasharray = pbc_style))
                    bonds.add(line(xy1 = pos[2*h+1][1], xy2 = pos[2*h+1][0], stroke = col, stroke_dasharray = pbc_style))
                else: # normal
                    bonds.add(line(xy1 = p, xy2 = p_right, stroke = col))
            #------------------- vertical ------------------- 
            if grid[2*h + 1][2*l]:
                col = grid[2*h + 1][2*l]
                if h == H - 1: # pbc
                    bonds.add(line(xy1 = p, xy2 = pos[2*h+2][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                    bonds.add(line(xy1 = pos[1][2*l+1], xy2 = pos[0][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                else: # normal
                    bonds.add(line(xy1 = p, xy2 = p_down, stroke = col))
            if mod == "tri":
                p_ddown = pos[2*h+3][2*l+3]
                #------------------- diagonal ------------------- 
                if grid[2*h + 1][2*l+1]:
                    col = grid[2*h + 1][2*l+1]
                    if h == H - 1 or l == L - 1: # pbc
                        h_other = (h + 1) % H
                        l_other = (l + 1) % L
                        print(h, l, h_other, l_other)
                        bonds.add(line(xy1 = p, xy2 = pos[2*h+2][2*l+2], stroke = col, stroke_dasharray = pbc_style))
                        bonds.add(line(xy1 = pos[2*h_other + 1][2*l_other+1], xy2 = pos[2*h_other][2*l_other], stroke = col, stroke_dasharray = pbc_style))
                    else: # normal
                        bonds.add(line(xy1 = p, xy2 = p_ddown, stroke = col))
    
    #=================== add sites ===================
    sites = group(stroke = "black", stroke_width = 1)
    svg.add(sites)
    for h in range(H):
        for l in range(L):
            if grid[2*h][2*l]:
                sites.add(circle(cxy = pos[2*h+1][2*l+1], r = 10, fill = grid[2*h][2*l]))
        
    return svg, size

def tri_1():
    mod = "tri"
    
    U = "red"
    D = "blue"
    x = None
    m = "magenta"
    c = "cyan"
    r = "red"
    b = "blue"
    y = "yellow"
    u = "orange"
    
    grid = [
             [U,x,D,x,U,x,D,x]
    ,       [m,x,m,x,x,m,x,x]
    ,      [D,x,U,x,D,x,D,x]
    ,     [x,x,x,x,x,m,x,x]
    ,    [U,x,D,x,U,x,U,x]
    ,   [m,x,m,x,m,x,x,x]
    ,  [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,m,x]
    ]
    site_mark = [
          [0,0,0,0]
    ,    [0,0,0,0]
    ,   [0,0,0,0]
    ,  [0,0,0,0]
    ]
    tile_mark = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,u,0],[u,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,  [[0,0,0],[0,u,0],[0,u,0],[0,0,0]]
    ]
    sqr_tile_mark = [
          [0,r,r,0]
    ,    [0,r,r,0]
    ,   [0,0,0,0]
    ,  [0,0,0,0]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    return mod, grid, site_mark, tile_mark

def tri_neg():
    mod = "tri"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "magenta"
    c = "cyan"
    b = "green"
    u = "orange"
    
    grid = [
             [U,b,D,x,U,x,D,x]
    ,       [x,x,x,x,b,x,b,x]
    ,      [D,c,U,x,D,x,U,x]
    ,     [m,x,x,m,x,x,x,x]
    ,    [U,x,U,c,D,x,D,x]
    ,   [c,x,x,m,x,x,b,x]
    ,  [D,m,U,c,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ]
    site_mark = [
          [0,0,0,0]
    ,    [0,0,0,0]
    ,   [0,0,0,0]
    ,  [0,0,0,0]
    ]
    tile_mark = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    return mod, grid, site_mark, tile_mark

def sqr_1():
    mod = "sqr"
    
    U = "red"
    D = "blue"
    x = None
    m = "magenta"
    c = "cyan"
    r = "red"
    b = "blue"
    y = "yellow"
    u = "orange"
    #~ c = None
    #~ o = None
    grid = [
      [U,x,U,x,U,m,D,x]
    , [m,0,m,0,x,0,x,0]
    , [D,x,D,x,D,m,U,x]
    , [x,0,x,0,x,0,x,0]
    , [U,m,D,x,U,m,D,x]
    , [x,0,x,0,x,0,x,0]
    , [U,m,D,x,U,m,D,x]
    , [x,0,x,0,x,0,x,0]
    ]
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    tile_mark = [
      [u,0,u,0]
    , [0,0,u,0]
    , [u,0,0,0]
    , [0,0,0,0]
    ]
    
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    return mod, grid, site_mark, tile_mark
    
def example():
    #~ mod, grid, site_mark, tile_mark = sqr_1()
    mod, grid, site_mark, tile_mark = tri_1()
    #~ mod, grid, site_mark, tile_mark = tri_neg()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))

    frame1 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod)
    frame2 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod)
    frame3 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod)
    frame4 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod)
    
    pic.add([[frame1], [frame2]])
    
    
    pic.write_svg("grid.svg")
    pic.display()

def tri_base():
    #------------------- grid ------------------- 
    mod = "tri"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "magenta"
    
    
    grid = [
              [U,x,U,x,U,x,U,x,U,x]
    ,        [x,x,x,x,x,x,x,x,x,x]
    ,       [U,x,U,x,U,x,U,x,U,x]
    ,      [x,x,x,x,x,x,x,x,x,x]
    ,     [U,x,U,x,D,x,U,x,U,x]
    ,    [x,x,x,x,x,x,x,x,x,x]
    ,   [U,x,U,x,U,x,U,x,U,x]
    ,  [x,x,x,x,x,x,x,x,x,x]
    , [U,x,U,x,U,x,U,x,U,x]
    ,[x,x,x,x,x,x,x,x,x,x]
    ]
    site_mark = [
          [0,0,0,0,0]
    ,    [0,0,0,0,0]
    ,   [0,0,0,0,0]
    ,  [0,0,0,0,0]
    , [0,0,0,0,0]
    ]
    g = "purple"
    #~ g = None
    c = "green"
    c = None
    u = "orange"
    u = None
    tile_mark1 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[g,0,u],[0,0,u],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[g,0,0],[g,c,u],[0,c,u],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[g,c,0],[0,c,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    g = "purple"
    g = None
    c = "green"
    #~ c = None
    u = "orange"
    u = None
    tile_mark2 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[g,0,u],[0,0,u],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[g,0,0],[g,c,u],[0,c,u],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[g,c,0],[0,c,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    g = "purple"
    g = None
    c = "green"
    c = None
    u = "orange"
    #~ u = None
    tile_mark3 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[g,0,u],[0,0,u],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[g,0,0],[g,c,u],[0,c,u],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[g,c,0],[0,c,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    g = "purple"
    #~ g = None
    c = "green"
    #~ c = None
    u = "orange"
    b = "blue"
    #~ u = None
    tile_mark4 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[g,c,u],[0,0,0],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    
    tile_mark5 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[g,0,0],[g,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[b,0,0],[0,0,0],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[g,0,0],[g,0,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    tile_mark6 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[b,c,0],[0,c,0],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[0,0,0],[0,c,0],[0,c,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    tile_mark7 = [
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,0,0],[0,0,u],[0,0,u],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[b,0,u],[0,0,u],[0,0,0]]
    ,  [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    , [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark1, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark2, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark3, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark4, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame1 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark1, mod = mod)
    frame2 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark2, mod = mod)
    frame3 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark3, mod = mod)
    frame4 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark4, mod = mod)
    frame5 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark5, mod = mod)
    frame6 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark6, mod = mod)
    frame7 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark7, mod = mod)
    
    pic.add([[frame1, frame2, frame3]])
    
    pic.write_svg("svg/tri_all_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame4]])
    pic.write_svg("svg/tri_base_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame5, frame6, frame7]])
    pic.write_svg("svg/tri_nbr_tile.svg")

def sqr_base():
    #------------------- grid ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "magenta"
    
    
    grid = [
              [U,x,U,x,U,x,U,x,U,x]
    ,        [x,x,x,x,x,x,x,x,x,x]
    ,       [U,x,U,x,U,x,U,x,U,x]
    ,      [x,x,x,x,x,x,x,x,x,x]
    ,     [U,x,U,x,D,x,U,x,U,x]
    ,    [x,x,x,x,x,x,x,x,x,x]
    ,   [U,x,U,x,U,x,U,x,U,x]
    ,  [x,x,x,x,x,x,x,x,x,x]
    , [U,x,U,x,U,x,U,x,U,x]
    ,[x,x,x,x,x,x,x,x,x,x]
    ]
    site_mark = [
          [0,0,0,0,0]
    ,    [0,0,0,0,0]
    ,   [0,0,0,0,0]
    ,  [0,0,0,0,0]
    , [0,0,0,0,0]
    ]
    
    g = "green"
    b = "blue"
    tile_mark2 = [
          [0, 0, 0, 0, 0]
    ,    [0, 0, 0, 0, 0]
    ,   [0, 0, g, 0, 0]
    ,  [0, 0, 0, 0, 0]
    , [0, 0, 0, 0, 0]
    ]
    tile_mark1 = [
          [0, 0, 0, 0, 0]
    ,    [0, g, g, 0, 0]
    ,   [0, g, g, 0, 0]
    ,  [0, 0, 0, 0, 0]
    , [0, 0, 0, 0, 0]
    ]
    tile_mark3 = [
          [0, 0, 0, 0, 0]
    ,    [0, 0, g, 0, 0]
    ,   [0, g, b, g, 0]
    ,  [0, 0, g, 0, 0]
    , [0, 0, 0, 0, 0]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark1, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark2, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark3, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame1 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark1, mod = mod)
    frame2 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark2, mod = mod)
    frame3 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark3, mod = mod)
    
    pic.add([[frame1]])
    
    pic.write_svg("svg/sqr_all_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame2]])
    pic.write_svg("svg/sqr_base_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame3]])
    pic.write_svg("svg/sqr_nbr_tile.svg")

def hex_base():
    #------------------- grid ------------------- 
    mod = "hex"
    
    U = "Red"
    D = "Blue"
    B = "Gray"
    x = None
    m = "magenta"
    m = None
    
    grid = [
      [U,x,U,m,U,x,U,m]
    , [m,x,m,x,m,x,m,x]
    , [U,m,U,x,U,m,U,x]
    , [m,x,m,x,m,x,m,x]
    , [U,x,B,m,U,x,U,m]
    , [m,x,m,x,m,x,m,x]
    , [U,m,D,x,U,m,U,x]
    , [m,x,m,x,m,x,m,x]
    , [U,x,U,m,U,x,U,m]
    , [m,x,m,x,m,x,m,x]
    , [U,m,U,x,U,m,U,x]
    , [m,x,m,x,m,x,m,x]
    , [U,x,U,m,U,x,U,m]
    , [m,x,m,x,m,x,m,x]
    , [U,m,U,x,U,m,U,x]
    , [m,x,m,x,m,x,m,x]
    ]
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    
    g = "green"
    b = "blue"
    tile_mark1 = [
      [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [g, 0, 0, 0]
    , [0, g, 0, 0]
    , [g, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    ]
    tile_mark2 = [
      [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, g, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    ]
    tile_mark3 = [
      [0, 0, 0, 0]
    , [0, g, 0, 0]
    , [g, 0, g, 0]
    , [0, b, 0, 0]
    , [g, 0, g, 0]
    , [0, g, 0, 0]
    , [0, 0, 0, 0]
    , [0, 0, 0, 0]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark1, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark2, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark3, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame2 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark2, mod = mod)
    grid[4][2] = U
    frame1 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark1, mod = mod)
    frame3 = generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark3, mod = mod)
    
    pic.add([[frame1]])
    
    pic.write_svg("svg/hex_all_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add(lingrad(id = "Gray", color1 = "black", color2 = "gray", dir = "d"))
    pic.add([[frame2]])
    pic.write_svg("svg/hex_base_tile.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame3]])
    pic.write_svg("svg/hex_nbr_tile.svg")

def tri_transitions():
    #------------------- grid ------------------- 
    mod = "tri"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "red"
    p = "purple"
    g = "green"
    o = "orange"
    
    # p/g/o
    
    grid = []
    
    grid.append([
              [U,x,D,x,D,m,U,x]
    ,        [m,x,m,x,x,x,x,x]
    ,       [D,x,U,x,U,x,U,x]
    ,      [x,x,x,x,m,x,m,x]
    ,     [D,m,U,x,D,x,D,x]
    ,    [x,x,x,x,x,x,x,x]
    ])
    grid.append([
              [U,m,D,x,D,m,U,x]
    ,        [x,x,x,x,x,x,x,x]
    ,       [D,m,U,x,U,x,D,x]
    ,      [x,x,x,x,m,x,m,x]
    ,     [D,m,U,x,D,x,U,x]
    ,    [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
              [U,x,D,x,D,m,U,x]
    ,        [m,x,x,m,x,x,x,x]
    ,       [D,x,U,x,U,x,D,x]
    ,      [x,x,x,m,x,x,m,x]
    ,     [D,m,U,x,D,x,U,x]
    ,    [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
              [U,x,D,x,D,m,U,x]
    ,        [m,x,m,x,x,x,x,x]
    ,       [D,x,U,x,U,m,U,x]
    ,      [x,x,x,x,x,x,x,x]
    ,     [D,m,U,x,D,m,D,x]
    ,    [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = [
          [0,0,0,0]
    ,    [0,0,0,0]
    ,   [0,0,0,0]
    ]
    tile_mark = []
    tile_mark.append([
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,g,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[g,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
          [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,    [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    ,   [[0,0,0],[0,0,0],[0,0,0],[0,g,0]]
    ])
    tile_mark.append(tile_mark[-1])
    
    shiftX = 0
    shiftY = 0
    for g, i_g in zipi(grid):
        grid[i_g] = shift(g, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    for t, i_t in zipi(tile_mark):
        tile_mark[i_t] = shift(t, shiftX, shiftY, mult = 1)
    
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame = []
    for g, i_g in zipi(grid):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[i_g], mod = mod))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    pfeil = [arw, (40, 30)]
    
    pic.add([[frame[0], pfeil, frame[1]]])
    
    pic.write_svg("svg/tri_trans_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2], pfeil, frame[3]]])
    pic.write_svg("svg/tri_trans_2.svg")
    
    arw.add(line(xy1 = (0, 0), xy2 = (40, 30), stroke = "red", stroke_width = 3))
    arw.add(line(xy1 = (40, 0), xy2 = (0, 30), stroke = "red", stroke_width = 3))
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[4], pfeil, frame[5]]])
    pic.write_svg("svg/tri_trans_3.svg")
    
def sqr_transitions():
    #------------------- grid ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "red"
    p = "purple"
    g = "green"
    o = "orange"
    
    grid = []
    
    grid.append([
      [U,x,D,x,D,x,U,x]
    , [m,x,m,x,m,x,m,x]
    , [D,x,U,x,U,x,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,D,m,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,m,D,x,D,x,U,x]
    , [x,x,x,x,m,x,m,x]
    , [D,m,U,x,U,x,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,D,m,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,x,D,x,D,m,U,x]
    , [m,x,m,x,x,x,x,x]
    , [D,x,U,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,D,m,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,x,D,m,D,x,U,x]
    , [m,x,x,x,x,x,m,x]
    , [D,x,U,m,U,x,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,D,m,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    tile_mark = []
    tile_mark.append([
        [g,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
        [0,0,g,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
        [0,g,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    
    shiftX = 0
    shiftY = 0
    for g, i_g in zipi(grid):
        grid[i_g] = shift(g, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    for t, i_t in zipi(tile_mark):
        tile_mark[i_t] = shift(t, shiftX, shiftY, mult = 1)
    
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame = []
    for g, i_g in zipi(grid):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[i_g], mod = mod))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    pfeil = [arw, (40, 30)]
    
    pic.add([[frame[0], pfeil, frame[1]]])
    
    pic.write_svg("svg/sqr_trans_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2], pfeil, frame[3]]])
    pic.write_svg("svg/sqr_trans_2.svg")
    
    arw.add(line(xy1 = (0, 0), xy2 = (40, 30), stroke = "red", stroke_width = 3))
    arw.add(line(xy1 = (40, 0), xy2 = (0, 30), stroke = "red", stroke_width = 3))
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[4], pfeil, frame[5]]])
    pic.write_svg("svg/sqr_trans_3.svg")
    
def hex_transitions():
    #------------------- grid ------------------- 
    mod = "hex"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "red"
    p = "purple"
    g = "green"
    o = "orange"
    
    grid = []
    
    grid.append([
      [U,x,D,m,U,x,D,m]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,x,U,x,D,x]
    , [m,x,m,x,m,x,m,x]
    , [U,x,D,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,D,m,U,x,D,m]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,U,x,D,x]
    , [x,x,x,x,m,x,m,x]
    , [U,x,D,x,D,x,U,x]
    , [m,x,m,x,x,x,x,x]
    , [D,x,U,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,x,D,m,U,x,D,m]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,x,U,m,D,x]
    , [m,x,m,x,x,x,x,x]
    , [U,x,D,x,D,x,U,x]
    , [x,x,x,x,m,x,m,x]
    , [D,m,U,x,U,x,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,x,D,x,U,x,D,m]
    , [x,x,m,x,m,x,x,x]
    , [D,x,U,x,U,x,D,x]
    , [m,x,x,x,x,x,m,x]
    , [U,x,D,m,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [D,m,U,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    tile_mark = []
    tile_mark.append([
        [0,0,0,0]
    ,   [0,0,0,0]
    ,   [g,0,0,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
        [0,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,g,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    tile_mark.append([
        [0,0,0,0]
    ,   [0,g,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    tile_mark.append(tile_mark[-1])
    
    shiftX = 0
    shiftY = 0
    for g, i_g in zipi(grid):
        grid[i_g] = shift(g, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    for t, i_t in zipi(tile_mark):
        tile_mark[i_t] = shift(t, shiftX, shiftY, mult = 1)
    
    
    #------------------- paint ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    
    frame = []
    for g, i_g in zipi(grid):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[i_g], mod = mod))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    pfeil = [arw, (40, 30)]
    
    pic.add([[frame[0], pfeil, frame[1]]])
    
    pic.write_svg("svg/hex_trans_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2], pfeil, frame[3]]])
    pic.write_svg("svg/hex_trans_2.svg")
    
    arw.add(line(xy1 = (0, 0), xy2 = (40, 30), stroke = "red", stroke_width = 3))
    arw.add(line(xy1 = (40, 0), xy2 = (0, 30), stroke = "red", stroke_width = 3))
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[4], pfeil, frame[5]]])
    pic.write_svg("svg/hex_trans_3.svg")
    
def sqr_transition_graph():
    #------------------- grid ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    o = "orange"
    
    grid = []
    
    grid.append([
      [U,x,D,x,U,x,D,x]
    , [m,x,m,x,m,x,m,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,m,D,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,D,x,U,x,D,x]
    , [b,x,b,x,b,x,b,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,D,x,U,x,D,x]
    , [p,x,p,x,p,x,p,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,p,D,x,U,p,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,x,D,b,U,x]
    , [b,x,b,x,x,x,x,x]
    , [U,x,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [m,x,m,x,m,x,m,x]
    , [D,x,U,x,D,b,U,x]
    , [b,x,b,x,x,x,x,x]
    , [U,m,D,x,U,p,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append(grid[0])
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,b,D,x,U,x]
    , [b,x,x,x,x,x,b,x]
    , [U,x,D,b,U,x,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [m,x,m,x,m,x,m,x]
    , [D,x,U,b,D,x,U,x]
    , [b,x,x,x,x,x,b,x]
    , [U,m,D,b,U,m,D,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    tile_mark = []
    tile_mark.append([
        [0,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    
    shiftX = 0
    shiftY = 0
    for g, i_g in zipi(grid):
        grid[i_g] = shift(g, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    for t, i_t in zipi(tile_mark):
        tile_mark[i_t] = shift(t, shiftX, shiftY, mult = 1)
    
    
    #------------------- paint ------------------- 
    
    frame = []
    for g, i_g in zipi(grid):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[0], mod = mod))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    
    pfeil = [arw, (40, 30)]
    plus = [plus_i, (30, 30)]
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0], plus, frame[1], pfeil, frame[2]]])
    pic.write_svg("svg/sqr_trans_graph_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[3], plus, frame[4], pfeil, frame[5]]])
    pic.write_svg("svg/sqr_trans_graph_2.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[6], plus, frame[7], pfeil, frame[8]]])
    pic.write_svg("svg/sqr_trans_graph_3.svg")

def sqr_mapping():
    #------------------- grid ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    x = None
    m = "red"
    r = m
    b = "blue"
    p = "#FF00FF"
    g = "green"
    o = "orange"
    
    grid = []
    
    grid.append([
      [x,x,x,m,x,x]
    , [m,x,x,x,x,x]
    , [x,x,x,m,x,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [x,x,x,b,x,x]
    , [b,x,x,x,x,x]
    , [x,x,x,b,x,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,U,x,D,x]
    , [x,x,x,x,x,x]
    , [D,x,D,x,U,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,U,b,D,x]
    , [b,x,x,x,x,x]
    , [D,x,D,b,U,x]
    , [x,x,x,x,x,x]
    ])
    #------------------- inner product ------------------- 
    grid.append([
      [U,x,U,x,D,x]
    , [x,x,x,x,x,x]
    , [D,x,D,x,U,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,U,x,D,x]
    , [x,x,x,x,x,x]
    , [D,x,D,x,U,x]
    , [x,x,x,x,x,x]
    ])
        #------------------- Zc ------------------- 
        #------------------- Vc ------------------- 
    grid.append([
      [U,x,U,r,D,x]
    , [r,x,x,x,x,x]
    , [D,x,D,r,U,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,U,x,D,x]
    , [b,x,b,x,b,x]
    , [D,x,D,x,U,x]
    , [x,x,x,x,x,x]
    ])
    grid.append([
      [U,x,U,r,D,x]
    , [p,x,b,x,b,x]
    , [D,x,D,r,U,x]
    , [x,x,x,x,x,x]
    ])
    
    
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    tile_mark = []
    tile_mark.append([
        [0,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    
    #------------------- paint ------------------- 
    
    frame = []
    for g, i_g in zipi(grid[:4]):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[0], mod = mod, space = 0.3))
    
    for g, i_g in zipi(grid[4:]):
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark[0], mod = mod))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    
    same_i = group(stroke_width = 7, stroke = "black")
    same_i.add(line(xy1 = (0, 10), xy2 = (30, 10)))
    same_i.add(line(xy1 = (0, 20), xy2 = (30, 20)))
    
    pfeil = [arw, (40, 30)]
    same = [same_i, (30, 30)]
    plus = [plus_i, (30, 30)]
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0]]])
    pic.write_svg("svg/sqr_mapping_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[1]]])
    pic.write_svg("svg/sqr_mapping_2.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2]]])
    pic.write_svg("svg/sqr_mapping_3.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[3]]])
    pic.write_svg("svg/sqr_mapping_4.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[4], same, frame[5]]])
    pic.write_svg("svg/sqr_mapping_5.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[6], plus, frame[7], pfeil, frame[8]]])
    pic.write_svg("svg/sqr_mapping_6.svg")

def shift_growth():
    #------------------- grid ------------------- 
    mod = "sqr"
    
    U = "gray"
    D = "Blue"
    x = None
    m = "red"
    r = m
    b = "blue"
    p = "#FF00FF"
    g = "green"
    o = "orange"
    
    grid = []
    
    grid.append([
      [U,x,U,x,U,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,x,U,x,U,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,x,U,x,U,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = []
    
    site_mark.append([
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ])
    site_mark.append([
      [0,0,0,0]
    , [0,0,0,0]
    , [r,0,0,0]
    ])
    site_mark.append([
      [0,0,0,0]
    , [r,0,0,0]
    , [r,0,0,0]
    ])
    site_mark.append([
      [r,0,0,0]
    , [r,0,0,0]
    , [r,0,0,0]
    ])
    site_mark.append([
      [r,0,0,0]
    , [r,0,0,0]
    , [r,r,0,0]
    ])
    site_mark.append([
      [r,0,0,0]
    , [r,r,0,0]
    , [r,r,0,0]
    ])
    site_mark.append([
      [r,r,0,0]
    , [r,r,0,0]
    , [r,r,0,0]
    ])
    site_mark.append([
      [r,r,0,0]
    , [r,r,0,0]
    , [r,r,r,0]
    ])
    site_mark.append([
      [r,r,0,0]
    , [r,r,r,0]
    , [r,r,r,0]
    ])
    site_mark.append([
      [r,r,r,0]
    , [r,r,r,0]
    , [r,r,r,0]
    ])
    site_mark.append([
      [r,r,r,0]
    , [r,r,r,0]
    , [r,r,r,r]
    ])
    site_mark.append([
      [r,r,r,0]
    , [r,r,r,r]
    , [r,r,r,r]
    ])
    site_mark.append([
      [r,r,r,r]
    , [r,r,r,r]
    , [r,r,r,r]
    ])
    
    tile_mark = []
    tile_mark.append([
        [0,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ,   [0,0,0,0]
    ])
    
    #------------------- paint ------------------- 
    
    frame = []
    for s, i_s in zipi(site_mark):
        frame.append(generate_grid(grid = grid[0], site_mark = s, tile_mark = tile_mark[0], mod = mod, space = 1))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    
    plus = [plus_i, (30, 30)]
    
    for s, i_s in zipi(site_mark):
        pic = canvas()
        pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
        pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
        pic.add([[frame[i_s]]])
        pic.write_svg("svg/sqr_shift_growth_{}.svg".format(i_s))
    
    
if __name__ == "__main__":
    #~ example()
    # Produces each three plots: base tile, tile affected by spin change, tile affected by tile update
    tri_base()
    sqr_base()
    hex_base()
    
    # Produces plots about legal transitions
    tri_transitions()
    sqr_transitions()
    hex_transitions()
    
    sqr_transition_graph()
    
    sqr_mapping()
    
    shift_growth()
