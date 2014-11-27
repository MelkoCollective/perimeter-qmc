#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    27.08.2014 09:46:09 CEST
# File:    svg_text.py

from addon import *
from addon.svg import *
import math
import copy

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
    
    H = int(len(grid) / 2)
    L = int(len(grid[0]) / 2)
    
    sp = 50
    offH = sp
    offL = sp
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
    bonds = group(stroke_width = 8)
    svg.add(bonds)
    pbc_style = "{}, {}, {}, {}, {}".format(sp/3, 2, 2, 2, 2)
    
    for h in range(H):
        for l in range(L):
            p = pos[2*h+1][2*l+1]
            p_right = pos[2*h + 1][2*l + 3] #no overflow since pos larger than H, L
            p_down = pos[2*h + 3][2*l + 1]
            
            def draw_helper(h, l, **kwargs):
                bonds.add(line(opacity = 0.5, **kwargs))
                xy1 = kwargs["xy1"]
                xy2 = kwargs["xy2"]
                bcol = kwargs["stroke"]
                arrow = kwargs.get("arrow", None)
                pos = "#FF8000"
                neg = "#0080FF"
                if grid[2*h][2*l] == "Red":
                    col = neg
                else:
                    col = pos
                
                if xy1[1] != xy2[1]:
                    temp = xy2
                    xy2 = xy1
                    xy1 = temp
                    if col == neg:
                        col = pos
                    else:
                        col = neg
                
                if bcol == "#FF00FF":
                    col = pos
                
                
                diff = [xy2[0], xy2[1]]
                diff[0] -= xy1[0]
                diff[1] -= xy1[1]
                
                if arrow == "double":
                    diff[0] *= 2
                    diff[1] *= 2
                
                if arrow == "double_invert":
                    diff[0] *= -2
                    diff[1] *= -2
                    xy1 = xy2

                s = math.sqrt(pow(diff[0], 2) + pow(diff[1], 2))
                diff[0]/=s
                diff[1]/=s
                
                lp = 0.3*s
                sp = (s - lp)/2
                op = 5
                
                orth = [diff[1]*op, -diff[0]*op]
                l = [diff[0]*lp, diff[1]*lp]
                
                start = [xy1[0] + diff[0]*sp, xy1[1] + diff[1]*sp]
                
                p1 = [start[0] + orth[0], start[1] + orth[1]]
                p2 = [start[0] - orth[0], start[1] - orth[1]]
                p3 = [start[0] + l[0], start[1] + l[1]]
                
                if arrow != "no":
                    bonds.add(polyline(stroke = "black", stroke_width = 2, fill = col, opacity = 1, points = "{0},{1} {2},{3} {4},{5}, {0},{1} {2},{3}".format(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1])))
            
            #------------------- horizontal ------------------- 
            if grid[2*h][2*l + 1]:
                col = grid[2*h][2*l + 1]
                if l == L - 1: # pbc
                    draw_helper(h, l, xy1 = p, xy2 = pos[2*h+1][2*l+2], stroke = col, stroke_dasharray = pbc_style, arrow = "double")
                    draw_helper(h, l, xy1 = pos[2*h+1][1], xy2 = pos[2*h+1][0], stroke = col, stroke_dasharray = pbc_style, arrow = "no")
                else: # normal
                    draw_helper(h, l, xy1 = p, xy2 = p_right, stroke = col)
            #------------------- vertical ------------------- 
            if grid[2*h + 1][2*l]:
                col = grid[2*h + 1][2*l]
                if h == H - 1: # pbc
                    draw_helper(h, l, xy1 = p, xy2 = pos[2*h+2][2*l+1], stroke = col, stroke_dasharray = pbc_style, arrow = "no")
                    draw_helper(h, l, xy1 = pos[1][2*l+1], xy2 = pos[0][2*l+1], stroke = col, stroke_dasharray = pbc_style, arrow = "double_invert")
                else: # normal
                    draw_helper(h, l, xy1 = p, xy2 = p_down, stroke = col)
            if mod == "tri":
                p_ddown = pos[2*h+3][2*l+3]
                #------------------- diagonal ------------------- 
                if grid[2*h + 1][2*l+1]:
                    col = grid[2*h + 1][2*l+1]
                    if h == H - 1 or l == L - 1: # pbc
                        h_other = (h + 1) % H
                        l_other = (l + 1) % L
                        draw_helper(h, l, xy1 = p, xy2 = pos[2*h+2][2*l+2], stroke = col, stroke_dasharray = pbc_style, arrow = "no")
                        draw_helper(h, l, xy1 = pos[2*h_other + 1][2*l_other+1], xy2 = pos[2*h_other][2*l_other], stroke = col, stroke_dasharray = pbc_style, arrow = "double")
                    else: # normal
                        draw_helper(h, l, xy1 = p, xy2 = p_ddown, stroke = col)
    
    #=================== add sites ===================
    sites = group(stroke = "black", stroke_width = 1)
    svg.add(sites)
    for h in range(H):
        for l in range(L):
            sites.add(circle(cxy = pos[2*h+1][2*l+1], r = 10, fill = grid[2*h][2*l]))
        
    return svg, size

def tri_sign():
    #------------------- color ------------------- 
    mod = "tri"
    
    U = "Red"
    D = "Blue"
    U_grad = lingrad(id = U, color1 = "red", color2 = "orange", dir = "u")
    D_grad = lingrad(id = D, color1 = "blue", color2 = "#0096FF", dir = "d")
    W = "Orange"
    S = "Green"
    U2_grad = lingrad(id = W, color1 = "yellow", color2 = "orange", dir = "u")
    D2_grad = lingrad(id = S, color1 = "green", color2 = "#0096FF", dir = "d")
    x = None
    r = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    
    #------------------- grid ------------------- 
    grid = [
             [U,p,D,x,U,p,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,r,U,x,U,p,D,x]
    ,     [b,x,b,x,x,x,x,x]
    ,    [U,x,D,r,U,x,D,x]
    ,   [r,x,x,x,b,x,p,x]
    ,  [D,b,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ]
    grid_2 = [
             [U,p,D,x,U,p,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,r,U,x,U,p,D,x]
    ,     [b,x,x,b,x,x,x,x]
    ,    [U,x,U,r,D,x,D,x]
    ,   [r,x,x,b,x,x,p,x]
    ,  [D,b,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ]
    grid_3 = [
             [U,p,D,x,U,p,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,x,U,x,U,p,D,x]
    ,     [b,r,x,p,x,x,x,x]
    ,    [U,x,U,x,D,x,D,x]
    ,   [r,x,x,b,x,x,p,x]
    ,  [D,b,U,r,D,x,U,x]
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
    if mod == "sqr":
        tile_mark = [
              [0,0,0,0]
        ,    [0,0,0,0]
        ,   [0,0,0,0]
        ,  [0,0,0,0]
        ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    frame = []
    frame.append(generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
    frame.append(generate_grid(grid = grid_2, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
    frame.append(generate_grid(grid = grid_3, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
    
    arw_i = group()
    arw_i.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    arw = [arw_i, (40, 30)]
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    plus = [plus_i, (30, 30)]
    
    #------------------- assembly ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0]]])
    pic.write_svg("svg/tri_sign_1.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[1]]])
    pic.write_svg("svg/tri_sign_2.svg")
    #~ pic.display()
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2]]])
    pic.write_svg("svg/tri_sign_3.svg")
    #~ pic.display()

def sqr_sign():
    #------------------- color ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    U_grad = lingrad(id = U, color1 = "red", color2 = "orange", dir = "u")
    D_grad = lingrad(id = D, color1 = "blue", color2 = "#0096FF", dir = "d")
    W = "Orange"
    S = "Green"
    U2_grad = lingrad(id = W, color1 = "yellow", color2 = "orange", dir = "u")
    D2_grad = lingrad(id = S, color1 = "green", color2 = "#0096FF", dir = "d")
    x = None
    r = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    
    #------------------- grid ------------------- 
    grid = [
      [U,p,D,x,U,p,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,r,U,x,U,p,D,x]
    , [b,x,b,x,x,x,x,x]
    , [U,x,D,r,U,x,D,x]
    , [r,x,x,x,b,x,p,x]
    , [D,b,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ]
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    
    tile_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    frame = []
    frame.append(generate_grid(grid = grid, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
    
    arw_i = group()
    arw_i.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    arw = [arw_i, (40, 30)]
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    plus = [plus_i, (30, 30)]
    
    #------------------- assembly ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0]]])
    pic.write_svg("svg/sqr_sign_1.svg")
    #~ pic.display()

def svg():
    #------------------- color ------------------- 
    mod = "sqr"
    
    U = "Red"
    D = "Blue"
    U_grad = lingrad(id = U, color1 = "red", color2 = "orange", dir = "u")
    D_grad = lingrad(id = D, color1 = "blue", color2 = "#0096FF", dir = "d")
    W = "Orange"
    S = "Green"
    U2_grad = lingrad(id = W, color1 = "yellow", color2 = "orange", dir = "u")
    D2_grad = lingrad(id = S, color1 = "green", color2 = "#0096FF", dir = "d")
    x = None
    r = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    
    #------------------- grid ------------------- 
    grid = []
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,b,U,x,D,b,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,b,U,x,D,b,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    grid.append([
      [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,b,D,x,U,b]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,b,U,x,D,b,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
      [U,b,D,x,U,x,D,x]
    , [x,x,x,x,b,x,x,x]
    , [D,x,U,x,D,x,U,b]
    , [x,x,b,x,x,x,x,x]
    , [U,x,D,x,U,b,D,x]
    , [b,x,x,x,x,x,x,x]
    , [D,x,U,b,D,x,U,x]
    , [x,x,x,x,x,x,b,x]
    ])
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    
    tile_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ]
    
    shiftX = 0
    shiftY = 0
    grid = shift(grid, shiftX, shiftY, mult = 2)
    site_mark = shift(site_mark, shiftX, shiftY, mult = 1)
    tile_mark = shift(tile_mark, shiftX, shiftY, mult = 1)
    
    #------------------- paint ------------------- 
    frame = []
    for g in grid:
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
        frame[-1][0].add(line(xy1 = (75, 25), xy2 = (75, 225), stroke = "black", opacity = 1, stroke_width = 1))
        frame[-1][0].add(line(xy1 = (25, 75), xy2 = (225, 75), stroke = "black", opacity = 1, stroke_width = 1))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    plus = [plus_i, (30, 30)]
    
    #------------------- assembly ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0]]])
    pic.write_svg("svg/winding_number_sqr_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[1]]])
    pic.write_svg("svg/winding_number_sqr_2.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2]]])
    pic.write_svg("svg/winding_number_sqr_3.svg")

def winding_number_tri():
    #------------------- color ------------------- 
    mod = "tri"
    
    U = "Red"
    D = "Blue"
    U_grad = lingrad(id = U, color1 = "red", color2 = "orange", dir = "u")
    D_grad = lingrad(id = D, color1 = "blue", color2 = "#0096FF", dir = "d")
    W = "Orange"
    S = "Green"
    U2_grad = lingrad(id = W, color1 = "yellow", color2 = "orange", dir = "u")
    D2_grad = lingrad(id = S, color1 = "green", color2 = "#0096FF", dir = "d")
    x = None
    r = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    
    #------------------- grid ------------------- 
    grid = []
    grid.append([
             [U,r,D,x,U,r,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,x,U,x,D,r,U,x]
    ,     [x,r,x,r,x,x,x,x]
    ,    [U,x,D,x,U,x,D,x]
    ,   [r,x,x,x,x,x,r,x]
    ,  [D,x,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
             [U,p,D,x,U,p,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,r,U,x,U,p,D,x]
    ,     [b,x,x,b,x,x,x,x]
    ,    [U,x,U,r,D,x,D,x]
    ,   [r,x,x,b,x,x,p,x]
    ,  [D,b,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    grid.append([
             [U,p,D,x,U,p,D,x]
    ,       [x,x,x,x,x,x,x,x]
    ,      [D,x,U,x,U,p,D,x]
    ,     [b,r,x,p,x,x,x,x]
    ,    [U,x,U,x,D,x,D,x]
    ,   [r,x,x,b,x,x,p,x]
    ,  [D,b,U,r,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    ])
    
    site_mark = [
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
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
    
    #------------------- paint ------------------- 
    frame = []
    for g in grid:
        frame.append(generate_grid(grid = g, site_mark = site_mark, tile_mark = tile_mark, mod = mod))
        #~ frame[-1][0].add(line(xy1 = (75, 25), xy2 = (75, 225), stroke = "black", opacity = 1, stroke_width = 1))
        #~ frame[-1][0].add(line(xy1 = (25, 75), xy2 = (225, 75), stroke = "black", opacity = 1, stroke_width = 1))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    plus = [plus_i, (30, 30)]
    
    #------------------- assembly ------------------- 
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[0]]])
    pic.write_svg("svg/winding_number_tri_1.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[1]]])
    pic.write_svg("svg/winding_number_tri_2.svg")
    
    pic = canvas()
    pic.add(lingrad(id = "Red", color1 = "red", color2 = "orange", dir = "u"))
    pic.add(lingrad(id = "Blue", color1 = "blue", color2 = "#0096FF", dir = "d"))
    pic.add([[frame[2]]])
    pic.write_svg("svg/winding_number_tri_3.svg")

if __name__ == "__main__":    
    tri_sign()
    sqr_sign()
    svg()
    winding_number_tri()
