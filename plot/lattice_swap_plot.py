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

global_bond_opacity = 0.6

def shift(mat, dx, dy, mult = 1):
    dx = -dx
    dy = -dy
    res = mat[mult*dy:] + mat[:mult*dy]
    
    for i in range(len(res)):
        res[i] = res[i][mult*dx:] + res[i][:mult*dx]
    
    return res

def generate_grid_2(**kwargs):
    grid = kwargs["grid"]
    grid_2 = kwargs.get("grid_2", grid)
    
    mod = kwargs["mod"]
    site_mark = kwargs["site_mark"]
    site_mark_2 = kwargs.get("site_mark_2", site_mark)
    tile_mark = kwargs["tile_mark"]
    tile_mark_2 = kwargs.get("tile_mark_2", tile_mark)
    
    swap = kwargs["swap"]
    swap_inv = []
    
    for r in swap:
        swap_inv.append([])
        for e in r:
            swap_inv[-1].append(1-e)
    
    svg = kwargs.get("svg", group())
    
    H = int(len(grid) / 2)
    L = int(len(grid[0]) / 2)
    
    sp = 50
    shear = 0.5
    offH = sp
    offL = sp
    offM = sp*2
    #=================== sqr ===================
    if mod == "sqr":
        dH = math.sqrt(pow(1, 2) - pow(shear, 2)) * sp
        dL = sp
        sL = shear*sp
        size = ((L-1) * dL + (H-1) * sL + 2 * offL, 2*(H-1) * dH + 2 * offH + offM)
    
    pos1 = []
    pos2 = []
    for h in range(2*(H+1)):
        pos1.append([])
        pos2.append([])
        for l in range(2*(L+1)):
            if mod == "sqr":
                pos1[-1].append((offL +  (l-1) * dL/2.0 + (H +2 - h) * sL/2.0, offH + (h-1) * dH/2.0))
                pos2[-1].append((offL +  (l-1) * dL/2.0 + (H +2 - h) * sL/2.0, offH + (h-1) * dH/2.0 + offM + (H-1) * dH))
    
    pos = [pos1, pos2]
    
    #=================== mark tile/site ===================
    sites = group(stroke = "none", fill_opacity = .5)
    svg.add(sites)
    
    for h in range(H):
        for l in range(L):
            s = swap[h][l]
            if site_mark[h][l]:
                sites.add(polygon(points = "{},{} {},{} {},{} {},{}".format(pos[s][2*h][2*l][0], pos[s][2*h][2*l][1]
                                                                     , pos[s][2*h+2][2*l][0], pos[s][2*h+2][2*l][1]
                                                                     , pos[s][2*h+2][2*l+2][0], pos[s][2*h+2][2*l+2][1]
                                                                     , pos[s][2*h][2*l+2][0], pos[s][2*h][2*l+2][1]
                                                                     ), fill = site_mark[h][l]))
            s = swap_inv[h][l]
            if site_mark_2[h][l]:
                sites.add(polygon(points = "{},{} {},{} {},{} {},{}".format(pos[s][2*h][2*l][0], pos[s][2*h][2*l][1]
                                                                     , pos[s][2*h+2][2*l][0], pos[s][2*h+2][2*l][1]
                                                                     , pos[s][2*h+2][2*l+2][0], pos[s][2*h+2][2*l+2][1]
                                                                     , pos[s][2*h][2*l+2][0], pos[s][2*h][2*l+2][1]
                                                                     ), fill = site_mark_2[h][l]))
    
    tiles = group(stroke = "none", stroke_width = 0 , fill_opacity = .4) #width = 3
    svg.add(tiles)
    
    for h in range(H):
        for l in range(L):
            if mod == "sqr":
                s = swap[h][l]                
                if tile_mark[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[s][2*h+1][2*l+1][0], pos[s][2*h+1][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+1][0], pos[s][2*h+3][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+3][0], pos[s][2*h+3][2*l+3][1]
                                                                     , pos[s][2*h+1][2*l+3][0], pos[s][2*h+1][2*l+3][1]
                                                                     ), fill = tile_mark[h][l], stroke = tile_mark[h][l]))
                s = swap_inv[h][l]                
                if tile_mark_2[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[s][2*h+1][2*l+1][0], pos[s][2*h+1][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+1][0], pos[s][2*h+3][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+3][0], pos[s][2*h+3][2*l+3][1]
                                                                     , pos[s][2*h+1][2*l+3][0], pos[s][2*h+1][2*l+3][1]
                                                                     ), fill = tile_mark_2[h][l], stroke = tile_mark_2[h][l]))
            

    
    #=================== draw grid ===================
    grid_line = group(stroke = "black", stroke_width = 1, stroke_dasharray = "2, 2")
    svg.add(grid_line)
    
    start = []
    end = []
    
    if mod == "sqr":
        for h in range(H):
            grid_line.add(line(xy1 = pos1[2*h+1][1], xy2 = pos1[2*h+1][-3]))
            grid_line.add(line(xy1 = pos2[2*h+1][1], xy2 = pos2[2*h+1][-3]))
        for l in range(L):
            grid_line.add(line(xy1 = pos1[1][2*l+1], xy2 = pos1[-3][2*l+1]))
            grid_line.add(line(xy1 = pos2[1][2*l+1], xy2 = pos2[-3][2*l+1]))
        
    #=================== draw bonds ===================
    bonds = group(stroke_width = 8, opacity = global_bond_opacity, stroke_linecap = "round")
    svg.add(bonds)
    pbc_style = "{}, {}, {}, {}, {}".format(sp/3, 2, 2, 2, 2)
    
    def helper(swap, grid):
        def sub_helper(sw, g, h, l):
            s = sw[h][l]
            sr = sw[h][(l+1)%L]
            sd = sw[(h+1)%H][l]
            p = pos[s][2*h+1][2*l+1]
            p_right = pos[sr][2*h + 1][2*l + 3] #no overflow since pos larger than H, L
            p_down = pos[sd][2*h + 3][2*l + 1]
            
            #------------------- horizontal ------------------- 
            if g[2*h][2*l + 1]:
                col = g[2*h][2*l + 1]
                if l == L - 1: # pbc
                    bonds.add(line(xy1 = p, xy2 = pos[s][2*h+1][2*l+2], stroke = col, stroke_dasharray = pbc_style))
                    bonds.add(line(xy1 = pos[sr][2*h+1][1], xy2 = pos[sr][2*h+1][0], stroke = col, stroke_dasharray = pbc_style))
                else: # normal
                    bonds.add(line(xy1 = p, xy2 = p_right, stroke = col))
            #------------------- vertical ------------------- 
            if g[2*h + 1][2*l]:
                col = g[2*h + 1][2*l]
                if h == H - 1: # pbc
                    bonds.add(line(xy1 = p, xy2 = pos[s][2*h+2][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                    bonds.add(line(xy1 = pos[sd][1][2*l+1], xy2 = pos[sd][0][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                else: # normal
                    bonds.add(line(xy1 = p, xy2 = p_down, stroke = col))

        for h in range(H):
            for l in range(L):
                for i in range(len(swap)):
                    sub_helper(swap[i], grid[i], h, l)
    
    helper([swap, swap_inv], [grid, grid_2])
    
    #=================== add sites ===================
    sites = group(stroke = "black", stroke_width = 1)
    svg.add(sites)
    for h in range(H):
        for l in range(L):
            if grid[2*h][2*l]:
                sites.add(circle(cxy = pos[0][2*h+1][2*l+1], r = 10, fill = grid[2*h][2*l]))
            if grid_2[2*h][2*l]:
                sites.add(circle(cxy = pos[1][2*h+1][2*l+1], r = 10, fill = grid_2[2*h][2*l]))
        
    return svg, size

def generate_trans(**kwargs):
    bra = kwargs["bra"]
    ket = kwargs["ket"]
    
    ket_1 = ket[0][0]
    ket_2 = ket[0][1]
    bra_1 = bra[0][0]
    bra_2 = bra[0][1]
    
    mod = kwargs["mod"]
    site_mark = ket[1][0]
    site_mark_2 = ket[1][1]
    tile_mark = kwargs["tile_mark"]
    tile_mark_2 = kwargs.get("tile_mark_2", tile_mark)
    
    swap = ket[2]
    swap_inv = []
    bra.append([])
    for i in range(len(swap)):
        swap_inv.append([])
        bra[3].append([])
        for j in range(len(swap[0])):
            swap_inv[-1].append(1-swap[i][j])
            bra[3][-1].append(1-bra[2][i][j])
    
    svg = kwargs.get("svg", group())
    
    H = int(len(ket_1) / 2)
    L = int(len(ket_1[0]) / 2)
    
    sp = 50
    shear = 0.5
    offH = sp
    offL = sp
    offM = sp*2
    #=================== sqr ===================
    if mod == "sqr":
        dH = math.sqrt(pow(1, 2) - pow(shear, 2)) * sp
        dL = sp
        sL = shear*sp
        size = ((L-1) * dL + (H-1) * sL + 2 * offL, 2*(H-1) * dH + 2 * offH + offM)
    
    pos1 = []
    pos2 = []
    for h in range(2*(H+1)):
        pos1.append([])
        pos2.append([])
        for l in range(2*(L+1)):
            if mod == "sqr":
                pos1[-1].append((offL +  (l-1) * dL/2.0 + (H +2 - h) * sL/2.0, offH + (h-1) * dH/2.0))
                pos2[-1].append((offL +  (l-1) * dL/2.0 + (H +2 - h) * sL/2.0, offH + (h-1) * dH/2.0 + offM + (H-1) * dH))
    
    pos = [pos1, pos2]
    
    #=================== mark tile/site ===================
    sites = group(stroke = "none", fill_opacity = .5)
    svg.add(sites)
    
    for h in range(H):
        for l in range(L):
            s = swap[h][l]
            if site_mark[h][l]:
                sites.add(polygon(points = "{},{} {},{} {},{} {},{}".format(pos[s][2*h][2*l][0], pos[s][2*h][2*l][1]
                                                                     , pos[s][2*h+2][2*l][0], pos[s][2*h+2][2*l][1]
                                                                     , pos[s][2*h+2][2*l+2][0], pos[s][2*h+2][2*l+2][1]
                                                                     , pos[s][2*h][2*l+2][0], pos[s][2*h][2*l+2][1]
                                                                     ), fill = site_mark[h][l]))
            s = swap_inv[h][l]
            if site_mark_2[h][l]:
                sites.add(polygon(points = "{},{} {},{} {},{} {},{}".format(pos[s][2*h][2*l][0], pos[s][2*h][2*l][1]
                                                                     , pos[s][2*h+2][2*l][0], pos[s][2*h+2][2*l][1]
                                                                     , pos[s][2*h+2][2*l+2][0], pos[s][2*h+2][2*l+2][1]
                                                                     , pos[s][2*h][2*l+2][0], pos[s][2*h][2*l+2][1]
                                                                     ), fill = site_mark_2[h][l]))
    
    tiles = group(stroke = "none", stroke_width = 0 , fill_opacity = .4) #width = 3
    svg.add(tiles)
    
    for h in range(H):
        for l in range(L):
            if mod == "sqr":
                s = swap[h][l]                
                if tile_mark[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[s][2*h+1][2*l+1][0], pos[s][2*h+1][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+1][0], pos[s][2*h+3][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+3][0], pos[s][2*h+3][2*l+3][1]
                                                                     , pos[s][2*h+1][2*l+3][0], pos[s][2*h+1][2*l+3][1]
                                                                     ), fill = tile_mark[h][l], stroke = tile_mark[h][l]))
                s = swap_inv[h][l]                
                if tile_mark_2[h][l]:
                    tiles.add(polygon(points = "{},{} {},{} {},{} {},{}".format(
                                                                       pos[s][2*h+1][2*l+1][0], pos[s][2*h+1][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+1][0], pos[s][2*h+3][2*l+1][1]
                                                                     , pos[s][2*h+3][2*l+3][0], pos[s][2*h+3][2*l+3][1]
                                                                     , pos[s][2*h+1][2*l+3][0], pos[s][2*h+1][2*l+3][1]
                                                                     ), fill = tile_mark_2[h][l], stroke = tile_mark_2[h][l]))
            

    
    #=================== draw grid ===================
    grid_line = group(stroke = "black", stroke_width = 1, stroke_dasharray = "2, 2")
    svg.add(grid_line)
    
    start = []
    end = []
    
    if mod == "sqr":
        for h in range(H):
            grid_line.add(line(xy1 = pos1[2*h+1][1], xy2 = pos1[2*h+1][-3]))
            grid_line.add(line(xy1 = pos2[2*h+1][1], xy2 = pos2[2*h+1][-3]))
        for l in range(L):
            grid_line.add(line(xy1 = pos1[1][2*l+1], xy2 = pos1[-3][2*l+1]))
            grid_line.add(line(xy1 = pos2[1][2*l+1], xy2 = pos2[-3][2*l+1]))
        
    #=================== draw bonds ===================
    bonds = group(stroke_width = 8, opacity = global_bond_opacity) #, stroke_linecap = "round"
    svg.add(bonds)
    
    pbc_style = "{}, {}, {}, {}, {}".format(sp/3, 2, 2, 2, 2)
    
    def helper(swap_arg, grid_arg):
        all_bonds = {}
        order = []
        def color_mix_helper(**kwargs):
            dims = (kwargs["xy1"], kwargs["xy2"])
            
            if dims in all_bonds.keys():
                order[all_bonds[dims]]["stroke"] = "#FF00FF"
            else:
                all_bonds[dims] = len(order)
                order.append(kwargs)
            
        def sub_helper(sw, g, h, l):
            s = sw[h][l]
            sr = sw[h][(l+1)%L]
            sd = sw[(h+1)%H][l]
            p = pos[s][2*h+1][2*l+1]
            p_right = pos[sr][2*h + 1][2*l + 3] #no overflow since pos larger than H, L
            p_down = pos[sd][2*h + 3][2*l + 1]
            
            #------------------- horizontal ------------------- 
            if g[2*h][2*l + 1]:
                col = g[2*h][2*l + 1]
                if l == L - 1: # pbc
                    #~ bonds.add(line(xy1 = p, xy2 = pos[s][2*h+1][2*l+2], stroke = col, stroke_dasharray = pbc_style))
                    #~ bonds.add(line(xy1 = pos[sr][2*h+1][1], xy2 = pos[sr][2*h+1][0], stroke = col, stroke_dasharray = pbc_style))
                    color_mix_helper(xy1 = p, xy2 = pos[s][2*h+1][2*l+2], stroke = col, stroke_dasharray = pbc_style)
                    color_mix_helper(xy1 = pos[sr][2*h+1][1], xy2 = pos[sr][2*h+1][0], stroke = col, stroke_dasharray = pbc_style)
                else: # normal
                    #~ bonds.add(line(xy1 = p, xy2 = p_right, stroke = col))
                    color_mix_helper(xy1 = p, xy2 = p_right, stroke = col)
            #------------------- vertical ------------------- 
            if g[2*h + 1][2*l]:
                col = g[2*h + 1][2*l]
                if h == H - 1: # pbc
                    #~ bonds.add(line(xy1 = p, xy2 = pos[s][2*h+2][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                    #~ bonds.add(line(xy1 = pos[sd][1][2*l+1], xy2 = pos[sd][0][2*l+1], stroke = col, stroke_dasharray = pbc_style))
                    color_mix_helper(xy1 = p, xy2 = pos[s][2*h+2][2*l+1], stroke = col, stroke_dasharray = pbc_style)
                    color_mix_helper(xy1 = pos[sd][1][2*l+1], xy2 = pos[sd][0][2*l+1], stroke = col, stroke_dasharray = pbc_style)
                else: # normal
                    #~ bonds.add(line(xy1 = p, xy2 = p_down, stroke = col))
                    color_mix_helper(xy1 = p, xy2 = p_down, stroke = col)
            
        for h in range(H):
            for l in range(L):
                for i in range(len(swap_arg)):
                    sub_helper(swap_arg[i], grid_arg[i], h, l)
        
        for v in order:
            bonds.add(line(**v))
        
    helper([swap, swap_inv, bra[2], bra[3]], [ket_1, ket_2, bra_1, bra_2])
    
    #=================== add sites ===================
    sites = group(stroke = "black", stroke_width = 1)
    svg.add(sites)
    
    for h in range(H):
        for l in range(L):
            if bra_1[2*h][2*l]:
                sites.add(circle(cxy = pos[bra[2][h][l]][2*h+1][2*l+1], r = 10, fill = bra_1[2*h][2*l]))
            if bra_2[2*h][2*l]:
                sites.add(circle(cxy = pos[bra[3][h][l]][2*h+1][2*l+1], r = 10, fill = bra_2[2*h][2*l]))
        
    return svg, size

def sqr_transition_graph():
    #------------------- grid ------------------- 
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
    m = "red"
    b = "blue"
    p = "#FF00FF"
    g = "green"
    
    
    bra = [
    #------------------- upper copy ------------------- 
    [[[U,x,D,x,U,x,D,x]
    , [m,x,m,x,m,x,m,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,m,D,x,U,m,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,x,S,x,W,x,S,x]
    , [m,x,m,x,m,x,m,x]
    , [S,x,W,x,S,x,W,x]
    , [x,x,x,x,x,x,x,x]
    , [W,m,S,x,W,m,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_make ------------------- 
    ,[[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ]
    ket = []
    ket.append([
    #------------------- upper copy ------------------- 
    [[[U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]
    , [D,x,U,b,D,x,U,x]
    , [b,x,x,x,x,x,b,x]
    , [U,x,D,b,U,x,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,x,S,x,W,b,S,x]
    , [b,x,b,x,x,x,x,x]
    , [S,x,W,x,S,x,W,x]
    , [x,x,x,x,b,x,b,x]
    , [W,b,S,x,W,x,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_make ------------------- 
    ,[[[g,g,g,g]
    , [g,g,g,g]
    , [g,g,g,g]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    #~ , [0,1,1,0]
    , [0,0,0,0]]
    #~ , [1,1,1,1]]
    ])
    ket.append(copy.deepcopy(ket[0]))
    ket[-1][2] = [[0,0,1,1]
                , [0,0,1,1]
                , [0,0,1,1]]
    ket.append(copy.deepcopy(ket[0]))
    ket[-1][2] = [[0,0,0,0]
                , [0,1,1,0]
                , [0,0,0,0]]
    ket.append(copy.deepcopy(ket[0]))
    ket[-1][2] = [[0,0,0,0]
                , [1,1,1,1]
                , [0,0,0,0]]
    tile_mark = []
    tile_mark.append([
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ])
    
    #------------------- paint ------------------- 
    
    frame = []
    for b in [bra] + ket:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    
    for k in ket:
        frame.append(generate_trans(bra = bra, ket = k, mod = mod, tile_mark = tile_mark[0]))
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    
    pfeil = [arw, (40, 30)]
    plus = [plus_i, (30, 30)]
    
    for i in range(1, len(ket) + 1):
        pic = canvas()
        pic.add(U_grad)
        pic.add(D_grad)
        pic.add(U2_grad)
        pic.add(D2_grad)
        pic.add([[frame[0], plus, frame[i], pfeil, frame[len(ket) + i]]])
        pic.write_svg("svg/sqr_swap_" + str(i) + ".svg")
    
def sqr_update_steps():
    #------------------- grid ------------------- 
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
    
    ### +---------------------------------------------------+
    ### |                   bra                             |
    ### +---------------------------------------------------+
    bra = []
    bra.append([
    #------------------- upper copy ------------------- 
    [[[U,x,D,x,U,x,D,x]
    , [r,x,r,x,r,x,r,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,r,D,x,U,r,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,x,S,x,W,x,S,x]
    , [r,x,r,x,r,x,r,x]
    , [S,x,W,x,S,x,W,x]
    , [x,x,x,x,x,x,x,x]
    , [W,r,S,x,W,r,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_rake ------------------- 
    ,[[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ])
    bra.append([
    #------------------- upper copy ------------------- 
    [[[U,r,D,x,U,x,D,x]
    , [x,x,x,x,r,x,r,x]
    , [D,r,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,x,D,r,U,x,D,r]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,r,S,x,W,r,S,x]
    , [x,x,x,x,x,x,x,x]
    , [S,x,W,x,S,r,W,x]
    , [r,x,r,x,x,x,x,x]
    , [W,x,S,x,W,r,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_rake ------------------- 
    ,[[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ])
    bra.append([
    #------------------- upper copy ------------------- 
    [[[x,r,x,x,x,x,x,x]
    , [x,x,x,x,r,x,r,x]
    , [x,r,x,x,x,x,x,x]
    , [x,x,x,x,x,x,x,x]
    , [x,x,x,r,x,x,x,r]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[x,r,x,x,x,r,x,x]
    , [x,x,x,x,x,x,x,x]
    , [x,x,x,x,x,r,x,x]
    , [r,x,r,x,x,x,x,x]
    , [x,x,x,x,x,r,x,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_rake ------------------- 
    ,[[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ])
    bra.append([
    #------------------- upper copy ------------------- 
    [[[D,r,U,x,D,x,D,x]
    , [x,x,x,x,r,x,r,x]
    , [U,r,D,x,U,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,x,D,r,U,x,D,r]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,r,S,x,W,r,S,x]
    , [x,x,x,x,x,x,x,x]
    , [W,x,S,x,S,r,W,x]
    , [r,x,r,x,x,x,x,x]
    , [S,x,W,x,W,r,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_rake ------------------- 
    ,[[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ])
    ### +---------------------------------------------------+
    ### |                   ket                             |
    ### +---------------------------------------------------+
    ket = []
    ket.append([
    #------------------- upper copy ------------------- 
    [[[U,x,D,x,U,x,D,x]
    , [b,x,b,x,b,x,b,x]
    , [D,x,U,x,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,x,S,x,W,x,S,x]
    , [b,x,b,x,b,x,b,x]
    , [S,x,W,x,S,x,W,x]
    , [x,x,x,x,x,x,x,x]
    , [W,b,S,x,W,b,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_make ------------------- 
    ,[[[g,g,g,g]
    , [g,g,g,g]
    , [g,g,g,g]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]
    ])

    
    ket.append(copy.deepcopy(ket[0]))
    ket[-1][2] = [
      [0,0,0,1]
    , [0,0,0,1]
    , [0,0,0,1]]
    
    ket.append(copy.deepcopy(ket[1]))
    ket[-1][0] = [
    #------------------- upper copy ------------------- 
     [[U,x,D,b,U,x,D,x]
    , [b,x,x,x,x,x,b,x]
    , [D,x,U,b,D,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,b,S,x,W,b,S,x]
    , [x,x,x,x,x,x,x,x]
    , [S,b,W,x,S,x,W,x]
    , [x,x,x,x,b,x,b,x]
    , [W,b,S,x,W,x,S,x]
    , [x,x,x,x,x,x,x,x]]]
    
    ket.append(copy.deepcopy(ket[1]))
    ket[-1][0] = [
    #------------------- upper copy ------------------- 
     [[x,x,x,b,x,x,x,x]
    , [b,x,x,x,x,x,b,x]
    , [x,x,x,b,x,x,x,x]
    , [x,x,x,x,x,x,x,x]
    , [x,b,x,x,x,b,x,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[x,b,x,x,x,b,x,x]
    , [x,x,x,x,x,x,x,x]
    , [x,b,x,x,x,x,x,x]
    , [x,x,x,x,b,x,b,x]
    , [x,b,x,x,x,x,x,x]
    , [x,x,x,x,x,x,x,x]]]
    ket[-1][2] = [
      [0,0,1,1]
    , [0,0,1,1]
    , [0,0,1,1]]
    
    ket.append([
    #------------------- upper copy ------------------- 
    [[[D,x,U,b,D,x,D,x]
    , [b,x,x,x,x,x,b,x]
    , [U,x,D,b,U,x,U,x]
    , [x,x,x,x,x,x,x,x]
    , [U,b,D,x,U,b,D,x]
    , [x,x,x,x,x,x,x,x]]
    #------------------- lower copy ------------------- 
    ,[[W,b,S,x,W,b,S,x]
    , [x,x,x,x,x,x,x,x]
    , [W,b,S,x,S,x,W,x]
    , [x,x,x,x,b,x,b,x]
    , [S,b,W,x,W,x,S,x]
    , [x,x,x,x,x,x,x,x]]]
    #------------------- upper site_make ------------------- 
    ,[[[g,g,g,g]
    , [g,g,g,g]
    , [g,g,g,g]]
    #------------------- lower site_make ------------------- 
    ,[[0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]]]
    #------------------- swap_map ------------------- 
    ,[[0,0,0,1]
    , [0,0,0,1]
    , [0,0,0,1]]
    ])
    
    tile_mark = []
    tile_mark.append([
      [0,0,0,0]
    , [0,0,0,0]
    , [0,0,0,0]
    ])
    #------------------- paint ------------------- 
    
    frame = []
    for b in [bra[0], ket[0]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[0], ket = ket[0], mod = mod, tile_mark = tile_mark[0]))
    
    for b in [ket[1]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[0], ket = ket[1], mod = mod, tile_mark = tile_mark[0]))
    
    for b in [bra[1]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[1], ket = ket[1], mod = mod, tile_mark = tile_mark[0]))
    
    for b in [ket[2]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[1], ket = ket[2], mod = mod, tile_mark = tile_mark[0]))
    
    for b in [bra[2], ket[3]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[2], ket = ket[3], mod = mod, tile_mark = tile_mark[0]))
    
    for b in [bra[3], ket[4]]:
        frame.append(generate_grid_2(swap = b[2], grid = b[0][0], grid_2 = b[0][1], site_mark = b[1][0], site_mark_2 = b[1][1], tile_mark = tile_mark[0], mod = mod))
    frame.append(generate_trans(bra = bra[3], ket = ket[4], mod = mod, tile_mark = tile_mark[0]))
    
    
    arw = group()
    arw.add(polyline(points = "{},{} {},{} {},{} {},{} {},{} {},{} {},{}".format(0, 10, 20, 10, 20, 0, 40, 15, 20, 30, 20, 20, 0, 20)))
    
    plus_i = group(stroke_width = 10, stroke = "black")
    plus_i.add(line(xy1 = (0, 15), xy2 = (30, 15)))
    plus_i.add(line(xy1 = (15, 0), xy2 = (15, 30)))
    
    pfeil = [arw, (40, 30)]
    pfeil2 = [copy.deepcopy(arw), (40, 30)]
    plus = [plus_i, (30, 30)]
    
    #------------------- init ------------------- 
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[0], plus, frame[1], pfeil, frame[2]]])
    pic.write_svg("svg/sqr_update_1.svg")
    
    #------------------- preswap ------------------- 
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[0], plus, frame[3], pfeil, frame[4]]])
    pic.write_svg("svg/sqr_update_2.svg")
    
    #------------------- update bonds ------------------- 
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[0], pfeil2, frame[5], pfeil, frame[6]]])
    pic.write_svg("svg/sqr_update_3.svg")
    
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[3], pfeil2, frame[7], pfeil, frame[8]]])
    pic.write_svg("svg/sqr_update_4.svg")
    
    #------------------- update spin ------------------- 
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[8], pfeil, frame[14]]])
    pic.write_svg("svg/sqr_update_5.svg")
    
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[14], pfeil, frame[12], plus, frame[13]]])
    pic.write_svg("svg/sqr_update_6.svg")
    
    #------------------- sample ------------------- 
    pic = canvas()
    pic.add(U_grad)
    pic.add(D_grad)
    pic.add(U2_grad)
    pic.add(D2_grad)
    pic.add([[frame[9], plus, frame[10], pfeil, frame[11]]])
    pic.write_svg("svg/sqr_update_7.svg")
    #~ pic.display()
    

if __name__ == "__main__":
    sqr_transition_graph()
    sqr_update_steps()
    
