#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    27.08.2014 15:46:06 CEST
# File:    svg.py

from .helper import *
from .parameter import bash
import numpy as np

def unpack_lazy(data):
    """
    Enables some shorthands for the interface.
    """
    if "dir" in data:
        if data["dir"] in ["u", "up"]:
            data["xy1"] = (0, 1)
            data["xy2"] = (0, 0)
        if data["dir"] in ["d", "down"]:
            data["xy1"] = (0, 0)
            data["xy2"] = (0, 1)
        if data["dir"] in ["r", "right"]:
            data["xy1"] = (0, 0)
            data["xy2"] = (1, 0)
        if data["dir"] in ["l", "left"]:
            data["xy1"] = (1, 0)
            data["xy2"] = (0, 0)
        
        if data["dir"] in ["dl", "ld"]:
            data["xy1"] = (1, 0)
            data["xy2"] = (0, 1)
        if data["dir"] in ["dr", "rd"]:
            data["xy1"] = (0, 0)
            data["xy2"] = (1, 1)
        if data["dir"] in ["ul", "lu"]:
            data["xy1"] = (1, 1)
            data["xy2"] = (0, 0)
        if data["dir"] in ["ur", "ru"]:
            data["xy1"] = (0, 1)
            data["xy2"] = (1, 0)
        
        del data["dir"]
        
    if "xy" in data:
        data["x"] = data["xy"][0]
        data["y"] = data["xy"][1]
        del data["xy"]
        
    if "cxy" in data:
        data["cx"] = data["cxy"][0]
        data["cy"] = data["cxy"][1]
        del data["cxy"]
        
    if "rxy" in data:
        data["rx"] = data["rxy"][0]
        data["ry"] = data["rxy"][1]
        del data["rxy"]
        
    if "xy1" in data:
        data["x1"] = data["xy1"][0]
        data["y1"] = data["xy1"][1]
        del data["xy1"]
        
    if "xy2" in data:
        data["x2"] = data["xy2"][0]
        data["y2"] = data["xy2"][1]
        del data["xy2"]
        
    if "size" in data:
        data["width"] = data["size"][0]
        data["height"] = data["size"][1]
        del data["size"]
    
    return data

ref_name = []
"""
All objects that can be refered to are in this list, since url(#REF) needs to br wrapped arount these calls.
"""
def parse_dict(data, no_print, only_print = None):
    """
    Prints key, val as key = "val" and caches special cases. All keys in no_print are ignored and all keys in only_print are printed if not in no_print.
    """
    
    def special_case(key, val):
        if is_number(val):
            val = "{:.3f}".format(val) #for nicer svg code
        
        if val in ref_name and key != "id":
            val =  "url(#" + val + ")"
        
        return key.replace("_", "-") + ' = "' + str(val) + '" '
    
    var = ""
    for k, v in data.items():
        if k not in no_print:
            if only_print:
                if k in only_print:
                    var += special_case(k, v)
            else:
                var += special_case(k, v)
    return var

class rect():
    """
    x, y, width, height, rx, ry # rx, ry stand for rounded edges
    http://www.w3.org/TR/SVG/shapes.html#RectElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <rect "
        var += parse_dict(self.data, self.internal)
        var += "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class circle():
    """
    cx, cy, r
    http://www.w3.org/TR/SVG/shapes.html#CircleElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <circle "
        var += parse_dict(self.data, self.internal)
        var += "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class ellipse():
    """
    cx, cy, rx, ry
    http://www.w3.org/TR/SVG/shapes.html#CircleElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <ellipse "
        var += parse_dict(self.data, self.internal)
        var += "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class line():
    """
    x1, y1, x2, y2, stroke_linecap = butt | round | square, stroke_dasharray = "sequence paint, no paint: 5, 5, 10, 5"
    http://www.w3.org/TR/SVG/shapes.html#LineElement
    http://www.w3schools.com/svg/svg_stroking.asp
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <line "
        var += parse_dict(self.data, self.internal)
        var += "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class polyline():
    """
    points = "50,375 150,375 150,325", stroke_linejoin = miter | round | bevel
    http://www.w3.org/TR/SVG/shapes.html#PolylineElement
    https://mdn.mozillademos.org/files/731/SVG_Stroke_Linejoin_Example.png
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <polyline " + parse_dict(self.data, self.internal) + "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class polygon():
    """
    points = "50,375 150,375 150,325", stroke_linejoin = miter | round | bevel
    www.w3.org/TR/SVG/shapes.html#PolygonElement
    https://mdn.mozillademos.org/files/731/SVG_Stroke_Linejoin_Example.png
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <polygon " + parse_dict(self.data, self.internal) + "/>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class text():
    """
    x, y, dx, dy (offset), text
    http://www.w3.org/TR/SVG/text.html#TextElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type", "text"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "    <text " + parse_dict(self.data, self.internal) + ">\n"
        var += "      " + self.data["text"] + "\n"
        var += "    </text>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class lingrad():
    """
    id, x1, y1, x2, y2, color1, color2, opacity1, opacity2, gradientUnits = userSpaceOnUse | objectBoundingBox, dir = l | r | u | d | dr | dl | ur | ul | rd | ru | lu | ld 
    http://www.w3.org/TR/SVG/pservers.html#LinearGradients
    http://www.w3schools.com/svg/svg_grad_linear.asp
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        global ref_name
        ref_name.append(kwargs["id"])
        
        self.default()
        self.internal = ["type"]
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "def"
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var  = "    <linearGradient " + parse_dict(self.data, self.internal, ["x1", "y1", "x2", "y2", "id"]) + ">\n"
        var += '      <stop offset = "{}" stop-color = "{}" stop-opacity = "{}"/>\n'.format(0, self["color1"], self.data.get("opacity1", 1))
        var += '      <stop offset = "{}" stop-color = "{}" stop-opacity = "{}"/>\n'.format(1, self["color2"], self.data.get("opacity2", 1))
        var += "    </linearGradient>\n"
        
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class group():
    """
    All properies defined in the group are carried over to all items in that group if not overwritten.
    http://www.w3.org/TR/SVG/struct.html#GElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = ["type"]
        self.defs = []
        self.items = []
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs:
            self[key] = kwargs[key]
    
    def default(self):
        self.data = {}
        self.data["type"] = "item"
    
    def add(self, svg_item):
        if svg_item["type"] == "def":
            self.defs.append(svg_item)
        elif svg_item["type"] == "item":
            self.items.append(svg_item)
    
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "  <g " + parse_dict(self.data, self.internal) + ">\n"
        
        if len(self.defs) != 0:
            var += "  <defs>\n"
            for d in self.defs:
                var += d.svg_parse()
            var += "  </defs>\n"
        
        for item in self.items:
            var += item.svg_parse()
            
        var += "  </g>\n"
        return var
        
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        
    def __getitem__(self, key):
        return self.data[key]

class canvas():
    """
    The outer most element.
    http://www.w3.org/TR/SVG/struct.html#SVGElement
    """
    #------------------- ctor ------------------- 
    def __init__(self, **kwargs):
        self.default()
        self.internal = []
        
        kwargs = unpack_lazy(kwargs)
        for key in kwargs.keys():
            self[key] = kwargs[key]
        
        self.items = []
        self.group = []
        self.defs = []
    
    def default(self):
        self.data = {}
        self["width"] = None
        self["height"] = None
    
    def add(self, svg_item):
        if is_list(svg_item):
            self.add_frame(svg_item)
        elif svg_item["type"] == "def":
            self.defs.append(svg_item)
        elif svg_item["type"] == "item":
            self.items.append(svg_item)
    
    def add_frame(self, matrix):
        """
        Special function that takes a 2D list with (svg_item, size) elements where size = (width, height).
        translates all svg_items in the matrix, s.t. all are centered in the minimal necessary space.
        """
        size = [[m[1] for m in mm] for mm in matrix]
        svg_item = [[m[0] for m in mm] for mm in matrix]
        max_y = [max([x[1] for x in mm]) for mm in size]
        max_x = [max([x[0] for x in mm]) for mm in [list(i) for i in zip(*size)]]
        
        acc_x = [0] + list(np.add.accumulate(max_x))
        acc_y = [0] + list(np.add.accumulate(max_y))
        
        size2 = [[(max_x[i], max_y[j]) for i in range(len(size[0]))] for j in range(len(size))]
        pos = [[(acc_x[i], acc_y[j]) for i in range(len(size[0]))] for j in range(len(size))]
        
        for svg_row, i in zipi(svg_item):
            for svg, j in zipi(svg_row):
                p = pos[i][j]
                s = size[i][j]
                s2 = size2[i][j]
                x = p[0] + (s2[0] - s[0]) / 2.0
                y = p[1] + (s2[1] - s[1]) / 2.0
                
                svg_item[i][j]["transform"] = "translate({}, {})".format(x, y)
                self.add(svg_item[i][j])
        
        self["size"] = (acc_x[-1], acc_y[-1])
        
    #------------------- const ------------------- 
    def __str__(self):
        return str(self.data)
    
    def svg_parse(self):
        var = "<?xml version = \"1.0\"?>\n<svg "
        var += parse_dict(self.data, self.internal) +">\n"
        
        if len(self.defs) != 0:
            var += "  <defs>\n"
            for d in self.defs:
                var += d.svg_parse()
            var += "  </defs>\n"
        
        for item in self.items:
            var += item.svg_parse()
        
        var += "</svg>\n"
        return "".join(var)
    
    def write_svg(self, filename = None):
        self.filename = filename
        os = open(self.filename,'w')
        os.write(self.svg_parse())
        os.close()
    
    def display(self):
        bash("{} {}".format("display", self.filename))
    #------------------- getter/setter ------------------- 
    def __setitem__(self, key, val):
        self.data[key] = val
        unpack_lazy(self.data)
        
    def __getitem__(self, key):
        return self.data[key]
