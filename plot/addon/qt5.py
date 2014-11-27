#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    11.08.2014 19:53:26 CEST
# File:    qt5.py

from .helper import *

qt_binding = "none"

# this tries to import PyQt5

try:
    from PyQt5.QtCore import *
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    qt_binding = "PyQt5"
    GREEN("PyQt5 loaded")
except ImportError:
    ERROR("no PyQt5 module found")
