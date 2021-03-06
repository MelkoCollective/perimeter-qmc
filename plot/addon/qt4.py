#!/usr/bin/python3
# -*- coding: utf-8 -*-
# 
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    11.08.2014 19:52:02 CEST
# File:    qt4.py

from .helper import *

# this tries to import PySide ot PyQt4 
# (they are almost always interchangeable)
# both are python binding to the c++ QT library
# I slightly prefer pyside
qt_binding = "none"

try:
    from PySide.QtCore import *
    from PySide.QtGui import *
    qt_binding = "PySide"
    GREEN("PySide loaded")
except ImportError:
    try:
        from PyQt4.QtCore import *
        from PyQt4.QtGui import *
        qt_binding = "PyQt4"
        GREEN("PyQt4 loaded")
    except ImportError:
        ERROR("no PySide or PyQt4 module found")
