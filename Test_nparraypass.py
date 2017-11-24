# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:48:56 2017

@author: Wendroff

"""

import numpy as np

def test(b):
    bb = b[0:5]
    bb[0] = -1.0

a = np.random.rand(5,1)
print(a)


test(a)
print(a)