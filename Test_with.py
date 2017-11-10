# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:48:50 2017

@author: lthpc
"""

class Sample:
    def __enter__(self):
        print("In __enter__()")
        return "Foo"
 
    def __exit__(self, type, value, trace):
        print ("In __exit__()")

def get_sample():
    print("In def get_sample")
    return Sample()
 
 
with get_sample() as sample:
    print ("sample:",sample)