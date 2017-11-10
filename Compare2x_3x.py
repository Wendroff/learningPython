# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:09:04 2017

@author: lthpc
"""

class Foo(object):
    pass
class Foo1:
    pass
print(type(Foo), type(Foo1))
print(dir(Foo))
print(dir(Foo1))
print(isinstance(Foo, object))
print(isinstance(Foo1, object))