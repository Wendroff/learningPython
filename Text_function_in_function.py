# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 20:00:46 2017

@author: lthpc
"""

#just a test about def-ing a function in anthor function


def mother_function(a):
    c = 1
    def test(a):
        return (a-c)

if __name__ == '__main__':
    a = 1
    print(mother_function(a))
    