# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:26:17 2017

@author: Wendroff
"""
import numpy as np

def DataCal(F3=0.5,G3=0.51):
    Row3 = np.array([0.47,0.45,0.455,0.44,0.48,F3,G3])
    Row4 = Aver_withoutM(Row3)
    Row5 = (Row3 - Row4)/Row4
    Row6 = []
    for C5 in Row5:
        Row6.append(Row6_fun(C5))
    #Row7 = np.array([1030820000.00,     1036000000.00,  1036000000.00, 	1036000000.00, 	1036000000.00, 	1025640000.00, 	1025640000.00])
    Row6 = np.array(Row6)
    Row8 = np.array([1036000000.00, 	1036000000.00, 	1036000000.00, 	1036000000.00, 	1036000000.00, 	1036000000.00, 	1036000000.00])
    Row7 = Row8.copy()
    Row7[0] *= 0.995
    Row7[5] *= 0.99
    Row7[6] *= 0.99
    
    Row9 = Row7.copy()
    Row10= Aver_withoutM(Row9)
    Row11= (Row9-Row10)/Row10
    Row12= []
    for C11 in Row11:
        Row12.append(Row12_fun(C11))
    Row14 = 1480000000.00*np.ones(7)
    Row13 = Row14*0.998
    Row15 = Row13*(1-Row3)
    Row16 = Row15.max()
    Row17 = 5.0*Row15/Row16
    Row18 = Row6 + Row12 + Row17
    Row19 = 0.5*np.ones(7)
    Row20 = Row18 * Row19
    print(Row20)
    H21 = Row20[0] - Row20[5] 
    I21 = Row20[0] - Row20[6] 
    return (H21,I21)

def Aver_withoutM(Row3):
    Row4 = Row3.sum()-Row3.max()-Row3.min()
    Row4 = Row4/5.0
    return Row4

def Row6_fun(Row5):
    if (Row5 < -0.05):
        return 77.5
    elif (Row5 < 0):
        return 80-np.floor(np.abs(Row5)*100)*0.5
    elif (Row5 < 0.05):
        return 79-np.floor(np.abs(Row5)*100)
    else:
        return 74

def Row12_fun(Row11):
    if (Row11 < -0.05):
        return 13.75
    elif (Row11 < 0):
        return 15-np.floor(np.abs(Row11)*100)*0.25
    elif (Row11 < 0.05):
        return 14.5-np.floor(np.abs(Row11)*100)*0.5
    else:
        return 12

if __name__ == '__main__':
    print(DataCal(0.4,0.51))