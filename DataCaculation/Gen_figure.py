# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:31:52 2017

@author: Wendroff
"""

import DataCal as dc
import pandas as pd


N   = 150
Big = -1e5
low_threshold = 0.4
high_threshold = 0.55
dx  = (high_threshold - low_threshold)/N
data=pd.DataFrame(columns=['F3','G3','H21','I21','DIF','DIF_wo_neg'])

output = open("tecout.dat","w")
output.write("VARIABLES = F3 G3 H21 I21 DIF DIF_without_negtive \n")
output.write("ZONE F = POINT I = " +  str(N+1) + '  J = ' +  str(N+1) + '\n')

F3___list  = []
G3___list  = []
H21__list  = []
I21__list  = []
DIF__list  = []
DIF2_list  = []

for i in range(N+1):
    for j in range(N+1):
        F3 = low_threshold + i*dx
        G3 = low_threshold + j*dx
        (H21, I21) = dc.DataCal(F3,G3)
        DIF = H21 - I21
        if ((H21>=0) & (I21>=0)):
            DIF2 = DIF
        else:
            DIF2 = Big
        
        output.write(" %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f \n" % (F3, G3, H21, I21, DIF, DIF2))
        F3___list.append(F3)
        G3___list.append(G3)
        H21__list.append(H21)
        I21__list.append(I21)
        DIF__list.append(DIF)
        DIF2_list.append(DIF2)
output.close()
pdout=pd.DataFrame({"F3":F3___list,"G3":G3___list,"H21":H21__list,"I21":I21__list,"DIF":DIF__list,"DIF_wo_neg":DIF2_list})
pdout.to_excel('pdout.xlsx')
