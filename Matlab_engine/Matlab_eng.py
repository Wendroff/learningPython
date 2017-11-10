# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 10:16:43 2017

@author: lthpc
"""

import matlab.engine
import os

#os.getcwd()
os.chdir('H:\\Work\\learningPython\\1_Pre_initailze')
#os.getcwd()

eng = matlab.engine.start_matlab()
eng.initial0(nargout=0)

eng.quit()