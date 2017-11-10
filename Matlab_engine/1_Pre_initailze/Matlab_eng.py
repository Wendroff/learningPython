# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 10:16:43 2017

@author: lthpc
"""

import matlab.engine
eng = matlab.engine.start_matlab()


#import os
#path = r"H:\Work\learningPython\1_Pre_initailze"
#os.chdir(path)

eng.initial0(nargout=0)

#eng.quit()