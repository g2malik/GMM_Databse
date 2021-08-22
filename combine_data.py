# --- 
# aims: Combines all the data files into one file
# modefication history: gmalik, August, 2021; 

# --------------------------------
# import libraries 

import numpy as np
import time
import os
import sys
# --------------------------------

data_types = ['Gaussian', 'Quadrant', 'Spatial']
start_time = 0
no_files = 0
for type in data_types:
    path = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift'''
    wfile = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift\%s.txt'''%(type)
    w = open (wfile, mode = 'w')
    
    for x in range(1,4):
        for z in range(1,53):
            fname = path + '\%s\%s %d %d %d.dat'%(type, type, start_time, x, z)
            f = open (fname, mode = 'r')
            no_files+=1
            for line in f:
                w.write(line)
        print("Done")
    print(no_files)
