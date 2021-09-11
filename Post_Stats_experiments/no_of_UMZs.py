# --- 
# aims: Gives the average number of UMZs in the data and its distribution
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

from types import ModuleType
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
# --------------------------------

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\gaussian.txt'''
f = open (fname, mode = 'r')

dt = 0.023

UMZ_order = []
six = 0

for line in f:
    if (line.startswith("z") or line.startswith("t"))  == True: #If heading analyse the array of # of UMZs
        pass

    else:
        lst = line.split()
        UMZs_str = int(lst[0])
        UMZ_order.append(int(UMZs_str)) #Build array of all the # of UMZs in this frame
        if UMZs_str == 1:
            six+=1



print(len(UMZ_order))
print(np.mean(UMZ_order))
UMZ_dist = plt.hist(UMZ_order, color = 'k', bins= 5)
plt.xlabel("# of UMZs",fontdict={'family' : 'Calibri', 'size':12})
print(UMZ_dist)
print(six)
plt.ylabel('Frequency')
plt.show()