# --- 
# aims: reads data on number of UMZs and calcualtes where the new UMZs are created
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
from num2words import num2words

# --------------------------------
# find closest values
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx] #returns tuple of index and value

# --------------------------------
# main

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\gaussian.txt'''
f = open (fname, mode = 'r')
spatial = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\spatial.txt'''
s = open (spatial, mode = 'r')
tol = 0.025
heights_dist = []
all_heights = []
wall = 0
middle1 = 0
middle2 = 0
freestream = 0
random = 0

UMZ_order = []
peaks_old = []
peaks_current = []
heights_current = []

wall_labels = []
middle1_labels = []
z_labels = []
for line in f:
    line2 = s.readline()
    if (line.startswith("z") or line.startswith("t"))  == True:
        UMZ_order = []
        peaks_old = []
        label = line
        z_coord = line.split()[2]
        z_coord = int(z_coord[2:])


    else:
        lst = line.split()
        lst2 = line2.split()
        UMZs_str = int(lst[0])
        UMZs_str2 = int(lst2[0])
        if UMZs_str2 != UMZs_str:
            print("Wrong Allignment")
        std = lst[1]
        for ii in range(3,3+UMZs_str):
            peaks_current.append(float(lst[ii]))
            heights_current.append(float(lst2[ii]))
            all_heights.append(float(lst2[ii]))
        UMZ_order.append(int(UMZs_str))
        new_peaks = []
        new_heights = []
        nearest_old_peaks = []
        if (len(peaks_current) - len(peaks_old)) == 1 and len(peaks_old)!=0:
            for j in range(len(peaks_current)):
                nearest_old_peaks.append(find_nearest(peaks_old, peaks_current[j])[1])
                if np.abs(peaks_current[j] - nearest_old_peaks[j])> tol: #doesnt include new peaks that are close to old ones
                    new_peaks.append(j)
                    new_heights.append(heights_current[j])
            #if len(new_peaks)<1:
                #new_peaks.append(np.abs(np.asarray(peaks_current) - np.asarray(nearest_old_peaks)).argmax()) #includes new close peaks
            
            #print(new_peaks)
            if len(new_peaks)==1:
                heights_dist.append(new_heights[0])
                z_labels.append(z_coord)
                
                if new_peaks[0] == 0:
                    wall+=1
                    wall_labels.append(label)
                    wall_labels.append(UMZs_str)
                    wall_labels.append(std)
                elif new_peaks[0] == (len(peaks_current)-1):
                    freestream+=1
                elif new_peaks[0] == 1:
                    middle1+=1
                    middle1_labels.append(label)
                    middle1_labels.append(UMZs_str)
                    middle1_labels.append(std)
                elif new_peaks[0] == 2:
                    middle2+=1
            else:
                random+=1
        
        peaks_old = peaks_current.copy()
        peaks_current = []
        heights_current = []

print("# of UMZs created near wall: ", wall)
print("# of UMZs created at middle: ",middle1)
print("# of UMZs created at middle: ",middle2)
print("# of UMZs created near freastream: ",freestream)
print(random)

bins_edge = np.linspace(0,1.1,12)
bar_edge = np.arange(0.1,1.05,0.1)
#print(bins_edge)
#print(bar_edge)
test= np.ones(27)

plt.subplot(2, 2, 1)
new_hist = plt.hist(heights_dist, bins=bins_edge)
plt.xlabel("Heights of new UMZs")
plt.ylabel("Frequency")
new_frequens = new_hist[0] #Gets the frequencies of the bins
#print(new_frequens)

plt.subplot(2, 2, 2)
all_hist = plt.hist(all_heights, bins=bins_edge)
plt.xlabel("Heights of all UMZs")
all_frequens = all_hist[0] #Gets the frequencies of the bins
#print(all_frequens)

"""
percent_frequens = np.divide(new_frequens, all_frequens, out=np.zeros_like(new_frequens), where=all_frequens!=0) #Divides the frequencies except at /0
#print(percent_frequens)
plt.subplot(2, 2, 3)
plt.bar(bar_edge, percent_frequens, align='edge', width = 0.1 )
plt.xlabel("Heights of new UMZs")
plt.ylabel("'%' of all UMZs that are new")
width=(bins_edge[1] - bins_edge[0])
"""

plt.show()
plt.close()


#print("frames with wall creation: ", wall_labels)
#print("frames with middle1 creation: ", middle1_labels)