# --- 
# aims: takes the first frame with new UMZs as creation event rather than before the new UMZs is seen
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
from num2words import num2words
import copy 

# --------------------------------
# find closest values
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx] #returns tuple of index and value

# --------------------------------
# update array with new value
def update_array(value, array):
    for i in len(array):
        if i == (len(array) - 1):
            array[i] = value
        else:
            array[i] = array[i+1]
    return array
    
# --------------------------------
# find next available counter space
def find_counter_space():
    for index in range(len(tracker)):
        if tracker[index] == 0:
            return index

    print("no space in counter")
    
# --------------------------------
# main

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift_band2\gaussian.txt'''
f = open (fname, mode = 'r')
spatial = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift_band2\spatial.txt'''
s = open (spatial, mode = 'r')
quadrant = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift_band2\quadrant.txt'''
q = open (quadrant, mode = 'r')

tol = 0.025
dt = 0.023

heights_dist = []
all_heights = []


UMZ_order = []
peaks_old = []
peaks_current = []
heights_current = []
spanwise_current = []
quadrant_current = 0
event_current = 0

counter = np.zeros((7,2))  #0: first recorded velocity / 1: # of frames
tracker = np.zeros(7)
coherence_dist = []

before_array = []
during_array =  []
after_array = []

before_2_array = []
during_2_array =  []
after_2_array = []

before_4_array = []
during_4_array =  []
after_4_array = []

before =  0

creation = False
creation_event = 0

for line in f:
    line_spat = s.readline()
    line_quad = q.readline()
    if (line.startswith("z") or line.startswith("t"))  == True:
        UMZ_order = []
        peaks_old = []
        label = line
        z_coord = line.split()[2]
        z_coord = int(z_coord[2:])

    else:
        lst = line.split()
        lst_spat = line_spat.split()
        lst_quad = line_quad.split()
        UMZs_str = int(lst[0])
        UMZs_str2 = int(lst_spat[0])
        UMZs_str3 = int(lst_quad[0])
        if (UMZs_str and UMZs_str2) != UMZs_str3:
            print("Wrong Allignment")
        std = lst[1]
        for ii in range(3,3+UMZs_str):
            peaks_current.append(float(lst[ii]))
            heights_current.append(float(lst_spat[ii]))
            all_heights.append(float(lst_spat[ii]))
        for ii in range(3+UMZs_str, 3+(UMZs_str*2)):
            spanwise_current.append(float(lst_spat[ii]))
        
        quadrant_current = float(lst_quad[3])
        event_current = lst_quad[5][1]
        if event_current == 'o':
            event_current = 0
        else:
            event_current = int(event_current)

            
        if creation_event == 2:
            after_2_array.append(quadrant_current)
            after_array.append(quadrant_current)
            creation_event = 0
        
        if creation_event == 4:
            after_4_array.append(quadrant_current)
            after_array.append(quadrant_current)
            creation_event = 0

        new_peaks = []
        new_heights = []
        nearest_old_peaks = []
        if (len(peaks_current) - len(peaks_old)) == 1 and len(peaks_old)!=0:
            for j in range(len(peaks_current)):
                nearest_old_peaks.append(find_nearest(peaks_old, peaks_current[j])[1])
                if np.abs(peaks_current[j] - nearest_old_peaks[j])> tol: #doesnt include new peaks that are close to old ones
                    new_peaks.append(peaks_current[j])
                    new_heights.append(heights_current[j])
            #if len(new_peaks)<1:
                #new_peaks.append(np.abs(np.asarray(peaks_current) - np.asarray(nearest_old_peaks)).argmax()) #includes new close peaks
            
            if len(new_peaks)==1: #If its a creation event if all other are coherent
                creation = True
                creation_event = event_current
                if creation_event == 2:
                    before_2_array.append(before)
                    during_2_array.append(quadrant_current)
                    before_array.append(before)
                    during_array.append(quadrant_current)
                if creation_event == 4:
                    before_4_array.append(before)
                    during_4_array.append(quadrant_current)
                    before_array.append(before)
                    during_array.append(quadrant_current)
        
        peaks_old = peaks_current.copy()
        peaks_current = []
        heights_current = []
        before = copy.deepcopy(quadrant_current)

coherence_dist = np.array(coherence_dist) * dt

print(np.mean(before_array))
print(np.mean(during_array))
print(np.mean(after_array))

y_array = np.abs([np.mean(before_array), np.mean(during_array), np.mean(after_array)])
y_2_array = np.abs([np.mean(before_2_array), np.mean(during_2_array), np.mean(after_2_array)])
y_4_array = np.abs([np.mean(before_4_array), np.mean(during_4_array), np.mean(after_4_array)])

print(len(before_array))
print(len(during_array))
print(len(after_array))

#plt.subplot(2,3,5)
#plt.scatter(['before', 'during', 'after'], y_array, marker = '_', s = 8000)
#plt.ylabel('Magnitude of average Quadrant event')

plt.figure(dpi = 150)
#plt.subplot(2,3,1)
plt.scatter(['before', 'during', 'after'], y_2_array, marker = '_', s = 8000)
plt.ylabel('Magnitude of Q2 event')


#plt.figure(dpi = 150)
#plt.subplot(2,3,5)
plt.scatter(['before', 'during', 'after'], y_4_array, marker = '_', s = 8000)
plt.ylabel('Magnitude of Q4 event')


plt.show()
plt.close()
