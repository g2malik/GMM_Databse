# --- 
# aims: takes the frame with the larger number of UMZs as the frame of during and takes the frame with less as after
#This is different to creation which takes the frame with new as during rather than the one with less
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
    for i in range(len(array)):
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

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift\gaussian.txt'''
f = open (fname, mode = 'r')
spatial = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift\spatial.txt'''
s = open (spatial, mode = 'r')
quadrant = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\no_x_shift\quadrant.txt'''
q = open (quadrant, mode = 'r')
tol = 0.025
heights_dist = []
all_heights = []


UMZ_order = []
peaks_before = []
peaks_current = []
heights_current = []
heights_before = []

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

before =  np.zeros(2)

destruction = False
destruction_event = 0

for line in f:
    line_spat = s.readline()
    line_quad = q.readline()
    if (line.startswith("z") or line.startswith("t"))  == True:
        UMZ_order = []
        peaks_before = []
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

        """ 
        if destruction_event == 2:
            after_2_array.append(quadrant_current)
            after_array.append(quadrant_current)
            destruction_event = 0
        
        if destruction_event == 4:
            after_4_array.append(quadrant_current)
            after_array.append(quadrant_current)
            destruction_event = 0
        """

        destroyed_peaks = []
        destroyed_heights = []
        nearest_current_peaks = []
        if (len(peaks_before) - len(peaks_current)) == 1 and len(peaks_current)!=0:
            for j in range(len(peaks_before)):
                nearest_current_peaks.append(find_nearest(peaks_current, peaks_before[j])[1])
                if np.abs(peaks_before[j] - nearest_current_peaks[j])> tol: #doesnt include new peaks that are close to old ones
                    destroyed_peaks.append(peaks_before[j])
                    destroyed_heights.append(heights_before[j])
            #if len(destroyed_peaks)<1:
                #destroyed_peaks.append(np.abs(np.asarray(peaks_current) - np.asarray(nearest_current_peaks)).argmax()) #includes new close peaks
            
            if len(destroyed_peaks)==1: #If its a destruction event if all other are coherent
                destruction = True
                destruction_event = event_current

                before_array.append(before[0])
                during_array.append(before[1])
                after_array.append(quadrant_current)

                if destruction_event == 2:
                    before_2_array.append(before[0])
                    during_2_array.append(before[1])
                    after_2_array.append(quadrant_current)
                    
                if destruction_event == 4:
                    before_4_array.append(before[0])
                    during_4_array.append(before[1])
                    after_4_array.append(quadrant_current)
        
        peaks_before = peaks_current.copy()
        peaks_current = []
        heights_before = heights_current.copy()
        heights_current = []
        before = update_array(quadrant_current, before)


print(np.mean(before_array))
print(np.mean(during_array))
print(np.mean(after_array))

y_array = np.abs([np.mean(before_array), np.mean(during_array), np.mean(after_array)])
y_2_array = np.abs([np.mean(before_2_array), np.mean(during_2_array), np.mean(after_2_array)])
y_4_array = np.abs([np.mean(before_4_array), np.mean(during_4_array), np.mean(after_4_array)])

print(len(before_array))
print(len(during_array))
print(len(after_array))

plt.subplot(2,3,5)
plt.scatter(['before', 'during', 'after'], y_array, marker = '_', s = 8000)
plt.ylabel('Magnitude of average Quadrant event')

plt.subplot(2,3,1)
plt.scatter(['before', 'during', 'after'], y_2_array, marker = '_', s = 8000)
plt.ylabel('Magnitude of Q2 event')

plt.subplot(2,3,3)
plt.scatter(['before', 'during', 'after'], y_4_array, marker = '_', s = 8000)
plt.ylabel('Magnitude of Q4 event')


plt.show()
plt.close()
