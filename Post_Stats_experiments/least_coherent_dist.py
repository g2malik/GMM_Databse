# --- 
# aims: Investigate the height distribution of only the least coherent newly created UMZs
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
from num2words import num2words
from matplotlib.pyplot import figure

# --------------------------------
# find closest values
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx] #returns tuple of index and value

# --------------------------------
# find index of same value
def find_tracker(value):
    idx = np.where(tracker[:,0] == value)
    if np.shape(idx)[1] == 1:
        return idx[0]
    if np.shape(idx)[1] == 0:
        return 7 #7 means not present in tracker
    if np.shape(idx)[1] >1:
        print("tracker error")    

# --------------------------------
# find next available counter space
def find_counter_space():
    for index in range(len(tracker)):
        if tracker[index] == 0:
            return index

    print("no space in counter")
    
# --------------------------------
# main

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\gaussian.txt'''
f = open (fname, mode = 'r')
spatial = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\spatial.txt'''
s = open (spatial, mode = 'r')

tol = 0.025
dt = 0.023

heights_dist = []
all_heights = []


UMZ_order = []
peaks_old = []
peaks_current = []
heights_current = []

counter = np.zeros((7,3))  #0: first recorded velocity / 1: # of frames / 3: Height of creation
tracker = np.zeros(7)
coherence_dist = []
least_coherent_heights = []

for line in f:
    line2 = s.readline()
    if (line.startswith("z") or line.startswith("t"))  == True:
        UMZ_order = []
        peaks_old = []
        label = line
        z_coord = line.split()[2]
        z_coord = int(z_coord[2:])

        for jj in range(len(tracker)):
            if counter[jj,1] > 0:
                coherence_dist.append(counter[jj,1])
                if counter[jj,1] == 1:
                    least_coherent_heights.append(counter[jj,2])
        counter[:,:] = 0
        tracker[:] = 0

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

        #Put code here to check if anything in counter and check if the velocity closes to the value in the tracker is
        #within tolerance. Then add 1 otherwise end counter and record number of frames

        for jj in range(len(tracker)):
            if tracker[jj] != 0: #If tracking at this place holder
                if np.abs(tracker[jj] - find_nearest(peaks_current, tracker[jj])[1]) < tol: #If peak in current coherent with tracker
                    counter[jj,1] += 1
                    tracker[jj] = find_nearest(peaks_current, tracker[jj])[1] #Update tracker  
                else:
                    if counter[jj,1] > 0: #If tracking and no close peaks, stop counting and record total time
                        coherence_dist.append(counter[jj,1])
                        if counter[jj,1] == 1:
                            least_coherent_heights.append(counter[jj,2])
                            pass #Record height of UMZ

                    counter[jj,:] = 0
                    tracker[jj] = 0


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
                heights_dist.append(new_heights[0]) #Records the height of creation event
                new_space = find_counter_space()
                counter[new_space,2] = new_heights[0]
                counter[new_space,0] = new_peaks[0] #Init counter with first value
                counter[new_space,1] = 1 #Start counting
                tracker[new_space] = new_peaks[0] #Update tracker

        
        peaks_old = peaks_current.copy()
        peaks_current = []
        heights_current = []


bins_edge = np.linspace(0.1,1.1,12)
bar_edge = np.arange(0,1.05,0.1)

#coherence_dist = np.array(coherence_dist) * dt

print("The average coherence for new UMZs is ", np.mean(coherence_dist))

plt.subplot(2, 2, 1)
least_coherent_hist = plt.hist(least_coherent_heights, bins=bins_edge)
plt.xlabel("Heights of least coherent newly created UMZs")
plt.ylabel("Frequency")
print(len(least_coherent_heights))
least_frequens = least_coherent_hist[0]


plt.subplot(2, 2, 2)
new_hist = plt.hist(heights_dist, bins=bins_edge)
plt.xlabel("Heights newly created UMZs")
plt.ylabel("Frequency")
new_frequens = new_hist[0] #Gets the frequencies of the bins


plt.subplot(2, 2, 3)
all_hist = plt.hist(all_heights, bins=bins_edge)
plt.xlabel("Heights of all UMZs")
plt.ylabel("Frequency")
all_frequens = all_hist[0] #Gets the frequencies of the bins
plt.show()
plt.close()


""" Plot of the creation coherence
coherence_hist = plt.hist(coherence_dist, bins = np.arange(0.023,0.368,0.023))
plt.xlabel(r'$\delta / u_{\tau}$')
plt.ylabel("Frequency")
print(coherence_hist[0])
plt.show()
plt.close()

"""


#print("frames with wall creation: ", wall_labels)
#print("frames with middle1 creation: ", middle1_labels)