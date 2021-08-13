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

gaussian = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\gaussian.txt'''
f = open (gaussian, mode = 'r')
spatial = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\spatial.txt'''
s = open (spatial, mode = 'r')
quadrant = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\quadrant.txt'''
q = open (quadrant, mode = 'r')


tol = 0.025
dt = 0.023

all_vel = []
all_heights = []
all_spanwise = []
all_quadrant = []
all_event = []


UMZ_order = []
peaks_old = []
peaks_current = []
heights_current = []
spanwise_current = []
quadrant_current = 0
event_current = 0

counter = np.zeros((7,2))  #0: first recorded velocity / 1: # of frames
tracker = np.zeros(7)
hq_counter = np.zeros((7,4)) #0:Height 1:Spanwise 2:Quadrant creation events

coherence_dist = []

new_vel_dist = []
new_height_dist = []
new_spanwise_dist = []
new_quadrant_dist = []
new_event_dist = []

for line in f:
    line_spat = s.readline()
    line_quad = q.readline()
    if (line.startswith("z") or line.startswith("t"))  == True:
        UMZ_order = []
        peaks_old = []
        label = line
        z_coord = line.split()[2]
        z_coord = int(z_coord[2:])

        for jj in range(len(tracker)):
            if counter[jj,1] > 0:
                new_vel_dist.append(counter[jj, 0])
                coherence_dist.append(counter[jj,1])
                new_height_dist.append(hq_counter[jj, 0])
                new_spanwise_dist.append(hq_counter[jj,1])
                new_quadrant_dist.append(hq_counter[jj,2])
                new_event_dist.append(hq_counter[jj,3])
        counter[:,:] = 0
        tracker[:] = 0
        hq_counter[:,:] = 0

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
            all_vel.append(float(lst[ii]))
            heights_current.append(float(lst_spat[ii]))
            all_heights.append(float(lst_spat[ii]))
        for ii in range(3+UMZs_str, 3+(UMZs_str*2)):
            spanwise_current.append(float(lst_spat[ii]))
            all_spanwise.append(float(lst_spat[ii]))
        
        quadrant_current = float(lst_quad[3])
        all_quadrant.append(quadrant_current)
        event_current = lst_quad[5][1]
        if event_current == 'o':
            event_current = 0
        else:
            event_current = int(event_current)
        all_event.append(event_current)

        UMZ_order.append(int(UMZs_str))

        #Put code here to check if anything in counter and check if the velocity closes to the value in the tracker is
        #within tolerance. Then add 1 otherwise end counter and record number of frames

        for jj in range(len(tracker)):
            if tracker[jj] != 0: #If tracking at this place holder
                if np.abs(tracker[jj] - find_nearest(peaks_current, tracker[jj])[1]) < tol: #If peak in current coherent with tracker
                    counter[jj,1] += 1
                    tracker[jj] = find_nearest(peaks_current, tracker[jj])[1]
                else:
                    if counter[jj,1] > 0:
                        new_vel_dist.append(counter[jj,0])
                        coherence_dist.append(counter[jj,1])
                        new_height_dist.append(hq_counter[jj,0])
                        new_spanwise_dist.append(hq_counter[jj,1])
                        new_quadrant_dist.append(hq_counter[jj,2])
                        new_event_dist.append(hq_counter[jj,3])
                    counter[jj,:] = 0
                    tracker[jj] = 0
                    hq_counter[jj,:] = 0


        new_peaks = []
        new_heights = []
        new_spanwise = []
        nearest_old_peaks = []
        if (len(peaks_current) - len(peaks_old)) == 1 and len(peaks_old)!=0:
            for j in range(len(peaks_current)):
                nearest_old_peaks.append(find_nearest(peaks_old, peaks_current[j])[1])
                if np.abs(peaks_current[j] - nearest_old_peaks[j])> tol: #doesnt include new peaks that are close to old ones
                    new_peaks.append(peaks_current[j])
                    new_heights.append(heights_current[j])
                    new_spanwise.append(spanwise_current[j])
            """
            if len(new_peaks)<1:
                far_peak_index = np.abs(np.asarray(peaks_current) - np.asarray(nearest_old_peaks)).argmax()
                new_peaks.append(peaks_current[far_peak_index]) #includes new close peaks
                new_heights.append(heights_current[far_peak_index])
                new_spanwise.append(spanwise_current[far_peak_index])
                if len(new_peaks) != 1:
                        print("Close peaks error")
            """
            if len(new_peaks)==1: #If all other peaks are coherent
                new_space = find_counter_space()
                counter[new_space,0] = new_peaks[0] #To find average of velocity
                counter[new_space,1] = 1 #Start counting
                tracker[new_space] = new_peaks[0] #Update tracker
                hq_counter[new_space,0] = new_heights[0] #Init counter with first value of x parameter
                hq_counter[new_space,1] = new_spanwise[0]
                hq_counter[new_space,2] = quadrant_current
                hq_counter[new_space,3] = event_current

        
        peaks_old = peaks_current.copy()
        peaks_current = []
        heights_current = []
        spanwise_current = []
        quadrant_current = 0
        event_current = 0


bins_edge = np.linspace(0,1.1,12)
bar_edge = np.arange(0.5,1.05,0.1)

coherence_dist = np.array(coherence_dist) * dt 
print("The average coherence for new UMZs is ", np.mean(coherence_dist))


#plt.subplot(3, 5, 1)
coherence_hist = plt.hist(coherence_dist, bins = np.arange(1,17,1))
#plt.xlabel("Coherence of new UMZs")
plt.ylabel("Frequency")
print(coherence_hist[0])
 
plt.subplot(3, 4, 1)
all_vel_hist = plt.hist(all_vel)
#plt.xlabel("Streamwise velocity of new UMZs")
plt.ylabel("Total Frequency")

plt.subplot(3, 4, 2)
all_height_hist = plt.hist(all_heights, bins = bins_edge)
#plt.xlabel("Heights of new UMZs")


plt.subplot(3, 4, 3)
all_spanwise_hist = plt.hist(all_spanwise)
#plt.xlabel("Spanwise velocity of new UMZs")


plt.subplot(3, 4, 4)
all_quad_hist = plt.hist(all_quadrant)
#plt.xlabel("Quadrant event mag. of new UMZs")


#plt.subplot(3, 5, 5)
#all_event_hist = plt.hist(all_event)
#plt.xlabel("Quadrant event name")

plt.subplot(3, 4, 5)
new_vel_hist = plt.hist(new_vel_dist)
#plt.xlabel("Streamwise velocity of new UMZs")
plt.ylabel("Creation Frequency")

plt.subplot(3, 4, 6)
new_height_hist = plt.hist(new_height_dist, bins = bins_edge)
#plt.xlabel("Heights of new UMZs")

plt.subplot(3, 4, 7)
new_vel_hist = plt.hist(new_spanwise_dist)
#plt.xlabel("Spanwise velocity of new UMZs")

plt.subplot(3, 4, 8)
new_vel_hist = plt.hist(new_quadrant_dist)
#plt.xlabel("Quadrant event mag. of new UMZs")

#plt.subplot(3, 5, 10)
#new_vel_hist = plt.hist(new_event_dist)
#plt.xlabel("Quadrant event name of new UMZs")



plt.subplot(3, 4, 9)
plt.scatter(new_vel_dist, coherence_dist)
plt.xlabel(r'$U/U_{\infty}$')
plt.ylabel(r'$\delta / u_{\tau}$')

plt.subplot(3, 4, 10)
plt.scatter(new_height_dist, coherence_dist)
plt.xlabel(r'$y/\delta$')

plt.subplot(3, 4, 11)
plt.scatter(new_spanwise_dist, coherence_dist)
plt.xlabel(r'$V/V_{\infty}$')

plt.subplot(3, 4, 12)
plt.scatter(new_quadrant_dist, coherence_dist)
plt.xlabel("Quadrant event mag.")

#plt.subplot(3, 5, 15)
#plt.scatter(new_event_dist, coherence_dist)
#plt.xlabel("Quadrant event name")

plt.show()
plt.close()