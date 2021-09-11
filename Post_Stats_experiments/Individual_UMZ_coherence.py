# --- 
# aims: calculates how long all invidivual UMZs are coherent vs. their height
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

from itertools import count
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import figure 
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
    for index in range(np.shape(tracker[:,0])[0]):
        if tracker[index,0] == 0:
            return index

    print("no space in counter")
    
# --------------------------------
# main

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\band2\gaussian.txt'''
f = open (fname, mode = 'r')
tol = 0.02
dt = 0.023

coherent_velocities = []
coherent_times = []

peaks_old = []
peaks_current = []
nearest_old_peaks = []


counter = np.zeros((7,2))  #0: first recorded velocity / 1: # of frames
tracker = np.array([(0.0,False), (0.0,False), (0.0,False), (0.0,False), (0.0,False), (0.0,False), (0.0,False)]) #0: last recorded velocity / 1: if recently changed

labels = []
for line in f:
    #line = f.readline()
    if (line.startswith("z") or line.startswith("t"))  == True: # store label and counter and reset everything
        label  = line

        for jj in range(np.shape(counter[:,0])[0]): #Check the tracker to find peaks that werent updated
            if (counter[jj,1] > 5): #Adds data with atleast 3 subsequent frames
                #print("Adding to data")
                coherent_velocities.append(counter[jj,0])
                coherent_times.append(counter[jj,1])
        counter[: :] = 0
        tracker[:, :] = 0  
        
        peaks_old = []

    else:
        lst = line.split()
        UMZs_str = int(lst[0])
        prom = lst[1]
        for ii in range(3,3+UMZs_str):
            peaks_current.append(float(lst[ii]))
        nearest_old_peaks = []

        if len(peaks_old) != 0:
            for j in range(len(peaks_current)):
                
                nearest_old_peaks.append(find_nearest(peaks_old, peaks_current[j])[1])

                if np.abs(peaks_current[j] - nearest_old_peaks[j])< tol: #If it a coherent peak
                    count_index = find_tracker(nearest_old_peaks[j]) #Checks if already counting
                    if count_index != 7: #If counter is present
                        counter[count_index,1] +=1 #Add 1 to the location of counter
                        tracker[count_index,0] = peaks_current[j] #Put in current velocity in the tracker
                        tracker[count_index,1] = True #Says that the counter was updated
                    
                    elif count_index == 7: #Creates new counter and start counting
                        new_space = find_counter_space()
                        counter[new_space,0] = nearest_old_peaks[j]
                        #print(counter[new_space,0])
                        counter[new_space,1] = 2
                        tracker[new_space,0] = peaks_current[j]
                        tracker[new_space,1] = True
                
                else: #If not a coherent peak do nothing?
                    pass
            
            #print(counter)
            #print(coherent_times)
            for jj in range(np.shape(tracker[:,0])[0]): #Check the tracker to find peaks that werent updated
                if (tracker[jj,1] == False and tracker[jj,0] != 0):
                    if counter[jj,1] > 5: #Only add to data if present in atleast 3 frames
                        #print("Adding to data")
                        coherent_velocities.append(counter[jj,0])
                        coherent_times.append(counter[jj,1])
                    counter[jj,:] = 0
                    tracker[jj, 0] = 0
            
            tracker[:,1] = False
            
            #if len(new_peaks)<1:
                #new_peaks.append(np.abs(np.asarray(peaks_current) - np.asarray(nearest_old_peaks)).argmax()) #includes new close peaks
        peaks_old = peaks_current.copy()
        peaks_current = []

coherent_times = np.array(coherent_times) * dt

print(np.mean(coherent_velocities))
print(np.mean(coherent_times))

figure(dpi = 200)
#plt.scatter(coherent_velocities, coherent_times, s=30)
plt.hist2d(coherent_velocities, coherent_times, bins = 12, cmap = 'Reds')
plt.xlabel(r'$U/U_{\infty}$',fontdict={'family' : 'Calibri', 'size':12})
plt.ylabel(r'$\delta / u_{\tau}$',fontdict={'family' : 'Calibri', 'size':12})
#plt.title(" Temporal Coherence vs. Mean Velocity of UMZ ",fontdict={'family' : 'Calibri', 'size':12})
plt.show()