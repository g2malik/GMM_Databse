# --- 
# aims: reads data on number of UMZs and gives a list of the number of frames a particular number of UMZs lasts
# calls: none
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

from types import ModuleType
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
from num2words import num2words
# --------------------------------

fname = r'''C:\Users\gagan\Documents\Work\Results\GMM Database\gaussian.txt'''
f = open (fname, mode = 'r')

dt = 0.023

one = []
two = []
three = []
four = []
five = []
six = [] #

most_coherent = []
UMZ_order = []

for line in f:
    if (line.startswith("z") or line.startswith("t"))  == True: #If heading analyse the array of # of UMZs
        time = 1
        for i in range(len(UMZ_order)):
            if i!=0 and UMZ_order[i] == UMZ_order[i-1]:
                time+=1
            if ((UMZ_order[i] != UMZ_order[i-1] or i==(len(UMZ_order)-1)) and time!=1):
                eval(num2words(UMZ_order[i-1])).append(time) #If not the same # anymore or last value, add value to array 
                if time>=7:
                    most_coherent.append(old_snap) #Identifies the locations where the number of UMZs stays constant for 7 frames
                time = 1
        old_snap = line
        UMZ_order = [] 
    else:
        lst = line.split()
        UMZs_str = lst[0]
        UMZ_order.append(int(UMZs_str)) #Build array of all the # of UMZs in this frame

#print(np.mean(two))
#print(np.mean(three))
#print(np.mean(four))
#print(np.mean(five))

two = np.array(two) * dt
three = np.array(three) * dt
four = np.array(four) * dt
five = np.array(five) * dt
six = np.array(six) * dt


plt.plot([2,3,4,5, 6], [np.mean(two), np.mean(three), np.mean(four), np.mean(five), np.mean(six)], color = 'k')
plt.xlabel("# of UMZs",fontdict={'family' : 'Calibri', 'size':12})
plt.xticks([2,3,4,5, 6])
plt.ylabel(r'$\delta / u_{\tau}$',fontdict={'family' : 'Calibri', 'size':12})
#plt.title(" Temporal Coherence vs. # of UMZ ",fontdict={'family' : 'Calibri', 'size':12})
plt.show()