# --- 
# aims: implement the GMM and the KDE and Q2 to snapshots while changing the frame of reference
# calls: kde_plus_Q2
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl 
from scipy.interpolate import griddata 
from scipy import interpolate
#from pylab import *
from kde_plus_gmm import kde_plus_gmm
import sys
# --------------------------------

def main(argv):
    
    BUFFER = 105
    LOG = 140
    DELTA = 265

    # ------ Read data from file ------
    path = '/gpfs/fs0/scratch/j/jphickey/jphickey/Boundary_Layer_PNAS_2017/'
    testing_file = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/test.dat' # temporary file for debugging
    test = open(testing_file, "w")

    u_avg_file = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/Q2_Detection/u_profiles/final_avg_u.dat' # file with avg profile
    open_u = open(u_avg_file, "r")

    y_target = np.loadtxt(u_avg_file)[0]
    y_target = y_target[BUFFER-1:DELTA-1] # 104 because skip the buffer layer and until 275 which is approx 1 delta
    u_avg = np.loadtxt(u_avg_file)[1]
    u_avg = u_avg[BUFFER-1:DELTA-1]

    start_t = 0
    x_plane = int(argv[1])
    z_plane = int(argv[2])

    gaussian = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/Cluster Results/Gaussian/gaussian %d %d %d.dat'%(start_t, x_plane+1, z_plane)
    g = open(gaussian, "w") #Includes properties of gaussain curves

    spatial = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/Cluster Results/Spatial/spatial %d %d %d.dat'%(start_t, x_plane+1, z_plane)
    s = open(spatial, "w") #Includes height and spanwise velocity of spatial UMZ

    quadrant = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/Cluster Results/Quadrant/quadrant %d %d %d.dat'%(start_t, x_plane+1, z_plane)
    q = open(quadrant, "w") #Includes stuff about the quadrant event

    g.write("t:%d "%(start_t))
    g.write("x:%d "%(x_plane+1))
    g.write("z:%d "%(z_plane))
    g.write("\n")

    s.write("t:%d "%(start_t))
    s.write("x:%d "%(x_plane+1))
    s.write("z:%d "%(z_plane))
    s.write("\n")

    q.write("t:%d "%(start_t))
    q.write("x:%d "%(x_plane+1))
    q.write("z:%d "%(z_plane))
    q.write("\n")


        
    x_shift = 0 #Amount the whole BL moves in the next snapshot

    for snap in range (start_t,37):
        
        time_stamp = '%02d' %(snap)
        fname = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Code/Slices/restart_010' + time_stamp + '_ydelta_slices.dat'
        f = open (fname, mode = 'r') 

        # ----- choose the starting position/parameters -----
        btfti = 0.04 # see Wu, Wallace and Hickey 2019, PoF 
        binbin = 50 

        first_z = 1     # start from the position at Z-label = ? {1...513}
        last_z = 51
        first_x = 1     # start from the position at X-label = ? {0...2326}
        last_x = 2326
        first_y = BUFFER   # [1,400] skip the buffer layer -- 100 wall units 
        last_y = 400    

        y_interval = last_y - first_y + 1

        wall_length = 2000 
        x_interval = int(0.22 * wall_length) #0.29 ACCORDING TO SUMMARY BUT 440 POINTS IN CODE
        total_x_planes = int((last_x-first_x)/x_interval)
        # --------------------------------

        # ------ Calculate x_cordinates -----
        #z_plane = 1 # {1...51}
        #x_plane = 0 # {0...4}
        
        stax = first_x + x_shift + (x_interval * x_plane)
        skpx = last_x - x_interval 
        endx = stax + x_interval 
        # --------------------------------

        # ----skip / go to certain index -----
        #for i in range(3):
            #data = f.readline()
        for ii in range(stax-1):
            data = f.readline()
        for jj in range(first_y-1):
            for ii in range(last_x):
                data = f.readline()
        for kk in range(z_plane-1):
            for jj in range(last_y):
                for ii in range(last_x):
                    data = f.readline()
        # --------------------------------

        # ------- get velocity etc -------
        xy    = [[],[]] # xod and yod grid
        xx    = [] #Raw x
        yy    = [] #Raw y
        uu    = []      # computational grid
        vv    = []      #Actual y cordinate
        ww    = []
        delta = []
        tt    = []
        small = 1e-15

        for j in range(y_interval):
            ylb = first_y + j + 1 # y-label grid 
            for i in range(x_interval):
                data = f.readline()
                lst = data.split()
                x = float(lst[0]) - 10842.4  # minus the starting position at re_theta = 1800 
                y = float(lst[1])
                yod = float(lst[3])
                xod = x / (y+small) * yod
                if yod != 0:
                    d = y/yod #Calculates boundary layer thickness
                    delta.append(d)
                tl = float(lst[5]) # passive scalar, index of btfti 
                u = float(lst[7]) # Streamwise velocity
                v = float(lst[8]) # Wall-normal velocity
                w = float(lst[9]) # Spanwise velocity
                yy.append(y)
                xx.append(x)
                xy[0].append(xod)
                xy[1].append(yod)
                uu.append(u)
                vv.append(v)
                ww.append(w)
                tt.append(tl)
            for i in range(skpx): # skip the x that we don't need 
                data = f.readline()
        # --------------------------------

        # --------------------------------
        # interpolation to uniform grid    
        xmax = xy[0][-1]
        xmin = xy[0][0]
        ymax = xy[1][-1]
        ymin = xy[1][0]

        interpx = 400  # number of pts 
        interpy = 400
        xi=np.linspace(xmin,xmax,interpx)
        yi=np.linspace(ymin,ymax,interpy)

        XY  = np.meshgrid(xi,yi)
        UU  = griddata((xy[0],xy[1]), uu, (XY[0],XY[1]), method =  'cubic')
        VV  = griddata((xy[0],xy[1]), vv, (XY[0],XY[1]), method =  'cubic')
        WW  = griddata((xy[0],xy[1]), ww, (XY[0],XY[1]), method =  'cubic')
        TT  = griddata((xy[0],xy[1]), tt, (XY[0],XY[1]), method =  'cubic')
        if ymax ==0:
            break
        delta_points = int(interpy/ymax) #Finds the number of points that gives a yod value of 1
        #delta_points = DELTA-BUFFER #Largest number of delta_points

        # ----------------------------------
        # prepare the data for the histogram 
        uhis = []
        u_free = []

        for i in range(interpx):
            for j in range(interpy):
                # detection of the turbultent region; BTFTI 
                lll = TT[i][j]
                if lll > btfti :
                    uhis.append(UU[i][j]) #uhis doesnt include turbulent region but UU does
        
        # ----------------------------------
        # Calculate next frame
        x_next = 1.316070556640625 #x distance between adjascent indices
        dt = 112.5

        u_shift = np.mean(uhis)

        displacement = u_shift * dt

        x_shift = x_shift + int(displacement/x_next)
        
        # ----------------------------------
        # Process raw Data
        
        y_raw = np.reshape(yy,(y_interval,x_interval))
        y_delta_raw = np.reshape(xy[1],(y_interval,x_interval))
        x_delta_raw = np.reshape(xy[0],(y_interval,x_interval))
        u_raw = np.reshape(uu,(y_interval,x_interval))
        v_raw = np.reshape(vv,(y_interval,x_interval))
        
        # -----------------------------------
        # gmm - main
        
        model = kde_plus_gmm(XY,UU,VV,WW,uhis,z_plane,binbin,x_plane+1,time_stamp, delta_points)
        g.write("%d "%(model.N_best))
        g.write("%s "%(model.prominence))
        g.write("%s "%(x_delta_raw[0,0]))
        s.write("%d "%(model.N_best))
        s.write("%s "%(model.prominence))
        s.write("%s "%(x_delta_raw[0,0]))
        q.write("%d "%(model.N_best))
        q.write("%s "%(model.prominence))
        q.write("%s "%(x_delta_raw[0,0]))

        g.writelines(["%s " %mean for mean in model.means_g])
        g.writelines(["%s " %std for std in model.std_g])
        g.writelines(["%s " %weight for weight in model.weights_g])
        g.writelines(["%s " %peak for peak in model.peaks_g])
        g.write("\n")

        s.writelines(["%s " %height for height in model.height])
        s.writelines(["%s " %spanwise for spanwise in model.spanwise])
        s.write("\n")
        
        # --------------------------------
        # Quadrant Detection

        u_event_profile = np.mean(u_raw[:DELTA-BUFFER,:], axis = 1) #np.mean if averaging othwerise either amin or amax for largest Q2/Q4
        v_event_profile = np.mean(v_raw[:DELTA-BUFFER,:], axis = 1) #np.mean if averaging othwerise either amin or amax for largest Q2/Q4
        y_event_profile = y_raw[:DELTA-BUFFER,0]
        yod_event_profile = y_delta_raw[:DELTA-BUFFER,int(x_interval/2)] #Take the profile of the middle 

        u_diff_prof = u_event_profile - u_avg #If total profile is being subtracted then mean otherwise event
        v_diff_prof = v_event_profile
        y_diff_prof = y_event_profile - y_target

        uv_diff_prof = u_diff_prof * v_diff_prof

        i_min = np.argmin(uv_diff_prof ) #Index of the strongest Q2/Q4 event or least u'v' value
        uv_min = np.amin(uv_diff_prof) #u'v' value of the strongest event
        y_min = yod_event_profile[i_min] #y location of the strongest Q2/Q4


        if u_diff_prof[i_min] < 0 and uv_min < 0:
            main_event = 'Q2'
        elif u_diff_prof[i_min] > 0 and uv_min < 0:
            main_event = 'Q4'
        elif uv_min > 0:
            main_event = 'none'
            

        q.write("%s "%(uv_min))
        q.write("%s "%(y_min))
        q.write("%s "%(main_event))
        q.write("\n")


    test.close()
    g.close()
    s.close()
    q.close()


if __name__ == "__main__":
    main(sys.argv)