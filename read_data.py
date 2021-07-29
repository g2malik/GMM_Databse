# --- 
# aims: implement the GMM and the KDE and Q2 to snapshots while changing the frame of reference
# calls: kde_plus_Q2
# modefication history: gmalik, July, 2021; 

# --------------------------------
# import libraries 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
import time 
from scipy.interpolate import griddata 
from scipy import interpolate
from pylab import *
from kde_plus_Q2 import kde_plus_Q2
# --------------------------------

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
y_target = y_target[BUFFER:DELTA] # 104 because skip the buffer layer and until 275 which is approx 1 delta
u_avg = np.loadtxt(u_avg_file)[1]
u_avg = u_avg[BUFFER:DELTA]

start_t = 0
xx = 0
for zz in range(1,513,10):  #Jinyuan used 10, I had 50 before
    gaussian = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/gaussian.dat' # temporary file for debugging
    g = open(gaussian, "w") #Includes properties of gaussain curves

    spatial = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/spatial.dat' # temporary file for debugging
    s = open(spatial, "w") #Includes height and spanwise velocity of spatial UMZ

    quadrant = '/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/GMM_Database/quadrant.dat' # temporary file for debugging
    q = open(quadrant, "w") #Includes stuff about the quadrant event

    g.write("t:%d "%(start_t))
    g.write("z:%d "%(zz))
    g.write("x:%d"%(xx+1))
    g.write("\n")

    s.write("t:%d "%(start_t))
    s.write("z:%d "%(zz))
    s.write("x:%d"%(xx+1))
    s.write("\n")

    q.write("t:%d "%(start_t))
    q.write("z:%d "%(zz))
    q.write("x:%d"%(xx+1))
    q.write("\n")

    u_event_profile = [] #peak streamwise velocity profile
    v_event_profile = [] #peak wall_normal velocity profile
    y_event_profile = [] #y_profile at each frame
    y_delta_profile = [] #y_profile at each frame

    y_main_event = [] # y location of the peak event in BL
    mag_main_event = [] #magnitude of the peak event in BL
    mag_log_event = [] #magnitude of the peak event in the log layer
    main_event = [] #Name of event in BL
    log_event = [] #Name of event in log layer
        
    x_shift = 0 #Amount the whole BL moves in the next snapshot

    for time in range (start_t,37):
        time_stamp = '%02d' %(time)
        fname = path + 'restart_010' + time_stamp + '_ydelta_adrian_scalar_omega_uvw_08240_10565.dat'
        f = open (fname, mode = 'r') 
        print('--- the %d th snapshot in time ---'%time)

        # ----- choose the starting position/parameters -----
        btfti = 0.04 # see Wu, Wallace and Hickey 2019, PoF 
        binbin = 50 

        first_z = 1     # start from the position at Z-label = ? {1...513}
        last_z = 513
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
        z_plane = zz # {1...513}
        x_plane = xx # {0...4}

        stax = first_x + x_shift + (x_interval * x_plane)
        skpx = last_x - x_interval 
        endx = stax + x_interval 
        # --------------------------------

        # ----skip / go to certain index -----
        for i in range(3):
            data = f.readline()
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
        xy    = [[],[]] # computational grid
        xx    = []
        yy    = []
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
                xy[0].append(x)
                xy[1].append(y)
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

        # Calculate next frame)
        x_next = 1.316070556640625 #x distance between adjascent indices
        dt = 112.5

        u_shift = np.mean(uhis)

        displacement = u_shift * dt

        x_shift = x_shift + int(displacement/x_next)
        
        # Quadrant Detection
        # --------------------------------
        x_raw = np.reshape(xx,(y_interval,x_interval))
        y_raw = np.reshape(yy,(y_interval,x_interval))
        y_delta_raw = np.reshape(xy[1],(y_interval,x_interval))
        u_raw = np.reshape(uu,(y_interval,x_interval))
        v_raw = np.reshape(vv,(y_interval,x_interval))
        u_event_profile.append(np.mean(u_raw[:DELTA-BUFFER,:], axis = 1))
        v_event_profile.append(np.mean(v_raw[:DELTA-BUFFER,:], axis = 1))
        y_event_profile.append(y_raw[:DELTA-BUFFER,0])
        y_delta_profile = np.mean(y_delta_raw[:DELTA-BUFFER,:], axis = 1)      
        
        # --------------------------------
        # Q2/Q4 - main
        model = kde_plus_Q2(XY,UU,VV,WW,uhis,z_plane,binbin,x_plane+1,time_stamp,testing_file, delta_points)
        g.write("%d "%(model.N_best))
        g.write("%.7f"%(model.prominence))
        g.write("%s"%(xx[0,0]))
        s.write("%d "%(model.N_best))
        s.write("%.7f"%(model.prominence))
        s.write("%s"%(xx[0,0]))
        q.write("%d "%(model.N_best))
        q.write("%.7f"%(model.prominence))
        q.write("%s"%(xx[0,0]))

        print("Starting x is ",xx[0,0])
        print("Displacement is ",displacement)

        g.writelines(["%s " %mean for mean in model.means_g])
        g.writelines(["%s " %std for std in model.std_g])
        g.writelines(["%s " %weight for weight in model.weights_g])
        g.writelines(["%s " %peak for peak in model.peaks_g])
        g.write("\n")

        s.writelines(["%s " %height for height in model.height])
        s.writelines(["%s " %spanwise for spanwise in model.spanwise])
        s.write("\n")
        
        # --------------------------------
          

    u_diff_prof = np.array(u_event_profile) - u_avg #If total profile is being subtracted then mean otherwise event
    v_diff_prof = np.array(v_event_profile)
    y_diff_prof = np.array(y_event_profile) - y_target
    print('Mean of y_profile difference is ', np.mean(y_diff_prof))
    print('Std of y profile is ', np.std(y_diff_prof))

    log_layer_i = LOG-BUFFER #Number of indices until the end of the log layer so this is causing the problems

    uv_diff_prof = u_diff_prof * v_diff_prof

    i_min = np.argmin(uv_diff_prof, axis = 1 ) #Index of the strongest Q2/Q4 event or least u'v' value
    uv_min = np.amin(uv_diff_prof, axis = 1 ) #u'v' value of the strongest event
    y_min = y_delta_profile[i_min] #y location of the strongest Q2/Q4


    i_log_min = np.argmin(uv_diff_prof[:,:log_layer_i], axis = 1) #Index of the strongest event in the log layer
    uv_log_min = np.amin(uv_diff_prof[:,:log_layer_i], axis = 1 ) #u'v' value of the strongest event
    y_log_min = y_delta_profile[i_log_min]


    for ii in range(len(u_event_profile)):
        # See if each frame has a dominant Q2 or Q4 event in total and in the log region
        if u_diff_prof[ii,i_min[ii]] < 0 and uv_min[ii] < 0:
            main_event.append('Q2')
        elif u_diff_prof[ii, i_min[ii]] > 0 and uv_min[ii] < 0:
            main_event.append('Q4')
        elif uv_min[ii] > 0:
            main_event.append('none')
            
        if u_diff_prof[ii,i_log_min[ii]] < 0 and uv_log_min[ii] < 0:
            log_event.append('Q2')
        elif u_diff_prof[ii, i_log_min[ii]] > 0 and uv_log_min[ii] < 0:
            log_event.append('Q4')
        elif uv_log_min[ii] > 0:
            log_event.append('none')


        q.write("%s "%(uv_min[ii]))
        q.write("%s "%(y_min[ii]))
        q.write("%s "%(main_event[ii]))

        q.write("%s "%(uv_log_min[ii]))
        q.write("%s "%(y_log_min[ii]))
        q.write("%s "%(log_event[ii]))
        q.write("\n")

        """
        print("Magnitude of event is ", uv_min[ii])
        print("Postion of event is ", y_min[ii])
        print("Event name is ", main_event[ii])
        
        print("Magnitude of log event is ", uv_log_min[ii])
        print("Postion of log event is ", y_log_min[ii])
        print("Log event name is ", log_event[ii])
        
        plt.plot(u_diff_prof[ii,:], y_delta_profile)
        plt.xlabel("Streamwise Velocity")
        plt.ylabel("y")
        #plt.xlim([-0.045,0.045])
        plt.title("010%02d mean U X#%d at Zlabel = %d"%(ii,x_plane+1,z_plane))
        
        #plt.savefig('/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/Q2_Detection/Results/mean/010%02d mean U x#%d Zlabel %d.png'%(ii,x_plane+1,z_plane), facecolor='w')
        plt.show()
        plt.close()
        
        plt.plot(v_diff_prof[ii,:], y_delta_profile)
        plt.xlabel("Wall-Normal Velocity")
        plt.ylabel("y")
        #plt.xlim([-0.02, 0.02])
        plt.title("010%02d mean V X#%d at Zlabel = %d"%(ii,x_plane+1,z_plane))
        #plt.savefig('/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/Q2_Detection/Results/mean/010%02d mean V x#%d Zlabel %d.png'%(ii,x_plane+1,z_plane), facecolor='w')
        plt.show()
        plt.close()
        
        plt.plot(uv_diff_prof[ii,:], y_delta_profile)
        plt.xlabel("u'v'")
        plt.ylabel("y")
        #plt.xlim([-0.02, 0.02])
        plt.title("010%02d U'V' X#%d at Zlabel = %d"%(ii,x_plane+1,z_plane))
        #plt.savefig('/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/Q2_Detection/Results/mean/010%02d mean V x#%d Zlabel %d.png'%(ii,x_plane+1,z_plane), facecolor='w')
        plt.show()
        plt.close()
        """




test.close()
print ('--- End of all ---')