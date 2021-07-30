# sub function called by read_reference 
# last modifed: gmalik, June, 2021 

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import h5py

from sklearn.mixture import GaussianMixture
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from scipy.signal import peak_prominences
#from matplotlib.pyplot import figure
 
class kde_plus_gmm():
    def __init__(self,XY,U_grid,V_grid,W_grid,U_hist,staz,binbin,jy,time_stamp,testing_file, delta_points):      
        
        data = np.reshape(U_hist,(1,-1))
        xx = np.linspace(np.min(data),np.max(data),500) #Evenly distributed velocity for the distribution curves only
        xx_interval = (np.max(data) - np.min(data))/500
        
                
        # ============ kde ============  
        bandwidth = (data.shape[1]**(-1/5))/1.65 #Slightly smaller than Scotts
        dist = int(0.045/xx_interval)
        kernel = gaussian_kde(data,bw_method=bandwidth) # Create an object with methods and atrributes of the distribution

        pks, props = find_peaks(kernel.evaluate(xx), distance = dist, height = 0.5, prominence = 0.01) #Evaluate Method gives the y values for given x values so can plot the x vs evaluted values to plot the function
        xpk = [xx[jj] for jj in pks] #Finds x-cordinates of the peaks. The indexes of the peaks were found above
        self.prominence = np.mean(props["prominences"])
                                
        
        nop = len(xpk) #Finds number of peaks and adds more if needed
        ypk = [kernel.evaluate(xx[jj]) for jj in pks] #Finds the the y value of the peaks
        xpk_sorted = -np.sort(-np.array(xpk)) # in decsending order   
        
        # ============ gmm ============  
        N = np.arange(1, nop+1) #An array of possible number of UMZs based on no.of peaks from KDE 
        models = [None for i in range(len(N))]
        X = np.array(U_hist).T.reshape(-1,1) #Reshapes the velocity array for GMM algorithm
        for i in range(len(N)):
            mini = np.reshape(xpk_sorted[:N[i]],(-1,1)) #column array of x-cordinates of first N[i] peaks for initilization
            wini = np.array([1/(i+1) for ik in range(i+1)]) # equal weights initially 
            models[i] = GaussianMixture( n_components = N[i], means_init = mini,  weights_init = wini, tol = 1e-4 ).fit(X) #Creates an array of model objects each with different no of peaks

        # compute the AIC and the BIC criterion 
        AIC = [m.aic(X) for m in models]
        BIC = [m.bic(X) for m in models]

        M_best = models[np.argmin(BIC)] # or AIC (Chooses the model with has the lowest BIC)
        
        U_wall = U_grid[:delta_points,:]
        W_wall = W_grid[:delta_points,:]
        Y_wall = XY[1][:delta_points,:]
        labels = M_best.predict(np.reshape(U_wall,(-1,1))) #Labels for the velcoities based on the best model
        labels = np.reshape(labels,U_wall.shape)
        
                
        self.N_best = N[np.argmin(BIC)]  ## number of kernels/gaussian   components in best model
       
        self.means_g = [] #Found means from GMM algorithm not the spatial UMZ
        self.std_g = [] #covariances or width of gaussians
        self.weights_g = list(M_best.weights_) #Total Probability of being in any of the UMZs
        
        self.spanwise = [] #Average Spanwise velocity of spatial UMZ rather than the Gaussian model
        self.height = []
       
        for ij in range(self.N_best):
            self.std_g.append(np.sqrt(M_best.covariances_[ij][0][0])) #Converting the resulting multidimensional array with empty dimensions to 1D
            self.means_g.append(M_best.means_[ij][0]) #The result is a multidimensional array but only has values in one of the dimensions

            self.spanwise.append(np.mean(W_wall[labels == ij]))
            self.height.append(np.mean(Y_wall[labels == ij]))
        
        self.means_g.reverse()
        self.std_g.reverse()
        self.weights_g.reverse()
        self.spanwise.reverse()
        self.height.reverse()
        

        
        #self.peaks_g = M_best.predict_proba(np.reshape(self.means_g, (-1, 1))) * np.exp(M_best.score_samples(np.reshape(self.means_g, (-1, 1))))[:, np.newaxis] #Does the same thing as finding pdf_individual in one line but for the mean values only
        self.peaks_g = M_best.predict_proba(np.reshape(self.means_g, (-1, 1)))
        self.peaks_g = self.peaks_g.max(axis = 1)
        
        """
        rsl = 500
        bins = binbin
        figure(figsize=(6.0, 4.0), dpi=rsl)
        plt.hist(data.T,bins,density=True, edgecolor='grey',facecolor='none')# PLots the velocity pdf
        plt.plot(xx,kernel.evaluate(xx),'k-',label='kernal density estimation(KDE)') #Plots the calculated pdf curve
        plt.scatter(xpk,ypk,c="r",s=80,marker="o",label='peaks of KDE')  #Plots the peaks from the KDE


        logprob = M_best.score_samples(xx.reshape(-1, 1))
        responsibilities = M_best.predict_proba(xx.reshape(-1, 1))
        pdf = np.exp(logprob)
        pdf_individual = responsibilities * pdf[:, np.newaxis] #The newaxis just makes the pdf it a column vector because you need to multiply it with the probabilities which are arranged in a column fashion

        
        #plt.plot(xx,pdf,'-',label='Gaussian Mixture Estimation')
        plt.plot(xx,pdf_individual,'--',label='individual Gaussian component')
        #plt.figtext(1,0,"std: %.2f prom: %.2f %s"%(np.std(props["peak_heights"]), np.mean(props["prominences"]), self.coherent), ha="center", fontsize=7)
        
       
        
        # ============ Graphing ============
        my_x_ticks = np.arange(0.4, 1.2, 0.1)
        my_y_ticks = np.arange(0, 13, 1)
        plt.xticks(my_x_ticks)
        plt.yticks(my_y_ticks)
    
        plt.xlabel("Streamwise Velocity",fontdict={'family' : 'Calibri', 'size':12})
        plt.ylabel("Frequncy",fontdict={'family' : 'Calibri', 'size':12})
        plt.title("restart 010%s p.d.f. X#%d at Zlabel = %d "%(time_stamp,jy,staz),fontdict={'family' : 'Calibri', 'size':12})
        plt.legend(loc='upper left', prop={'family':'Calibri', 'size':10},frameon=False)
        
        ax=plt.gca();# get the handle of the axis 
        ax.spines['bottom'].set_linewidth(0.5);
        ax.spines['left'].set_linewidth(0.5);
        ax.spines['right'].set_linewidth(0.5);
        ax.spines['top'].set_linewidth(0.5);
        
        
        #plt.savefig('/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/Q2_Detection/Results/010%s p.d.f x#%d Zlabel %d.png'%( time_stamp,jy,staz), facecolor='w')
        plt.show()
        plt.close()
        
        
        # ===== Contour snapshot of flow ======
        figure(figsize=(11.25, 4.5), dpi=1000)
        levs = np.linspace(0.55,1.1,250)
        thre =[np.arange(0.039,0.040,0.001)][0]
        stream = plt.contourf(XY[0][:delta_points,:],XY[1][:delta_points,:],U_wall,levs, cmap=plt.cm.jet)
        plt.contour(XY[0][:delta_points][:],XY[1][:delta_points][:delta_points],labels,2,colors='k', linewidths = 0.6, linestyles = 'dashed', alpha=0.5)

        #plt.plot([min(XY[0][0][:]), max(XY[0][0][:])], [0.17, 0.17])
        #isoline = plt.contour(XY[0][:delta_points,:],XY[1][:delta_points,:],TT[:delta_points,:],thre,colors=['white','white']) # plot the BTFTI 
        plt.clim(0, 1.03)
        cbar=plt.colorbar(stream, orientation='vertical', shrink=0.8)
        cbar.set_label('velocity U')
        
        plt.title("010%s stream snapshot X#%d at Zlabel = %d"%(time_stamp,jy,staz),fontdict={'family' : 'Calibri', 'size':12})
        #plt.savefig('/gpfs/fs0/scratch/j/jphickey/g2malik/working_code/temporal_coherence/Results/0 x#%d/010%s stream snapshot X#%d Zlabel %d.png'%(jy,time_stamp,jy,staz), facecolor='w')
        plt.show()
        plt.close()
        """
        
        """
        #Gives a graph for the profile with respect to x at a specific height. So can find a peak rather than the average
        figure(figsize=(6.0, 2.0), dpi=rsl)
        y_event_height = int(2)
        plt.plot(XY[0][y_event_height][:], U_grid[y_event_height,:])
        plt.xlabel("x")
        plt.ylabel("Streamwise Velocity")
        plt.show()
        plt.close()

        
        figure(figsize=(6.0, 2.0), dpi=rsl)
        plt.plot(XY[0][y_event_height][:], V_grid[y_event_height,:])
        plt.xlabel("x")
        plt.ylabel("Wall-normal Velocity")
        plt.show()
        plt.close()
        
        

        # ============ plot BIC/AIC ============
        figure(num=10, figsize=(7.0, 4.0), dpi=rsl)        
        plt.plot(AIC)
        plt.plot(BIC)
        plt.title("AIC/BIC plot")
        plt.show()
        plt.close()      
     
        """