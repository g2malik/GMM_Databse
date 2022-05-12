# GMM_Databse
Summary
Uniform momentum zones(UMZs) are a fluid "structure" oberserbed in turbulent regimes. They are created by multiple side-by-side vortices that cause a region of a uniform
streamwise velocity. The purpose of my research was to use the data science concepts such as the gaussian mixture model (unsupervised learning) to detect these structures and
in 2TB of fluid simulation data. Then the detection algorithm could be used to gain new insight into the dynamics of UMZs and contribute to existing literature. The simulation
data is located on the Niagra computing cluster.  

kde_plus_gmm.py contains the code for the unsupervised learning algorithm  
read_data_slices.py / read_data_slicesSLURM.py contains the code to make the data from the simulation usable for the detection algorithm  
Post_Stats & Post_Stats_experiments contain the scripts for post processing the results and extracting statistics from the detected UMZs  


![FLuid Snapshot](https://user-images.githubusercontent.com/85200064/167990364-0f0b9617-7893-440e-83ba-9232c18eef1f.png)
![UMZ plot](https://user-images.githubusercontent.com/85200064/167990379-2b24e6f4-57bb-4aaf-bf98-8e2002d95800.png)

