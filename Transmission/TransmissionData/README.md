This folder contains the files with the transmission coefficient data to be used to determine b0 and b1 terms for the
stochastic charge project, besides for the principal part integral in b0 (which is done with another code).  
  
It contains data for the following masses:  
1e21 Planck masses  
2e21 Planck masses  
4e21 Planck masses  
8e21 Planck masses  
  
The parameters are as follows:  
k=[-5,-4,-3,-2,-1,1,2,3,4,5]  
Z=[-2,-1,0,1,2] (for Y=0)  
Y=[-2,-1,0,1,2] (for Z=0)  
rmax=10^6 * BH mass  
nstep=1.81*10^6  
e-folds=26.6  
  
For "non-fine" files:  
dh=0.1*TH (where TH is the Hawking temperature for each mass)  
hnum=200  
hmin=0.01*TH  
  
For "fine" files:  
dh=0.01*TH  
hnum=200  
For M=1e21: hmin=0.01*TH+10*0.1*TH  
For M=2e21: hmin=0.01*TH+20*0.1*TH  
For M=4e21: hmin=0.01*TH+41*0.1*TH  
For M=8e21: hmin=0.01*TH+83*0.1*TH
