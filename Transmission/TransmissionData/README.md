This folder contains the files with the transmission coefficient data to be used to determine b0 and b1 terms for the
stochastic charge project, besides for the principal part integral in b0 (which is done with another code).\n
\n
It contains data for the following masses:\n
1e21 Planck masses\n
2e21 Planck masses\n
4e21 Planck masses\n
8e21 Planck masses\n
\n
The parameters are as follows:\n
k=[-5,-4,-3,-2,-1,1,2,3,4,5]\n
Z=[-2,-1,0,1,2] (for Y=0)\n
Y=[-2,-1,0,1,2] (for Z=0)\n
rmax=10^6 * BH mass\n
nstep=1.81*10^6\n
e-folds=26.6\n
\n
For "non-fine" files:\n
dh=0.1*TH (where TH is the Hawking temperature for each mass)\n
hnum=200\n
hmin=0.01*TH\n
\n
For "fine" files:\n
dh=0.01*TH\n
hnum=200\n
For M=1e21: hmin=0.01*TH+10*0.1*TH\n
For M=2e21: hmin=0.01*TH+20*0.1*TH\n
For M=4e21: hmin=0.01*TH+41*0.1*TH\n
For M=8e21: hmin=0.01*TH+83*0.1*TH
