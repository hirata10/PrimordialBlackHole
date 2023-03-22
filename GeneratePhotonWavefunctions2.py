#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve


from sympy import symbols
from sympy.physics.wigner import wigner_3j
import module1 as m1 # module1 contains the function to convert r_star to r\n",

import I_functions_class as Inp
from astropy.io import fits


# In[8]:


nu = 1.
k = 3.
h = 1.e-21
mu = 1.e-23 
lam = 1.
GC = 1.
M = 1.e21
c = 1.
tol = 1.e-10

r_initial = 2000.*M
r_final= -70.*M

#omega = 1.e-21
l = 2
#tol = 1.e-10

T = 1/(8*np.pi*M)
#omega = np.linspace(.01*T,10.*T,1000)
omega = np.linspace(.01*T,20.*T,2000)


# In[19]:

hdu = fits.Header()
hdu['type']='Photon Wavefunction'
hdu['r_init']=r_initial
hdu['r_final']=r_final
hdu['M']= M
hdu['l']=l
empty_primary = fits.PrimaryHDU(header=hdu)
hdul = fits.HDUList([empty_primary])

#hdul = fits.open('PhotonWaveFunctionFits/2.fits')

# In[20]:


for x in range(len(omega)):
    PWF = Inp.PhotonWaveFunction(M, omega[x], l, tol)
    r_gamma, F_points_gamma, z_points_gamma, f_points_gamma_prime = PWF.RK_4(r_initial, r_final, 240000, up = True)
    r_gamma_in, F_points_gamma_in, z_points_gamma_in, f_points_gamma_prime_in = PWF.RK_4(r_final, r_initial, 240000, up = False)
    Rgamma, Tgamma = PWF.get_R_and_T_coeff(r_gamma,F_points_gamma,z_points_gamma,r_gamma_in,F_points_gamma_in,z_points_gamma_in)
    
    #store every other point in the file
    #also save the reversed rin, and photon fields in so that you only have to store one r values
    
    #cutoff r to 2000M so that we can match analytically in the integrator later

    col1gammar = fits.Column(name='rpoints_up',format='D',array=r_gamma[::2])
    col2gammaF = fits.Column(name='F_points_up',format='M',array=F_points_gamma[::2])
    col2gammaz = fits.Column(name='z_points_up',format='M',array=z_points_gamma[::2])
    #col1gammarin = fits.Column(name='rpointsin',format='D',array=r_gamma_in)
    col2gammaFin = fits.Column(name='F_points_in',format='M',array=F_points_gamma_in[::-1][::2])
    col2gammazin = fits.Column(name='z_points_in',format='M',array=z_points_gamma_in[::-1][::2])
    
    cols = fits.ColDefs([col1gammar, col2gammaF,col2gammaz,col2gammaFin,col2gammazin])
    table_hdu = fits.BinTableHDU.from_columns(cols)
    table_hdu.header['omega']=omega[x]
    table_hdu.header['R']=Rgamma
    table_hdu.header['T']=Tgamma
    
    hdul.append(table_hdu)
    print(omega[x])


# In[21]:


hdul.writeto('/fs/scratch/PCON0003/emily/PhotonWaveFunctionFits/2ExtendedOmega.fits')





