import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve
import sys 
import os


from sympy import symbols
from sympy.physics.wigner import wigner_3j
import module1 as m1 # module1 contains the function to convert r_star to r\n",

import I_functions_class as Inp
from astropy.io import fits

l = int(sys.argv[1])

hdu_c = fits.open('/users/PCON0003/koivuemily/PrimordialBlackHole/Constants.fits')


nu =  hdu_c[0].header['nu']
mu = hdu_c[0].header['mu']
lam = hdu_c[0].header['lam']
GC = hdu_c[0].header['GC']
c = hdu_c[0].header['c']
tol = hdu_c[0].header['tol']

M = hdu_c[0].header['M']
r_initial = hdu_c[0].header['r_init']
r_final= hdu_c[0].header['r_final']
T = hdu_c[0].header['Temp']

Photon_direc = hdu_c[0].header['P_direc']



omega = np.linspace(.01*T,20.*T,2000)



hdu = fits.Header()
hdu['type']='Photon Wavefunction'
hdu['r_init']=r_initial
hdu['r_final']=r_final
hdu['M']= M
hdu['l']=l
empty_primary = fits.PrimaryHDU(header=hdu)
hdul = fits.HDUList([empty_primary])



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





hdul.writeto(Photon_direc + str(l)+'ExtendedOmega.fits')
