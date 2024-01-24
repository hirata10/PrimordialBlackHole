import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve

from sympy import symbols
from sympy.physics.wigner import wigner_3j
import module1 as m1 # module1 contains the function to convert r_star to r\n",

import cmath
from astropy.io import fits
import fitsio
import I_functions_class as Inp
from importlib import reload
import sys
reload(Inp)

x = int(sys.argv[1])
pathe = sys.argv[2]
patho = sys.argv[3]
omega_scale = float(sys.argv[4])
l= int(sys.argv[5])

hdu_c = fits.open('/users/PCON0003/koivuemily/PrimordialBlackHole/ConstantsM2.fits')
nu =  hdu_c[0].header['nu']
mu = hdu_c[0].header['mu']
lam = hdu_c[0].header['lam']
GC = hdu_c[0].header['GC']
c = hdu_c[0].header['c']
tol = hdu_c[0].header['tol']
direcPhoton = hdu_c[0].header['P_direc']
direcElectron = hdu_c[0].header['E_direc']

M = hdu_c[0].header['M']
T = hdu_c[0].header['Temp']
omega=omega_scale*T

X = 0
X_prime = 0
X_gamma = 1
parity = 1
n=1
lim1=0
lim2=1

ks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]
k_primes=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]
omega_index = round(omega*100*(8*np.pi*M))


hduX = fits.Header()
hduX['type']='IFunc'
hduX['parity'] = 'Even'
hduX['Subs'] = 'IN IN UP'
hduX['M']= M
hduX['omega']= omega
hduX['l']= l
empty_primary = fits.PrimaryHDU(header=hduX)
hdulX = fits.HDUList([empty_primary])

#ok so this is for a specific omega 
# need to span hs, and need to have different k,k',l

hduXo = fits.Header()
hduXo['type']='IFunc'
hduXo['parity'] = 'Odd'
hduXo['Subs'] = 'IN IN UP'
hduXo['M']= M
hduXo['omega']= omega
hduXo['l']= l
empty_primaryo = fits.PrimaryHDU(header=hduXo)
hdulXo = fits.HDUList([empty_primaryo])

h = np.linspace(.01*T,20*T,2000)

vals_e = np.zeros((len(ks),len(k_primes)),dtype=complex)
vals_o = np.zeros((len(ks),len(k_primes)),dtype=complex)

hduP = fitsio.FITS(direcPhoton+str(l)+'ExtendedOmega.fits')
r_points_gamma = hduP[omega_index]['rpoints_up'][:]
 
psi_gammalomega = hduP[omega_index]['F_points_up'][:]
psi_gammalomega_prime = hduP[omega_index]['z_points_up'][:]
                
        
rs = np.array([m1.r_star_to_r(x,M,tol) for x in r_points_gamma])


def prof_make_tables():
    
        
       
    h_index= round(h[x]*100*(8*np.pi*M))

    if (omega-h[x])>0.0095*T:
        for k in range(len(ks)):
           
            j=(np.abs(ks[k]) - 1/2)

            
            if ks[k]>0:
                hdu = fits.open(direcElectron+str(np.abs(ks[k]))+'ExtendedOmega.fits')
            else:
                hdu = fits.open(direcElectron+'min'+str(np.abs(ks[k]))+'ExtendedOmega.fits')
            
            F_points_xkh = hdu[h_index].data.field('F_points_in')
            G_points_xkh = hdu[h_index].data.field('G_points_in')


            for k_prime in range(len(k_primes)):
                j_prime = (np.abs(k_primes[k_prime]) - 1/2)
                tryA2=Inp.IfunctionsNoM(X,ks[k],X_prime,k_primes[k_prime],X_gamma,l,parity,h[x],omega-h[x],omega,M,n)

                
                vals_e[k][k_prime], vals_o[k][k_prime] = tryA2.IBarplusplusfunc(X,ks[k],X_prime,k_primes[k_prime],psi_gammalomega,psi_gammalomega_prime,l,h[x],omega-h[x],omega,M,rs,r_points_gamma,F_points_xkh,G_points_xkh,nu,mu,lam,GC,c,direcElectron)

    col0 = fits.Column(name='k val',format='D',array=ks)
    col1 = fits.Column(name='k_prime -10',format='M',array=vals_e[:,0])
    col2 = fits.Column(name='k_prime -9',format='M',array=vals_e[:,1])
    col3 = fits.Column(name='k_prime -8',format='M',array=vals_e[:,2])
    col4 = fits.Column(name='k_prime -7',format='M',array=vals_e[:,3])
    col5 = fits.Column(name='k_prime -6',format='M',array=vals_e[:,4])
    col6 = fits.Column(name='k_prime -5',format='M',array=vals_e[:,5])
    col7 = fits.Column(name='k_prime -4',format='M',array=vals_e[:,6])
    col8 = fits.Column(name='k_prime -3',format='M',array=vals_e[:,7])
    col9 = fits.Column(name='k_prime -2',format='M',array=vals_e[:,8])
    col10 = fits.Column(name='k_prime -1',format='M',array=vals_e[:,9])
    col11 = fits.Column(name='k_prime 1',format='M',array=vals_e[:,10])
    col12 = fits.Column(name='k_prime 2',format='M',array=vals_e[:,11])
    col13 = fits.Column(name='k_prime 3',format='M',array=vals_e[:,12])
    col14 = fits.Column(name='k_prime 4',format='M',array=vals_e[:,13])
    col15 = fits.Column(name='k_prime 5',format='M',array=vals_e[:,14])
    col16 = fits.Column(name='k_prime 6',format='M',array=vals_e[:,15])
    col17 = fits.Column(name='k_prime 7',format='M',array=vals_e[:,16])
    col18 = fits.Column(name='k_prime 8',format='M',array=vals_e[:,17])
    col19 = fits.Column(name='k_prime 9',format='M',array=vals_e[:,18])
    col20 = fits.Column(name='k_prime 10',format='M',array=vals_e[:,19])



    cols = fits.ColDefs([col0,col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20])
    table_hdu = fits.BinTableHDU.from_columns(cols)
    table_hdu.header['h']=h[x]

    hdulX.append(table_hdu)


    col0o = fits.Column(name='k val',format='D',array=ks)
    col1o = fits.Column(name='k_prime -10',format='M',array=vals_o[:,0])
    col2o = fits.Column(name='k_prime -9',format='M',array=vals_o[:,1])
    col3o = fits.Column(name='k_prime -8',format='M',array=vals_o[:,2])
    col4o = fits.Column(name='k_prime -7',format='M',array=vals_o[:,3])
    col5o = fits.Column(name='k_prime -6',format='M',array=vals_o[:,4])
    col6o = fits.Column(name='k_prime -5',format='M',array=vals_o[:,5])
    col7o = fits.Column(name='k_prime -4',format='M',array=vals_o[:,6])
    col8o = fits.Column(name='k_prime -3',format='M',array=vals_o[:,7])
    col9o = fits.Column(name='k_prime -2',format='M',array=vals_o[:,8])
    col10o = fits.Column(name='k_prime -1',format='M',array=vals_o[:,9])
    col11o = fits.Column(name='k_prime 1',format='M',array=vals_o[:,10])
    col12o = fits.Column(name='k_prime 2',format='M',array=vals_o[:,11])
    col13o = fits.Column(name='k_prime 3',format='M',array=vals_o[:,12])
    col14o = fits.Column(name='k_prime 4',format='M',array=vals_o[:,13])
    col15o = fits.Column(name='k_prime 5',format='M',array=vals_o[:,14])
    col16o = fits.Column(name='k_prime 6',format='M',array=vals_o[:,15])
    col17o = fits.Column(name='k_prime 7',format='M',array=vals_o[:,16])
    col18o = fits.Column(name='k_prime 8',format='M',array=vals_o[:,17])
    col19o = fits.Column(name='k_prime 9',format='M',array=vals_o[:,18])
    col20o = fits.Column(name='k_prime 10',format='M',array=vals_o[:,19])



    colso = fits.ColDefs([col0o,col1o,col2o,col3o,col4o,col5o,col6o,col7o,col8o,col9o,col10o,col11o,col12o,col13o,col14o,col15o,col16o,col17o,col18o,col19o,col20o])
    table_hduo = fits.BinTableHDU.from_columns(colso)
    table_hduo.header['h']=h[x]

    hdulXo.append(table_hduo)


prof_make_tables()

del vals_e
del vals_o



hdulX.writeto(pathe,overwrite=True)
hdulXo.writeto(patho,overwrite=True)

hdulX.close()
hdulXo.close()
