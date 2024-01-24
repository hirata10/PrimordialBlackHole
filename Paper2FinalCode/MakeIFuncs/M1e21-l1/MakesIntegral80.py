import numpy as np 
from astropy.io import fits
import os
import time
import sys 



output_path = str(sys.argv[1])
l = int(sys.argv[2])
omega_val = np.array(sys.argv[3:])
omega_val = omega_val.astype(float)


hdu_c = fits.open('/users/PCON0003/koivuemily/PrimordialBlackHole/Constants.fits')

M = hdu_c[0].header['M']
Temp = hdu_c[0].header['Temp']
alpha = hdu_c[0].header['alpha']
mu = hdu_c[0].header['mu']


direcElectron = hdu_c[0].header['E_direc']
direcPhoton = hdu_c[0].header['P_direc']



def getRandT(l,omega):
    hduP = fits.open(direcPhoton+str(l)+'ExtendedOmega.fits')
    omega_index = round(omega*100*(8*np.pi*M))
    
    R,T = hduP[omega_index].header['R'],hduP[omega_index].header['T']
    return R,T


print(Temp)
print(Temp/mu)



ks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]
k_primes=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]

def get_energy_index(en):
    return round(en*100*(8*np.pi*M))




#h and omega here dont include the T contribution, so the 8piM cancel in all exponentials 
def term1(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulto'+str(h_index)+'.fits')
        coeff = 2/(np.exp((omega+h))+1)
        print(coeff)
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                #print( I_hdu_e[1].data['k_prime -1'])
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(R)**2


def term2(R,T,h,omega,h_index):
    I= 0.
    
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 1/((np.exp((h))+1)*(np.exp((omega-h))+1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(R)**2

def term3(R,T,h,omega,h_index):
    I= 0.
    coeff =0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 2*np.exp(h)/((np.exp(omega+h)+1)*(np.exp(h)+1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2
            
    return I*coeff * np.abs(R)**2

def term4(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -1/(np.exp(omega)-1)
        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(T)**2

def term5(R,T,h,omega,h_index):
    I= 0.
    coeff=0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 2*np.exp(omega)/((np.exp((omega+h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(T)**2

def term6(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -2*np.exp(h)/((np.exp((h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(T)**2

def term7(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -2/((np.exp((h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                I += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 + np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return I*coeff * np.abs(T)**2

def term8(R,T,h,omega,h_index):
    I= 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
            coeff = 2*(2*np.exp(omega)-1) /((np.exp((h+omega))+1)*(np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulto'+str(h_index)+'.fits')

    
            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)

def term9(R,T,h,omega,h_index):
    I= 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
            coeff = 1/((np.exp((h-omega))+1)*(np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)

def term10(R,T,h,omega,h_index):
    I= 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = 2*np.exp(h)/((np.exp((h+omega))+1)*(np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)

def term11(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = -1/((np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)

def term12(R,T,h,omega,h_index):
    I= 0.
    coeff = 0.
    
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = 2*np.exp(h)/((np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)

def term13(R,T,h,omega,h_index):
    I= 0.
    coeff = 0. 
    
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')==True:
   
            coeff = 2/((np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    I += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    I += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(I*coeff *np.conjugate(T)*R)




def term1part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff = 0.
    
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulto'+str(h_index)+'.fits')
        #print(h)
        coeff = 2/(np.exp((omega+h))+1)
        #print(coeff)
 
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                #print( I_hdu_e[1].data['k_prime -1'])
                #print(str(k_primes[k_prime_ind]))
                Ie += np.real(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]))
                Io += np.real(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]))
    
    #return Ie*coeff *np.abs(R)**2,Io*coeff * np.abs(R)**2
    #return Ie,Io
    return Ie*coeff*np.abs(R)**2,Io*coeff*np.abs(R)**2
    #return Ie*coeff,Io*coeff


def term2part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 1/((np.exp((h))+1)*(np.exp((omega-h))+1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return Ie*coeff * np.abs(R)**2,Io*coeff * np.abs(R)**2

def term3part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff =0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 2*np.exp(h)/((np.exp(omega+h)+1)*(np.exp(h)+1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2
            
    return Ie*coeff * np.abs(R)**2, Io*coeff * np.abs(R)**2

def term4part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff = 0.

    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -1/(np.exp(omega)-1)

        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return Ie*coeff * np.abs(T)**2,Io*coeff * np.abs(T)**2

def term5part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff=0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = 2*np.exp(omega)/((np.exp((omega+h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return Ie*coeff * np.abs(T)**2,Io*coeff * np.abs(T)**2

def term6part(R,T,h,omega,h_index):
    Ie= 0.
    Io=0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -2*np.exp(h)/((np.exp((h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return Ie*coeff * np.abs(T)**2,Io*coeff * np.abs(T)**2

def term7part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        coeff = -2/((np.exp((h))+1)*(np.exp(omega)-1))
        I_hdu_e = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
        I_hdu_o = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
        for k_ind in range(len(ks)):
            for k_prime_ind in range(len(k_primes)):
                Ie += np.abs(I_hdu_e[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2 
                Io += np.abs(I_hdu_o[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])**2

    return Ie*coeff * np.abs(T)**2,Io*coeff * np.abs(T)**2

def term8part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
            coeff = 2*(2*np.exp(omega)-1) /((np.exp((h+omega))+1)*(np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/InUpIn/trialParallelresulto'+str(h_index)+'.fits')

    
            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)

def term9part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
            coeff = 1/((np.exp((omega-h))+1)*(np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)

def term10part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff = 0. 
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = 2*np.exp(h)/((np.exp((h+omega))+1)*(np.exp(h)+1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpUpIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)

def term11part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff = 0.
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = -1/((np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/InInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)

def term12part(R,T,h,omega,h_index):
    Ie= 0.
    Io= 0.
    coeff = 0.
    
    if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
            coeff = -2*np.exp(h)/((np.exp(h)+1)*(np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/PlusPlus/UpInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)

def term13part(R,T,h,omega,h_index):
    Ie= 0.
    Io = 0.
    coeff = 0. 
    
    if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')==True:
        if os.path.exists('Omega'+str(omega)+'T/MinusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')==True:
   
   
            coeff = -2/((np.exp(h)+1)*(np.exp(omega)-1))
            I_hdu_e1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o1 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInUp/trialParallelresulto'+str(h_index)+'.fits')
            I_hdu_e2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInIn/trialParallelresulte'+str(h_index)+'.fits')
            I_hdu_o2 = fits.open('Omega'+str(omega)+'T/MinusPlus/UpInIn/trialParallelresulto'+str(h_index)+'.fits')


            for k_ind in range(len(ks)):
                for k_prime_ind in range(len(k_primes)):
                    Ie += I_hdu_e1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_e2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])
                    Io += I_hdu_o1[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind]*np.conjugate(I_hdu_o2[1].data['k_prime '+str(k_primes[k_prime_ind])][k_ind])

    return np.real(Ie*coeff *np.conjugate(T)*R),np.real(Io*coeff *np.conjugate(T)*R)



def integral_contribution(R,T,h,omega,h_index): 
    I1 = term1(R,T,h,omega,h_index)
    I2 = term2(R,T,h,omega,h_index)
    I3 = term3(R,T,h,omega,h_index)
    I4 = term4(R,T,h,omega,h_index)
    I5 = term5(R,T,h,omega,h_index)
    I6 = term6(R,T,h,omega,h_index)
    I7 = term7(R,T,h,omega,h_index)
    I8 = term8(R,T,h,omega,h_index)
    I9 = term9(R,T,h,omega,h_index)
    I10 = term10(R,T,h,omega,h_index)
    I11 = term11(R,T,h,omega,h_index)
    I12 = term12(R,T,h,omega,h_index)
    I13 = term13(R,T,h,omega,h_index)
    """print(I1)
    print(I2)
    print(I3)
    print(I4)
    print(I5)
    print(I6)
    print(I7)
    print(I8)
    print(I9)
    print(I10)
    print(I11)
    print(I12)
    print(I13)"""
    dh = .01*Temp
    #e = need to figure out what e is in these units
    #e = 0.3028 #in natural units 
    #print(I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+I11+I12+I13)
    
    #print((I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+I11+I12+I13)*dh/(4*np.pi*np.pi))
    I_total = (I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+I11+I12+I13)*dh/(4*np.pi*np.pi)
    return I_total


# In[10]:


def integral_contribution_per_term(R,T,h,omega,h_index): 
    I1e,I1o = term1part(R,T,h,omega,h_index)
    I2e,I2o = term2part(R,T,h,omega,h_index)
    I3e,I3o = term3part(R,T,h,omega,h_index)
    I4e,I4o = term4part(R,T,h,omega,h_index)
    I5e,I5o = term5part(R,T,h,omega,h_index)
    I6e,I6o = term6part(R,T,h,omega,h_index)
    I7e,I7o = term7part(R,T,h,omega,h_index)
    I8e,I8o = term8part(R,T,h,omega,h_index)
    I9e,I9o = term9part(R,T,h,omega,h_index)
    I10e,I10o = term10part(R,T,h,omega,h_index)
    I11e,I11o = term11part(R,T,h,omega,h_index)
    I12e,I12o = term12part(R,T,h,omega,h_index)
    I13e,I13o = term13part(R,T,h,omega,h_index)
    
    """print(I1)
    print(I2)
    print(I3)
    print(I4)
    print(I5)
    print(I6)
    print(I7)
    print(I8)
    print(I9)
    print(I10)
    print(I11)
    print(I12)
    print(I13)"""
    dh = .01*Temp
    
    I_totale = np.array((I1e,I2e,I3e,I4e,I5e,I6e,I7e,I8e,I9e,I10e,I11e,I12e,I13e))*dh/(4*np.pi*np.pi)
    I_totalo = np.array((I1o,I2o,I3o,I4o,I5o,I6o,I7o,I8o,I9o,I10o,I11o,I12o,I13o))*dh/(4*np.pi*np.pi)
    return (I_totale,I_totalo)



dh = .01*Temp

import matplotlib.pyplot as plt 




def partinfoI(omega):
    start = time.time()
    ITotAllHe = np.zeros((2000,13))
    ITotAllHo = np.zeros((2000,13))
    Ie = 0.
    #ITotAllHe = np.zeros((2000))
    #ITotAllHo = np.zeros((2000))
    print(ITotAllHe.shape)
    R,T = getRandT(l,omega*Temp) 
    print(R)
    hs = np.linspace(.01, 20, 2000)
    for h_place in range(2000):
        print(h_place)
        h_index = get_energy_index(hs[h_place]*Temp)
        #print(test)
        #print(hs[h_place])
        results =  integral_contribution_per_term(R,T,hs[h_place],omega,h_place)
        ITotAllHe[h_place][:] = results[0][:]
        print(results[0][1])
        #ITotAllHe[h_place] = results[0]
        #print(ITotAllHe[h_place])
        ITotAllHo[h_place][:] = results[1][:]
        print(results[1][1])
        #ITotAllHo[h_place] = results[1]
        #Ie = results[2]
        #print(time.time()-start) 

    #print(ITotAllHe[:][0])
    I1vale= np.sum(ITotAllHe.T[0])*alpha
    I2vale= np.sum(ITotAllHe.T[1])*alpha
    I3vale= np.sum(ITotAllHe.T[2])*alpha
    I4vale= np.sum(ITotAllHe.T[3])*alpha
    I5vale= np.sum(ITotAllHe.T[4])*alpha
    I6vale= np.sum(ITotAllHe.T[5])*alpha
    I7vale= np.sum(ITotAllHe.T[6])*alpha
    I8vale= np.sum(ITotAllHe.T[7])*alpha
    I9vale= np.sum(ITotAllHe.T[8])*alpha
    I10vale= np.sum(ITotAllHe.T[9])*alpha
    I11vale= np.sum(ITotAllHe.T[10])*alpha
    I12vale= np.sum(ITotAllHe.T[11])*alpha
    I13vale= np.sum(ITotAllHe.T[12])*alpha
    
    
    I1valo= np.sum(ITotAllHo.T[0])*alpha
    I2valo= np.sum(ITotAllHo.T[1])*alpha
    I3valo= np.sum(ITotAllHo.T[2])*alpha
    I4valo= np.sum(ITotAllHo.T[3])*alpha
    I5valo= np.sum(ITotAllHo.T[4])*alpha
    I6valo= np.sum(ITotAllHo.T[5])*alpha
    I7valo= np.sum(ITotAllHo.T[6])*alpha
    I8valo= np.sum(ITotAllHo.T[7])*alpha
    I9valo= np.sum(ITotAllHo.T[8])*alpha
    I10valo= np.sum(ITotAllHo.T[9])*alpha
    I11valo= np.sum(ITotAllHo.T[10])*alpha
    I12valo= np.sum(ITotAllHo.T[11])*alpha
    I13valo= np.sum(ITotAllHo.T[12])*alpha
    
    with open(output_path, "a") as file:
    # write to file
        file.writelines(str(omega)+'\n')
        file.writelines(str(I1vale)+'\t'+ str(I2vale)+'\t'+ str(I3vale)+'\t'+ str(I4vale)+'\t'+ str(I5vale)+'\t'+ str(I6vale)+'\t'+ str(I7vale)+'\t'+ str(I8vale)+'\t'+ str(I9vale)+'\t'+ str(I10vale)+'\t'+ str(I11vale)+'\t'+ str(I12vale)+'\t'+ str(I13vale))
        #file.writelines(str(I1vale))
        file.writelines('\n')
        file.writelines(str(I1valo)+'\t'+ str(I2valo)+'\t'+ str(I3valo)+'\t'+ str(I4valo)+'\t'+ str(I5valo)+'\t'+ str(I6valo)+'\t'+ str(I7valo)+'\t'+ str(I8valo)+'\t'+ str(I9valo)+'\t'+ str(I10valo)+'\t'+ str(I11valo)+'\t'+ str(I12valo)+'\t'+ str(I13valo))
        #file.writelines(str(I1valo))
        file.writelines('\n')
        
    print(I1vale,I1valo)
    return ITotAllHe,ITotAllHo


for i in range(len(omega_val)):
    partinfoI(omega_val[i])


