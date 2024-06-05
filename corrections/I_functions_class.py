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

import cmath
from astropy.io import fits
import fitsio
import time 



# In[4]:


class ElectronWaveFunction:
    
    """
    
    Electron Wave Function class implements the parameters and RK-4 method for the electron's radial part of the 
    wave function.
    
    
    Parameters
    ----------
    
    nu: float
    
    k: float
       quantum number 
    
    h: float
       Energy eigenvalue 
    
    
    mu: float
        mass of the electron
    
    
    lam: #NOTE: Do we even need this here?
    
    
    GC: float
        graviational constant, being set to 1 here
    
    c: float
       speed of light, set to 1 
    
    M: float
       mass of the blackhole
    
    
    tol: float
         tolerance level for various functions implemented 
         
         
         
    Methods
    -------
    
    r_to_r_star(r,M)
        Returns the corresponding value of r_star for a given r value
    
    r_star_to_r_condition(r_star, r, M)
        Returns the condition for the root solver used in the next funciton
        
    r_star_to_r(r_star, M, tol)
        Returns the corresponding r value for a given r_star value
        
    diff_equations(r_star, F, G)
        Returns the right hand side of the differential equations for the F and G radial functions
        
    RK_4(r_initial, r_final, N)
        Returns the set of F and G points after running the RK 4 method to solve the differential equations given 
        in the diff_equations() method. 
        
    """

    def __init__(self, nu, h, k, mu, M, lam, GC, c, tol):
        
        self.nu = nu
        self.k = k
        self.h = h
        self.mu = mu
        self.lam = lam
        self.GC = GC
        self.c = c
        self.M = M 
        self.tol = tol
        self.v = np.sqrt(h**2 -mu**2)/h


    #This function returns the differential equation    
    def diff_eq(self, r_star, F, G):
        
        #this is w.r.t. to r_star coordinate
        
        r = m1.r_star_to_r(r_star, self.M, self.tol) #corresponding r value for the given r_star
        
        dG = ((self.h - self.mu * np.sqrt(1 - (2*self.M / r)))* F - 
            (self.k / r) * np.sqrt(1 - (2*self.M / r)) * G) 
        
        dF = ((-self.mu * np.sqrt(1 - (2*self.M / r)) - self.h) * G + 
            (self.k / r) * np.sqrt(1 - (2*self.M / r)) * F) 


        return np.asarray([dF, dG], dtype = np.complex128)
    
    
    def RK_4(self, r_initial, r_final, N, up):
        
        #this is in terms of r_star
        
        step_size = (r_final - r_initial) / (N-1)
        
        r_points = np.linspace(r_initial, r_final, N)        
        
        #we need to specify the r_points for initial and final so that they match the F and G that come out! 
        
        #Define initial values for F and G
        if self.h>self.mu:
        
            if up == True: #for up scattering states
                F = complex(np.sqrt(self.h + self.mu), 0)
                G = complex(0, -np.sqrt(self.h - self.mu))

            elif up == False: ###NOTE: this would be for the "in" funtions
                F = complex(1,0) ###NOTE: These are the eigenvalues for r_star -> -(infinity) 
                G = complex(0,1)
            
            #r_points = np.arange(r_initial, r_final, N)
        if self.h<self.mu:
            
            if up == True:
                F = complex((self.h+self.mu)/np.sqrt(self.mu**2 -self.h**2),0)
                G = complex(1.,0)
                r_points = np.linspace(np.maximum(r_initial,r_final),np.minimum(r_initial, r_final), N) 
            elif up == False:
                F = complex(0,0)
                G = complex(0,0)


        #Arrays to keep track of r, F and G points:

        F_points = []
        G_points = []
        
        
        for r in r_points:  
            #Update F_points and G_points arrays:
            F_points.append(F)
            G_points.append(G)


            #Calculate all the slopes (each variable has a separate k-slope):
            (k1F, k1G) = tuple(step_size * self.diff_eq(r, F, G))

            (k2F, k2G) = tuple(step_size * self.diff_eq(r + 0.5*step_size, F + 0.5*k1F, G + 0.5*k1G))

            (k3F, k3G) = tuple(step_size * self.diff_eq(r + 0.5*step_size, F + 0.5*k2F, G + 0.5*k2G))

            (k4F, k4G) = tuple(step_size * self.diff_eq(r + step_size, F + k3F, G + k3G))


            #Calculate next F_point and G_point
            F = F + (k1F + 2*k2F + 2*k3F + k4F) / 6
            G = G + (k1G + 2*k2G + 2*k3G + k4G) / 6
            
        #want to normalize the last points of things 
        
        if self.h>self.mu: 
            if up == True:
                # want F and g to match 41 in the paper 
                #A is the normalization here 

                # want this for rstar at negative infinity to get this
                # rstar for up is from positive infinity to negative infinity so rstar negative infinity should be the last point
                A = 2 * np.sqrt(self.h) * np.exp(complex(0,1)*self.h*r_points[-1])/ (F_points[-1]+ complex(0,1)*G_points[-1])

                #A* arijit's F and G = F_up and G_up from the paper
                #use this to isolate T! then do a similar thing for up==False 
                F_points = A*np.array(F_points)
                G_points = A*np.array(G_points)

            if up ==False:

                zeta = self.mu**2 * self.M/np.sqrt(self.h**2 -self.mu**2)

                #need a normalization here 

                #this assumes that r_initial and r_final are chosen so the integrater starts at r_points[0] and ends at r_points[N]
                #therefore, you have to pass the initial conditions appropriately for this 
                #ie rstar = infinity is the last point 

                B = 2*self.h*np.sqrt(self.v)*np.exp(complex(0,-1)*zeta*np.log(r_points[-1]/2/self.M))*np.exp(complex(0,-1)*np.sqrt(self.h**2-self.mu**2)*r_points[-1])/(np.sqrt(self.h-self.mu)*F_points[-1] -complex(0,1)*np.sqrt(self.h+self.mu)*G_points[-1])

                F_points = B*np.array(F_points)
                G_points = B*np.array(G_points)
                
        if self.h<self.mu:
            if up ==True:
                
                A = 2*complex(0,-1)*np.sqrt(self.h)*np.exp(complex(0,1)*self.h*r_points[-1])/(complex(0,-1)*F_points[-1]+G_points[-1])
                F_points = A*np.array(F_points)
                G_points = A*np.array(G_points)
                
        return (r_points, F_points, G_points)
    
    
    
    
    
    def get_R_and_T_coeff(self,r_up,F_points_up,G_points_up,r_in,F_points_in,G_points_in):
        
        #assumes the inputs to this are already normalized by A and B above 
        if self.h<self.mu:
            R_half_k_h=1.
            T_half_k_h =0. 
            
            delta = complex(0,-1)*.5*np.log((complex(0,1)*F_points_up[-1]+G_points_up[-1])/(2*complex(0,1)*np.sqrt(self.h)*np.exp(complex(0,-1)*self.h*r_up[-1])))
            #match delta to small solution! 
        
        if self.h>self.mu:
            delta = 0 
            zeta = self.mu**2 * self.M/np.sqrt(self.h**2 -self.mu**2)

            #want F_points_up as rstar goes to infinity and up has rstar starting at positive infinity and going to negative 
            a1=G_points_up[0]*np.sqrt(self.v)
            b1=np.exp(complex(0,-1)*zeta*np.log(r_up[0]/2/self.M))
            c1=np.exp(complex(0,-1)*np.sqrt(self.h**2 -self.mu**2)*r_up[0])
            d1=1/(complex(0,-1)*np.sqrt(self.h-self.mu))



            T_half_k_h = a1*b1*c1*d1




            #r_in integrates from negative infinity to positive infinity so that if we want rstar = infinity, that is the last argument 

            a1 = (np.sqrt(self.h-self.mu)*F_points_in[-1] +complex(0,1)*np.sqrt(self.h+self.mu)*G_points_in[-1])
            b1=np.exp(complex(0,-1)*zeta*np.log(r_in[-1]/2/self.M))
            c1=np.exp(complex(0,-1)*np.sqrt(self.h**2 - self.mu**2)*r_in[-1])
            d1=1/2/self.h/np.sqrt(self.v)

            R_half_k_h = a1*b1*c1*d1

        return R_half_k_h,T_half_k_h,delta


class PhotonWaveFunction:
    
    """
    Methods and their descriptions
    
    
    """
    
    def __init__(self, M, omega, l, tol):
        
        self.M = M
        self.omega = omega
        self.l = l
        self.tol = tol
    


    def diff_eq(self, r_star, z_vector):
    
        f = z_vector[0]
        z = z_vector[1]

        f_prime = z

        r = m1.r_star_to_r(r_star, self.M, self.tol)

        z_prime = (-self.omega**2 + ((self.l * (self.l+1) *(1 - (2*self.M)/r)) / (r*r))) * f

        return np.asarray([f_prime, z_prime], dtype = np.complex128)
    
    
    
    
    def RK_4(self, r_initial, r_final, N, up):
        
        step_size = (r_final - r_initial) / (N-1)
        
        r_points = np.linspace(r_initial, r_final, N)
        
        #Initial Conditions
        if up == True: #for up scattering states 
            f_i = complex(1.0, 0)
            z_i = complex(0, self.omega)    
            
        elif up == False: #for in scattering states
            f_i = complex(1.0, 0) 
            z_i = complex(0, -self.omega) #what arijit originally had but things work better if it is the other way 
            #z_i = complex(0, self.omega)
        
        
        f_points = [] 
        z_points = []
        
        z_vector = np.asarray([f_i, z_i])
        
        f_points_prime = []
        
        
        for r in r_points:
    
            f_points.append(z_vector[0])
            z_points.append(z_vector[1])
            
            prime = self.diff_eq(r, z_vector)
            #f_points_prime.append(prime[0])
            
            k1 = step_size * self.diff_eq(r, z_vector)

            k2 = step_size * self.diff_eq(r + 0.5*step_size, z_vector + 0.5*k1)

            k3 = step_size * self.diff_eq(r + 0.5*step_size, z_vector + 0.5*k2)

            k4 = step_size * self.diff_eq(r + step_size, z_vector + k3)

            z_vector += (k1 + 2*k2 + 2*k3 + k4) / 6
            
        if up ==True:
            
            #evaluated at rstar -> neg infinity
            b = (2*np.complex(0,1)*self.omega*np.exp(np.complex(0,1)*self.omega*r_points[-1]) )/(f_points[-1]*np.complex(0,1)*self.omega + z_points[-1])
            
            f_points = b*np.array(f_points)
            z_points = b*np.array(z_points)
            #f_points_prime = b*np.array(f_points_prime)
            
        
        if up == False:
            
            #evaluated at rstar -> infinity
            a = (2*np.complex(0,-1)*self.omega*np.exp(np.complex(0,-1)*self.omega*r_points[-1]))/(f_points[-1]*np.complex(0,-1)*self.omega + z_points[-1])
            
            f_points = a*np.array(f_points)
            z_points = a*np.array(z_points)
            #f_points_prime = a*np.array(f_points_prime)
            
        return(r_points, f_points, z_points, f_points_prime)
    
    def get_R_and_T_coeff(self,r_up,F_points_up,F_points_up_prime,r_in,F_points_in,F_points_in_prime):
        
        T = F_points_up[0]*np.exp(np.complex(0,-1)*self.omega*r_up[0])
        
        T2 = F_points_in[0]*np.exp(np.complex(0,1)*self.omega*r_in[0])
        
        R = (np.complex(0,1)*self.omega*F_points_in[-1] + F_points_in_prime[-1])/(2*np.complex(0,1)*self.omega*np.exp(np.complex(0,1)*self.omega*r_in[-1]))
        
        #R = (np.complex(0,-1)*self.omega*F_points_up[-1] + F_points_up_prime[-1]) / (2*np.complex(0,1)*self.omega* np.exp(np.complex(0,-1)*self.omega*r_up[-1])*(1+2*np.complex(0,1)*cmath.phase(T))) #bad!!!
        
        return R,T


class IfunctionsNoM:
    
    def __init__(self,X,k,X_prime,k_prime,X_gamma,l,parity,h,h_prime,omega,M,n):
        self.X=X
        self.k=k
        self.X_prime=X_prime
        self.k_prime=k_prime
        self.X_gamma=X_gamma
        self.l=l
        self.h=h
        self.h_prime=h_prime
        self.omega=omega
        self.M=M
        
        self.j=(np.abs(k) - 1/2)
        
        self.j_prime = (np.abs(k_prime) - 1/2)
        
        self.coeff_no_m_even = (complex(0,-1)/(np.sqrt(4*h*h_prime))) * self.Delta_no_m(k,k_prime,l)
        #self.coeff_m_even = ((-1)**(m+.5))*float(wigner_3j(self.j,self.j_prime,l,-m,m_prime,m_gamma))
        
        self.coeff_no_m_odd = -(1./np.sqrt(4*h*h_prime)) * self.Pi_no_m(k,k_prime,l)
        #self.step_size_for_inte = int(n*240000)
        
        
        self.omega_index = round(self.omega*100*(8*np.pi*self.M) )
        self.h_index = round(self.h*100*(8*np.pi*self.M) )
        self.h_prime_index = round(self.h_prime*100*(8*np.pi*self.M) )
    
    
    def Delta_no_m(self,k,k_prime,l):
        s=0
        s_prime=0
        
        if k !=0:
            s= (k/np.abs(k))
        if k_prime!=0:
            s_prime=k_prime/np.abs(k_prime)
        
        j=(np.abs(k) - 1/2)
        
        j_prime = (np.abs(k_prime) - 1/2)
        
        threej=float(wigner_3j(j,j_prime,l,1/2,-1/2,0))
        
        x=(1/2) * np.sqrt((2*j+1)*(2*j_prime +1)*(2*l+1)/(4*np.pi))*threej* (1 + s*s_prime*(-1)**(j-j_prime+l))

        return x
    
    def Pi_no_m(self,k,k_prime,l):
        
        s=0
        s_prime=0
        
        if k !=0:
            s= (k/np.abs(k))
        if k_prime!=0:
            s_prime=k_prime/np.abs(k_prime)
        
        j=(np.abs(k) - 1/2)
        
        j_prime = (np.abs(k_prime) - 1/2)
        
        
        w= float(wigner_3j(j,j_prime,l,-1/2,-1/2,1))
        x = s*np.sqrt((2*j+1)*(2*j_prime+1)*(2*l+1)/(4*np.pi))*w*(1-s*s_prime*(-1)**(j-j_prime+l))
        
        return x
    
    
    
    def IBarminusplusfunc(self,X,k,X_prime,k_prime,psi_gammalomega,psi_gammalomega_prime,l,h,h_prime,omega,M,r,r_points, F_points_xminkh, G_points_xminkh,nu,mu,lam,GC,c,direcElectron):
        
        #x=0 in, x=1 up
        
        
        #hdu_c = fits.open('/users/PCON0003/koivuemily/PrimordialBlackHole/Constants.fits')


        #nu =  hdu_c[0].header['nu']
        #mu = hdu_c[0].header['mu']
        #lam = hdu_c[0].header['lam']
        #GC = hdu_c[0].header['GC']
        #c = hdu_c[0].header['c']
        
        #direcElectron = hdu_c[0].header['E_direc']
        #direcPhoton = hdu_c[0].header['P_direc']


        
        if X==0:
            uptrue = False
        else:
            uptrue = True
        if X_prime==0:
            uptrue_Xprime = False
        else:
            uptrue_Xprime = True
                    
        
       # direcElectron = '/users/PCON0003/koivuemily/PrimordialBlackHole/ElectronWaveFunctionFits/'

        if k_prime<0:
            hdu = fitsio.FITS(direcElectron+str(np.abs(k_prime))+'ExtendedOmega.fits')
        else:
            hdu = fitsio.FITS(direcElectron+'min'+str(np.abs(k_prime))+'ExtendedOmega.fits')
        
        #time1= time.time()
        if uptrue_Xprime == True: 
            F_points_xminkh_prime = hdu[self.h_prime_index]['F_points_up'][:]
            G_points_xminkh_prime = hdu[self.h_prime_index]['G_points_up'][:]
        else: 
            F_points_xminkh_prime = hdu[self.h_prime_index]['F_points_in'][:]
            G_points_xminkh_prime = hdu[self.h_prime_index]['G_points_in'][:]
        #time2 = time.time()
        hdu.close()
        
        g_xminkh = np.array(G_points_xminkh)
        gstar_xprime_minkprime_hprime = np.conjugate(np.array(G_points_xminkh_prime))
        
        f_xminkh = np.array(F_points_xminkh)
        fstar_xprime_minkprime_hprime = np.conjugate(np.array(F_points_xminkh_prime))
            
        rdependentpart_top = (1-2*M/r)/(r**2*np.sqrt(2*omega**3))
        rdependentpart_bottom = np.sqrt(1-2*M/r)/(r*np.sqrt(2*omega))
        
        topline = (g_xminkh*fstar_xprime_minkprime_hprime - f_xminkh*gstar_xprime_minkprime_hprime)*np.sqrt(l*(l+1))*psi_gammalomega*rdependentpart_top
           
        botline = (g_xminkh*fstar_xprime_minkprime_hprime + f_xminkh*gstar_xprime_minkprime_hprime)*((k-k_prime)/(omega*np.sqrt(l*(l+1))))*psi_gammalomega_prime*rdependentpart_bottom
        
       
        inte_e = 0.
        rstar_step = r_points[0]-r_points[1]   #switched order to get positive steps in
        inte_e = (np.sum(topline[1:-1]+botline[1:-1])+.5*(topline[0]+topline[-1]+botline[0]+botline[-1]))*rstar_step
        
        
        rdependentpart = np.sqrt((1-2*M/r))/(2*r*np.sqrt(2*omega))
        
        line = (g_xminkh * fstar_xprime_minkprime_hprime + f_xminkh*gstar_xprime_minkprime_hprime)*psi_gammalomega*rdependentpart
        
        inte_o = 0.
        inte_o = (np.sum(line[1:-1])+.5*(line[0]+line[-1]))*rstar_step
        
        return (inte_e*self.coeff_no_m_even,inte_o*self.coeff_no_m_odd)#
        
    def IBarplusplusfunc(self,X,k,X_prime,k_prime,psi_gammalomega,psi_gammalomega_prime,l,h,h_prime,omega,M,r,r_points, F_points_xkh, G_points_xkh,nu,mu,lam,GC,c,direcElectron):
        
        #x=0 in, x=1 up
        
        #hdu_c = fits.open('/users/PCON0003/koivuemily/PrimordialBlackHole/Constants.fits')


        #nu =  hdu_c[0].header['nu']
        #mu = hdu_c[0].header['mu']
        #lam = hdu_c[0].header['lam']
        #GC = hdu_c[0].header['GC']
        #c = hdu_c[0].header['c']
        
        #direcElectron = hdu_c[0].header['E_direc']

        """
        #far from BH and then work towards that (up)
        #r_initial = 5000.*M #in rstar
        #r_final = -70.*M  #in rstar 
        
        #check
        r_initial = lim1*M #in rstar
        r_final = lim2*M  #in rstar 
        #
        #r_final_gamma = 5000.*M
        #r_initial_gamma = -70.*M
        
        
         #up for photon 
        r_initial_gamma = lim1*M
        #r_final_gamma = -70.*M
        
        
        #check 
        r_final_gamma=lim2*M
        
        #r_initial_gamma = 5000.*M
        #r_final_gamma = 1000.*M
     """   
        
        if X==0:
            uptrue = False
        else:
            uptrue = True
        if X_prime==0:
            uptrue_Xprime = False
        else:
            uptrue_Xprime = True
            
        
        #might pass this in instead. will ask 
        #try0 = ElectronWaveFunction(nu, h, k, mu, M, lam, GC, c, tol)
        #r_points,F_points_xkh, G_points_xkh = try0.RK_4(r_initial, r_final, self.step_size_for_inte, up = uptrue)
        
        #try1 = ElectronWaveFunction(nu, h_prime, -k_prime, mu, M, lam, GC, c, tol)
       # r_points_prime,F_points_xminkh_prime, G_points_xminkh_prime = try1.RK_4(r_initial, r_final, self.step_size_for_inte, up = uptrue_Xprime)
        
        
        #phot_try0= PhotonWaveFunction(M,omega,l,tol)
        #r_points_gamma, f_points_gammalomega, z_points_gammalomega, f_points_gammalomega_prime = phot_try0.RK_4(r_initial_gamma, r_final_gamma, self.step_size_for_inte, up = uptrue_Xgamma)
        
        #phot_try1= PhotonWaveFunction(M,omega,l,tol)
        #r_points_gamma_prime, f_points_gammalomega_prime, z_points_gammalomega_prime = phot_try0.RK_4(r_initial_gamma, r_final_gamma, N, up = uptrue_Xgamma_prime)
        
        #want to write a script that shows intermediate results, can put in plots about transmission coefficients 
        
        
        """direcPhoton = '/fs/scratch/PCON0003/emily/PhotonWaveFunctionFits/'
        
        #hdu = fits.open(direcPhoton+str(l)+'.fits')
        
        #r_points_gamma = hdu[self.omega_index].data.field('rpoints_up')
        #if uptrue_Xgamma == True: 
        #    psi_gammalomega = hdu[self.omega_index].data.field('F_points_up')
        #    psi_gammalomega_prime = hdu[self.omega_index].data.field('z_points_up')
        #else: 
        #    psi_gammalomega = hdu[self.omega_index].data.field('F_points_in')
        #    psi_gammalomega_prime = hdu[self.omega_index].data.field('z_points_up')"""
         
        #direcElectron = '/users/PCON0003/koivuemily/PrimordialBlackHole/ElectronWaveFunctionFits/'
        #if k>0:
        #    hdu = fits.open(direcElectron+str(np.abs(k))+'ExtendedOmega.fits')
        #else:
        #    hdu = fits.open(direcElectron+'min'+str(np.abs(k))+'ExtendedOmega.fits')
        
        #r_points = hdu[self.h_index].data.field('rpoints_up')
        #if uptrue == True: 
        #    F_points_xkh = hdu[self.h_index].data.field('F_points_up')
        #    G_points_xkh = hdu[self.h_index].data.field('G_points_up')
        #else: 
        #    F_points_xkh = hdu[self.h_index].data.field('F_points_in')
        #    G_points_xkh = hdu[self.h_index].data.field('G_points_in')
        if k_prime<0:
            hdu = fitsio.FITS(direcElectron+str(np.abs(k_prime))+'ExtendedOmega.fits')
        else:
            hdu = fitsio.FITS(direcElectron+'min'+str(np.abs(k_prime))+'ExtendedOmega.fits')
        
        #r_points = hdu[self.h_index+1]['rpoints_up']
        if uptrue_Xprime == True: 
            F_points_xminkh_prime = hdu[self.h_prime_index]['F_points_up'][:]
            G_points_xminkh_prime = hdu[self.h_prime_index]['G_points_up'][:]
        else: 
            F_points_xminkh_prime = hdu[self.h_prime_index]['F_points_in'][:]
            G_points_xminkh_prime = hdu[self.h_prime_index]['G_points_in'][:]
        
        #need analytic continuation of wavefunctions! 
        #can just add the analytic piece at the end! so need to do that also
            
        gstar_xkh = np.conjugate(np.array(G_points_xkh))
        gstar_xprime_minkprime_hprime = np.conjugate(np.array(G_points_xminkh_prime))
        
        fstar_xkh = np.conjugate(np.array(F_points_xkh))
        fstar_xprime_minkprime_hprime = np.conjugate(np.array(F_points_xminkh_prime))
            
        #psi_gammalomega = np.array(f_points_gammalomega)
        #psi_gammalomega_prime = np.array(z_points_gammalomega)
        
        #r = np.array([m1.r_star_to_r(x,M,tol) for x in r_points])
        rdependentpart_top = (1-2*M/r)/(r**2*np.sqrt(2*omega**3))
        rdependentpart_bottom = np.sqrt(1-2*M/r)/(r*np.sqrt(2*omega))
        
        
        topline = (-1*gstar_xkh*gstar_xprime_minkprime_hprime + fstar_xkh*fstar_xprime_minkprime_hprime)*np.sqrt(l*(l+1))*psi_gammalomega*rdependentpart_top
           
        botline = (fstar_xkh*fstar_xprime_minkprime_hprime + gstar_xkh*gstar_xprime_minkprime_hprime)*((k-k_prime)/(omega*np.sqrt(l*(l+1))))*psi_gammalomega_prime*rdependentpart_bottom
        
        
        #print(topline, botline)
        #return topline+botline
        
        #midpoint method
        intee = 0.
        rstar_step = r_points[0]-r_points[1] #switched order to get positive steps in
        
        #should we look at F and G for mixed up and in cases 
        
        #for index in range(len(topline)-1):
        #    intee += .5*(topline[index]+topline[index+1] + botline[index]+botline[index+1])*rstar_step
            
        intee = (np.sum(topline[1:-1]+botline[1:-1])+.5*(topline[0]+topline[-1]+botline[0]+botline[-1]))*rstar_step

        rdependentpart = np.sqrt((1-2*M/r))/(2*r*np.sqrt(2*omega))
        
        #print(r_points[0],r_points_gamma[0])

        
        line = (fstar_xkh * fstar_xprime_minkprime_hprime + gstar_xkh*gstar_xprime_minkprime_hprime)*psi_gammalomega*rdependentpart
        
        #midpoint method
        inteo = 0.
        
        inteo = (np.sum(line[1:-1])+.5*(line[0]+line[-1]))*rstar_step
 
    
        
        return (intee*self.coeff_no_m_even,inteo*self.coeff_no_m_odd)  

    

