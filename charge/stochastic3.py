# Calling format:
# python correction_factor.py [PBH mass] [step size]
# nominally step size = 1 but larger steps integrate faster

import sys
import numpy as np
from scipy.special import spherical_jn
import I_functions_class_charge as Inp

nu = 1.0
mu = 4.1796514508e-23
lam = 1.0
GC = 1.0
c = 1.0
tol = 1.0e-10
M = float(sys.argv[1])
r_initial = 4000*M
r_final = -2000*M
Temp = 1./(8*np.pi*M)

h = np.linspace(0.01*Temp, 30*Temp, 3000)
dh = 0.01*Temp
kmax=5

trunc = {'r': 3*np.array([r_final, r_final+1000*M, r_initial-1000*M, r_initial])} # avoid truncation
alpha = 0.0072973525643

NPT = 160000

NumTot = 0.
DenTot1 = 0.
DenTot2 = 0.
DenTot = 0.
TT = np.zeros(3)
firststep = True
for x in range(int(sys.argv[2])-1,3000,int(sys.argv[2])):
    if h[x]>mu + 9*mu*M/r_initial:
        dNum = dDen1 = dDen2 = 0.
        f_em = 1./(1.+np.exp(h[x]/Temp))
        dfdh = -1./(2*np.cosh(h[x]/Temp/2))**2/Temp
        v = np.sqrt(h[x]**2-mu**2)/h[x]
        #print('^', x, h[x], mu, v)
        for k in list(range(-kmax,0))+list(range(1,kmax+1)):
            for Z in range(-1,2):
                EWF = Inp.ElectronWaveFunction(nu, h[x], k, mu, M, lam, GC, c, tol, aZ=alpha*Z, trunc=trunc)
                r_points, F_points, G_points = EWF.RK_4(r_initial, r_final, NPT, up = True)
                r_points_in, F_points_in, G_points_in = EWF.RK_4(r_final, r_initial, NPT, up = False)
                R,T,delta = EWF.get_R_and_T_coeff(r_points,F_points,G_points,r_points_in,F_points_in,G_points_in)
                TT[Z+1] = np.abs(T)**2
            dTTdZ = (TT[2]-TT[0])/2.
            #for p in range(-1,2,2):
            #    EWF = Inp.ElectronWaveFunction(nu, h[x]+dh*p, k, mu, M, lam, GC, c, tol, aZ=0., trunc=trunc)
            #    r_points, F_points, G_points = EWF.RK_4(r_initial, r_final, NPT, up = True)
            #    r_points_in, F_points_in, G_points_in = EWF.RK_4(r_final, r_initial, NPT, up = False)
            #    R,T,delta = EWF.get_R_and_T_coeff(r_points,F_points,G_points,r_points_in,F_points_in,G_points_in)
            #    TT[p+1] = np.abs(T)**2
            #dTTdh = (TT[2]-TT[0])/2./dh
            print('# {:9.6f} {:3d} {:11.5E} {:12.5E}'.format(h[x]/Temp, k, TT[1], dTTdZ))
            sys.stdout.flush()
            dNum = dNum + 2./np.pi*np.abs(k)*TT[1]*f_em*(1.-TT[1]*f_em)
            dDen1 = dDen1 + 4./np.pi*np.abs(k)*(-alpha/(2*M)*TT[1]*dfdh)
            dDen2 = dDen2 + 4./np.pi*np.abs(k)*( - dTTdZ*f_em)
        print('{:9.6f} {:12.5e} {:12.5e} {:12.5e} {:12.5e}'.format(h[x]/Temp, dNum, dDen1, dDen2, dDen1+dDen2))
        hstep = hstep_ = int(sys.argv[2])*(h[1]-h[0])
        if firststep:
            hstep = hstep_/2 + h[x]-mu
        NumTot += dNum*hstep
        DenTot1 += dDen1*hstep
        DenTot2 += dDen2*hstep
        DenTot += (dDen1+dDen2)*hstep
        firststep = False

print('# M =', M)
print('# step size =', int(sys.argv[2]))
print('# total rates = ', NumTot, DenTot1, DenTot2, DenTot)
print('# <Z^2> =', NumTot/DenTot)
