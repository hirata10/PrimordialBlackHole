# Calling format:
# python correction_factor.py [PBH mass] [omega/Temp] [step size]
# nominally step size = 1 but larger steps integrate faster

import sys
import numpy as np
from scipy.special import spherical_jn
import I_functions_class as Inp

def ampcorr(xi, l):
   j = spherical_jn(l, xi, derivative=False)
   jp = spherical_jn(l, xi, derivative=True)
   xj = xi*j
   xjp = xi*jp+j
   return(xj**2+xjp**2)

# gets fraction of integral up to xi, and then the total
def frac_I(xi, l, v):
   K = 400; N = 10000
   z = np.linspace(1/N,K-1/N,N*K)*xi
   integrand = spherical_jn(l, z, derivative=False)/z*np.exp(-1j/v*z)
   return np.sum(integrand[:N])/np.sum(integrand), np.sum(integrand)

nu = 1.0
mu = 4.1796514508e-23
lam = 1.0
GC = 1.0
c = 1.0
tol = 1.0e-10
M = float(sys.argv[1])
r_initial = 2000*M
r_final = -70*M
Temp = 1./(8*np.pi*M)

h = np.linspace(0.01*Temp, 20*Temp, 2000)
kmax=2
lmax=5
totEm = np.zeros((lmax+1,))
cutEm = np.zeros((lmax+1,))
omega = float(sys.argv[2])*Temp

xi_ = r_initial*omega

for x in range(0,2000,int(sys.argv[3])):
    if h[x]>mu:
        v = np.sqrt(h[x]**2-mu**2)/h[x]
        print('^', x, h[x], mu, v)
        Gamma = 0.
        for k in list(range(-kmax,0))+list(range(1,kmax+1)):
            EWF = Inp.ElectronWaveFunction(nu, h[x], k, mu, M, lam, GC, c, tol)
            r_points, F_points, G_points = EWF.RK_4(r_initial, r_final, 240000, up = True)
            r_points_in, F_points_in, G_points_in = EWF.RK_4(r_final, r_initial, 240000, up = False)
            R,T,delta = EWF.get_R_and_T_coeff(r_points,F_points,G_points,r_points_in,F_points_in,G_points_in)
            Gamma += 2*np.abs(k)*np.abs(T)**2/(1+np.exp(h[x]/Temp))
            sys.stdout.flush()
        for L in range(1,lmax+1):
            a = ampcorr(xi_,L)
            thisfrac, integrand = frac_I(xi_,L,v)
            totEm[L] += np.abs(integrand)**2*Gamma*v**2
            cutEm[L] += np.abs(integrand)**2*Gamma*v**2*np.abs(thisfrac)**2/a
            print(':', L, np.abs(thisfrac)**2/a)
            sys.stdout.flush()

        print('>>> {:4d} {:6.4f} {:6.4f} {:12.5E}'.format(x, h[x]/Temp, v, Gamma))
        print(totEm)
        print(cutEm)
        print('')
        sys.stdout.flush()

print('# M =', M)
print('# omega/Temp =', omega/Temp)
print('# step size =', int(sys.argv[3]))
for L in range(1,lmax+1): print('{:2d} {:10.6f}'.format(L, cutEm[L]/totEm[L]))
