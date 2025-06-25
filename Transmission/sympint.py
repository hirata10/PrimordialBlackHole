"""Calling format:
python sympint.py <M*mu>

Generates a table of dN/dh dt as a function of Z (column) and h (row)
"""

import sys
import numpy as np

# in these routines, rsub = r-2M (to avoid underflow)

mu = 4.1796514508e-23
alpha = 0.0072973525643

# the matrix to step (F,G) from rsub_prev to rsub_next
def stepmatr(M, h, aZ, k, rsub_prev, rsub_next):
    rsub = (rsub_next+rsub_prev)/2.
    dr = rsub_next-rsub_prev
    r = rsub+2*M
    nu = np.sqrt(rsub/r) # sqrt{1-2M/r}

    dmat = np.array([[-k/r, -(h+aZ/r)/nu+mu],
        [(h+aZ/r)/nu+mu, k/r]])/nu

    return( (np.identity(2)+dr/2.*dmat)@np.linalg.inv(np.identity(2)-dr/2.*dmat) )

def get_T(M, h, aZ, k):

    r_i = 5000*M
    nstep = 1200000
    drr = -18./nstep
    tMat = np.identity(2)
    for i in range(nstep):
        tMat = stepmatr(M,h,aZ,k,r_i,r_i*(1+drr))@tMat
        r_i *= 1+drr

    FGp = np.linalg.solve(tMat, np.array([1j,1]))/2
    a = ((h+mu)/(h-mu))**.25
    Fs = FGp[0]*a
    Gs = FGp[1]/a
    Pout = np.abs(Fs-1j*Gs)**2
    Pin = np.abs(Fs+1j*Gs)**2
    return(1-Pout/Pin)

def get_em_dNdhdt(M, h, Z, kmax=5):
    Rate = 0.
    for k in range(-kmax,kmax+1):
        if k!=0:
            f = 1./(1.+np.exp(8*np.pi*M*h+4*np.pi*alpha*Z))
            Rate += f*get_T(M,h,Z*alpha,k)*2*np.abs(k)
    return(Rate/2./np.pi)

M = 1e21
h = 3.51/(8*np.pi*1e21)

M = float(sys.argv[1])/mu
dh = 0.005/M

Z = [-15,-5,0,5,15]
NZ = 5
Rate = np.zeros(NZ)
for x in range(1,201):
  h = x*dh
  if h>mu:
    dRate = np.zeros(NZ)
    s = '{:7.4f}'.format(M*h)
    for iZ in range(NZ):
        dRate[iZ] = get_em_dNdhdt(M,h,Z[iZ])
        s += ' {:11.5E}'.format(dRate[iZ])
    Rate = Rate + dRate*dh
    print(s)
print('#', Rate*M)
