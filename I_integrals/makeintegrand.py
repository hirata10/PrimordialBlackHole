import os
import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve
import module1 as m1 # module1 contains the function to convert r_star to r\n",
import cmath
from astropy.io import fits
import I_functions_class as Inp
from importlib import reload
import sys
import gc 
reload(Inp)

x = int(sys.argv[1])
y = int(sys.argv[2])


hdu_c = fits.open('/users/PCON0003/bowenchen12686/ondemand/data/sys/myjobs/constants/Constants4.fits')
nu =  hdu_c[0].header['nu']
mu = hdu_c[0].header['mu']
lam = hdu_c[0].header['lam']
GC = hdu_c[0].header['GC']
c = hdu_c[0].header['c']
tol = hdu_c[0].header['tol']
direcPhoton = hdu_c[0].header['P_direc']
direcElectron = hdu_c[0].header['E_direc']
alpha = hdu_c[0].header['alpha']



M = hdu_c[0].header['M']
T = hdu_c[0].header['Temp']
hdu_c.close()


ks=[-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10]

h = np.linspace(.01*T,20*T,2000)
h_index= round(h[x]*100*(8*np.pi*M))


h_prime = np.linspace(.01*T,20*T,2000)
h_prime_index= round(h_prime[y]*100*(8*np.pi*M))





temp_data21 = [[] for _ in range(11)] #F up+kh
temp_data22 = [[] for _ in range(11)] #F in+kh
temp_data23 = [[] for _ in range(11)] #G up+kh
temp_data24 = [[] for _ in range(11)] #G in+kh
temp_data25 = [[] for _ in range(11)] #R+k
temp_data26 = [[] for _ in range(11)] #T+k
temp_data27 = [[] for _ in range(11)] #F up+kh1
temp_data28 = [[] for _ in range(11)] #F in+kh1
temp_data29 = [[] for _ in range(11)] #G up+kh1
temp_data210 = [[] for _ in range(11)] 
temp_data211 = [[] for _ in range(11)] 
temp_data212 = [[] for _ in range(11)] 
temp_data213 = [[] for _ in range(11)] 
temp_data214 = [[] for _ in range(11)] 
temp_data215 = [[] for _ in range(11)] 
for k in range(1, 11):
    file_path = f"{direcElectron}{k}ExtendedOmega.fits"
    with fits.open(file_path) as hdu:
        data21 = hdu[h_index].data.field('F_points_in')
        data22 = hdu[h_index].data.field('G_points_in')
        data23 = hdu[h_index].data.field('F_points_up')
        data24 = hdu[h_index].data.field('G_points_up')
        data25 = hdu[h_index].header['R']
        data26 = hdu[h_index].header['T']
        data27 = hdu[h_prime_index].data.field('F_points_in')
        data28 = hdu[h_prime_index].data.field('G_points_in')
        data29 = hdu[h_prime_index].data.field('F_points_up')
        data210 = hdu[h_prime_index].data.field('G_points_up')
        data211 = hdu[h_prime_index].header['R']
        data212 = hdu[h_prime_index].header['T']
        data213 = hdu[h_index].header['delta']
        data214 = hdu[h_prime_index].header['delta']
        data215 = hdu[h_prime_index].data.field('rpoints_up')

        temp_data21[k].append(data21)
        temp_data22[k].append(data22)
        temp_data23[k].append(data23)
        temp_data24[k].append(data24)
        temp_data25[k].append(data25)
        temp_data26[k].append(data26)
        temp_data27[k].append(data27)
        temp_data28[k].append(data28)
        temp_data29[k].append(data29)
        temp_data210[k].append(data210)
        temp_data211[k].append(data211)
        temp_data212[k].append(data212)
        temp_data213[k].append(data213)
        temp_data214[k].append(data214)
        temp_data215[k].append(data215)

F_points_xkh_in = np.array([np.array(data, dtype=complex) for data in temp_data21 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_in = np.array([np.array(data, dtype=complex) for data in temp_data22 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_up = np.array([np.array(data, dtype=complex) for data in temp_data23 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_up = np.array([np.array(data, dtype=complex) for data in temp_data24 if len(data) > 0], dtype=object).squeeze(axis=1)
Rk = np.array([np.array(data, dtype=complex) for data in temp_data25 if len(data) > 0], dtype=object).squeeze(axis=1)
Tk = np.array([np.array(data, dtype=complex) for data in temp_data26 if len(data) > 0], dtype=object).squeeze(axis=1)
Dk = np.array([np.array(data, dtype=complex) for data in temp_data213 if len(data) > 0], dtype=object).squeeze(axis=1)

F_points_xkh_prime_in = np.array([np.array(data, dtype=complex) for data in temp_data27 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime_in = np.array([np.array(data, dtype=complex) for data in temp_data28 if len(data) > 0], dtype=object).squeeze(axis=1)
F_points_xkh_prime_up = np.array([np.array(data, dtype=complex) for data in temp_data29 if len(data) > 0], dtype=object).squeeze(axis=1)
G_points_xkh_prime_up = np.array([np.array(data, dtype=complex) for data in temp_data210 if len(data) > 0], dtype=object).squeeze(axis=1)
Rk_prime = np.array([np.array(data, dtype=complex) for data in temp_data211 if len(data) > 0], dtype=object).squeeze(axis=1)
Tk_prime = np.array([np.array(data, dtype=complex) for data in temp_data212 if len(data) > 0], dtype=object).squeeze(axis=1)
Dk_prime = np.array([np.array(data, dtype=complex) for data in temp_data214 if len(data) > 0], dtype=object).squeeze(axis=1)
r_points_gamma = np.array(temp_data215[1], dtype=float).squeeze()
rs = np.array([m1.r_star_to_r(x,M,tol) for x in r_points_gamma])

del temp_data21, temp_data22, temp_data23, temp_data24,temp_data25,temp_data26,temp_data27,temp_data28,temp_data29,temp_data210,temp_data211,temp_data212,temp_data213,temp_data214,temp_data215


temp_data31 = [[] for _ in range(11)] #F up-kh
temp_data32 = [[] for _ in range(11)] #F in-kh
temp_data33 = [[] for _ in range(11)] #G up-khh
temp_data34 = [[] for _ in range(11)] #G in-khh
temp_data35 = [[] for _ in range(11)] #R-k
temp_data36 = [[] for _ in range(11)] #T-k
temp_data37 = [[] for _ in range(11)] #F up-kh1
temp_data38 = [[] for _ in range(11)] #F in-kh1
temp_data39 = [[] for _ in range(11)] #G up-kh1
temp_data310 = [[] for _ in range(11)] #G in-kh1
temp_data311 = [[] for _ in range(11)] 
temp_data312 = [[] for _ in range(11)] 
temp_data313 = [[] for _ in range(11)] 
temp_data314 = [[] for _ in range(11)] 



for k in range(1, 11):
    file_path = f"{direcElectron}min{k}ExtendedOmega.fits"
    with fits.open(file_path) as hdu:
        data31 = hdu[h_index].data.field('F_points_in')
        data32 = hdu[h_index].data.field('G_points_in')
        data33 = hdu[h_index].data.field('F_points_up')
        data34 = hdu[h_index].data.field('G_points_up')
        data35 = hdu[h_index].header['R']
        data36 = hdu[h_index].header['T']
        data313 = hdu[h_index].header['delta']

        data37 = hdu[h_prime_index].data.field('F_points_in')
        data38 = hdu[h_prime_index].data.field('G_points_in')
        data39 = hdu[h_prime_index].data.field('F_points_up')
        data310 = hdu[h_prime_index].data.field('G_points_up')
        data311 = hdu[h_prime_index].header['R']
        data312 = hdu[h_prime_index].header['T']
        data314 = hdu[h_prime_index].header['delta']

        temp_data31[k].append(data31)
        temp_data32[k].append(data32)
        temp_data33[k].append(data33)
        temp_data34[k].append(data34)
        temp_data35[k].append(data35)
        temp_data36[k].append(data36)
        temp_data37[k].append(data37)
        temp_data38[k].append(data38)
        temp_data39[k].append(data39)
        temp_data310[k].append(data310)
        temp_data311[k].append(data311)
        temp_data312[k].append(data312)
        temp_data313[k].append(data313)
        temp_data314[k].append(data314)

F_points_xminkh_in = np.array([np.array(data, dtype=complex) for data in temp_data31 if len(data) > 0], dtype=object)
G_points_xminkh_in = np.array([np.array(data, dtype=complex) for data in temp_data32 if len(data) > 0], dtype=object)
F_points_xminkh_up = np.array([np.array(data, dtype=complex) for data in temp_data33 if len(data) > 0], dtype=object)
G_points_xminkh_up = np.array([np.array(data, dtype=complex) for data in temp_data34 if len(data) > 0], dtype=object)
Rmink = np.array([np.array(data, dtype=complex) for data in temp_data35 if len(data) > 0], dtype=object).squeeze(axis=1)
Tmink = np.array([np.array(data, dtype=complex) for data in temp_data36 if len(data) > 0], dtype=object).squeeze(axis=1)
Dmink = np.array([np.array(data, dtype=complex) for data in temp_data313 if len(data) > 0], dtype=object).squeeze(axis=1)

F_points_xminkh_prime_in = np.array([np.array(data, dtype=complex) for data in temp_data37 if len(data) > 0], dtype=object)
G_points_xminkh_prime_in = np.array([np.array(data, dtype=complex) for data in temp_data38 if len(data) > 0], dtype=object)
F_points_xminkh_prime_up = np.array([np.array(data, dtype=complex) for data in temp_data39 if len(data) > 0], dtype=object)
G_points_xminkh_prime_up = np.array([np.array(data, dtype=complex) for data in temp_data310 if len(data) > 0], dtype=object)
Rmink_prime = np.array([np.array(data, dtype=complex) for data in temp_data311 if len(data) > 0], dtype=object).squeeze(axis=1)
Tmink_prime = np.array([np.array(data, dtype=complex) for data in temp_data312 if len(data) > 0], dtype=object).squeeze(axis=1)
Dmink_prime = np.array([np.array(data, dtype=complex) for data in temp_data314 if len(data) > 0], dtype=object).squeeze(axis=1)


del temp_data31, temp_data32, temp_data33, temp_data34,temp_data35,temp_data36,temp_data37,temp_data38,temp_data39,temp_data310,temp_data311,temp_data312,temp_data313,temp_data314
###################################################

def FGsort(k,ksign,FGk,FGmink):
    if (k*ksign>0):
        return np.squeeze(FGk)
    else:
        return np.squeeze(FGmink)

def RTsort(k,RTk,RTmink):
    if (k>0):
        return np.array(RTk[k-1]).reshape(1).item()
    else:
        return np.array(RTmink[-k-1]).reshape(1).item()


def safe_phase(z: complex) -> complex:
    """Return e^{i arg z} even if |z| = 0 (then define it as 1)."""
    mag = abs(z)
    return 1.0 if mag == 0 else z / mag   


def analytic_tail(h, hp, Rh: complex, Th: complex,Dh: complex,Rh_prime: complex, Th_prime: complex,Dh_prime: complex, X, Xp, r_match=-70.0,eps: float = 1e-5*T,e_mass = mu):
    alpha = h  - hp           # h − h′
    beta  = -alpha
    e_plus  = np.exp((eps + 1j*beta)  * r_match)
    e_minus = np.exp((eps + 1j*alpha) * r_match)

    eT   = safe_phase(Th)          # e^{i arg Th}
    eTp  = safe_phase(Th_prime)    # e^{i arg Th'}

    R_up = Rh
    R_prime_up =  Rh_prime
    if (h<e_mass):
        eT = 1j*np.exp(1j*Dh)
        R_up = 1

    if (hp<e_mass):
        eTp = 1j*np.exp(1j*Dh_prime)
        R_prime_up = 1

    if   (X, Xp) == (1,1):
        term1 =  1j / (alpha  + 1j*eps) * e_plus
        term2 =  -R_up*np.conjugate(R_prime_up)*(eTp**2) / (eT**2)* 1j / (alpha - 1j*eps) * e_minus
        return  term1 + term2
    elif (X, Xp) == (1,0):
        return  (2*R_up/(eT**2)*Th_prime* 1j / (alpha - 1j*eps)*e_minus)
    
    elif (X, Xp) == (0,1):
        return (2*np.conjugate(Th)*np.conjugate(R_prime_up)*(eTp**2)*1j / (alpha - 1j*eps)*e_minus)
    
    else:
        return (-2*np.conjugate(Th)*Th_prime*1j / (alpha - 1j*eps) * e_minus)



def make_Is():     
    array_Is = np.zeros((12, len(ks)), dtype=np.complex128)
    for k in range(len(ks)):

        Jobj = Inp.IfunctionsNoM(ks[k], h[x], h_prime[y], M)
        I1_upup_n = Jobj.Ione(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,F_points_xkh_prime_up[abs(ks[k])-1],F_points_xminkh_prime_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_prime_up[abs(ks[k])-1],G_points_xminkh_prime_up[abs(ks[k])-1]))
        I1_upin_n = Jobj.Ione(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,F_points_xkh_prime_in[abs(ks[k])-1],F_points_xminkh_prime_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_prime_in[abs(ks[k])-1],G_points_xminkh_prime_in[abs(ks[k])-1]))
        I1_inup_n = Jobj.Ione(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,F_points_xkh_prime_up[abs(ks[k])-1],F_points_xminkh_prime_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_prime_up[abs(ks[k])-1],G_points_xminkh_prime_up[abs(ks[k])-1]))
        I1_inin_n = Jobj.Ione(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,F_points_xkh_prime_in[abs(ks[k])-1],F_points_xminkh_prime_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_prime_in[abs(ks[k])-1],G_points_xminkh_prime_in[abs(ks[k])-1]))

        I2_upup_n = Jobj.Itwo(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],-1,F_points_xkh_prime_up[abs(ks[k])-1],F_points_xminkh_prime_up[abs(ks[k])-1]),FGsort(ks[k],-1,G_points_xkh_prime_up[abs(ks[k])-1],G_points_xminkh_prime_up[abs(ks[k])-1]))
        I2_upin_n = Jobj.Itwo(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_up[abs(ks[k])-1],F_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_up[abs(ks[k])-1],G_points_xminkh_up[abs(ks[k])-1]),FGsort(ks[k],-1,F_points_xkh_prime_in[abs(ks[k])-1],F_points_xminkh_prime_in[abs(ks[k])-1]),FGsort(ks[k],-1,G_points_xkh_prime_in[abs(ks[k])-1],G_points_xminkh_prime_in[abs(ks[k])-1]))
        I2_inup_n = Jobj.Itwo(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],-1,F_points_xkh_prime_up[abs(ks[k])-1],F_points_xminkh_prime_up[abs(ks[k])-1]),FGsort(ks[k],-1,G_points_xkh_prime_up[abs(ks[k])-1],G_points_xminkh_prime_up[abs(ks[k])-1]))
        I2_inin_n = Jobj.Itwo(rs,r_points_gamma,FGsort(ks[k],1,F_points_xkh_in[abs(ks[k])-1],F_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],1,G_points_xkh_in[abs(ks[k])-1],G_points_xminkh_in[abs(ks[k])-1]),FGsort(ks[k],-1,F_points_xkh_prime_in[abs(ks[k])-1],F_points_xminkh_prime_in[abs(ks[k])-1]),FGsort(ks[k],-1,G_points_xkh_prime_in[abs(ks[k])-1],G_points_xminkh_prime_in[abs(ks[k])-1]))

        I1_upup_tail = analytic_tail (h[x],h_prime[y],RTsort(ks[k],Rk,Rmink), RTsort(ks[k],Tk,Tmink),RTsort(ks[k],Dk,Dmink),RTsort(ks[k],Rk_prime,Rmink_prime), RTsort(ks[k],Tk_prime,Tmink_prime),RTsort(ks[k],Dk_prime,Dmink_prime), 1, 1)
        I1_upin_tail = analytic_tail (h[x],h_prime[y],RTsort(ks[k],Rk,Rmink), RTsort(ks[k],Tk,Tmink),RTsort(ks[k],Dk,Dmink),RTsort(ks[k],Rk_prime,Rmink_prime), RTsort(ks[k],Tk_prime,Tmink_prime),RTsort(ks[k],Dk_prime,Dmink_prime), 1, 0)
        I1_inup_tail = analytic_tail (h[x],h_prime[y],RTsort(ks[k],Rk,Rmink), RTsort(ks[k],Tk,Tmink),RTsort(ks[k],Dk,Dmink),RTsort(ks[k],Rk_prime,Rmink_prime), RTsort(ks[k],Tk_prime,Tmink_prime),RTsort(ks[k],Dk_prime,Dmink_prime), 0, 1)
        I1_inin_tail = analytic_tail (h[x],h_prime[y],RTsort(ks[k],Rk,Rmink), RTsort(ks[k],Tk,Tmink),RTsort(ks[k],Dk,Dmink),RTsort(ks[k],Rk_prime,Rmink_prime), RTsort(ks[k],Tk_prime,Tmink_prime),RTsort(ks[k],Dk_prime,Dmink_prime), 0, 0)

        array_Is[0,k] = I1_upup_n + I1_upup_tail
        array_Is[1,k] = I1_upin_n + I1_upin_tail
        array_Is[2,k] = I1_inup_n + I1_inup_tail
        array_Is[3,k] = I1_inin_n + I1_inin_tail
        array_Is[4,k] = I2_upup_n
        array_Is[5,k] = I2_upin_n
        array_Is[6,k] = I2_inup_n
        array_Is[7,k] = I2_inin_n

        array_Is[8,k] = I1_upup_n
        array_Is[9,k] = I1_upin_n
        array_Is[10,k] = I1_inup_n
        array_Is[11,k] = I1_inin_n    
    return array_Is



# ── build complex-data cube and write as primary image ─────────────
array_zeros = make_Is()                          # (8, 20) complex128

cube = np.empty(array_zeros.shape + (2,), dtype=np.float64)  # (8,20,2)
cube[..., 0] = array_zeros.real    # real part
cube[..., 1] = array_zeros.imag    # imag part

primary = fits.PrimaryHDU(data=cube)
primary.header['xidx'] = x          # fixed-h index
primary.header['yidx'] = y          # ω index

primary.writeto(f'output_x{x}_omega{y}.fits', overwrite=True)

del F_points_xkh_in,G_points_xkh_in,F_points_xkh_up,G_points_xkh_up,Rk,Tk,Dk,Rk_prime,Tk_prime,Dk_prime
del F_points_xminkh_in,G_points_xminkh_in,F_points_xminkh_up,G_points_xminkh_up,Rmink,Tmink,Dmink,Rmink_prime,Tmink_prime,Dmink_prime
gc.collect()