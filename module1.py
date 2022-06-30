#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import math
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy.optimize import fsolve


# In[2]:


def r_to_r_star(r, M):
    return r + 2*M*np.log((r - 2*M)/(2*M))


def r_star_to_r_condition(r_star,r, M): 
    return  r_star - r_to_r_star(r, M) 


def r_star_to_r(r_star, M, tol): #will need to replace all the local variables using self 

    if r_star > 12*M:

        r_prev = (2+1e-10)*M

        r_next = r_star - 2*M * np.log((r_prev - 2*M)/(2*M))

        #Counter for the while loop
        i = 0

        while abs((r_next - r_prev)/r_prev) > tol and i < 50: 

            i = i + 1

            r_prev = r_next

            r_next = r_star - 2*M * np.log((r_prev - 2*M)/(2*M))


        return r_next


    if r_star < -3*M: 

        r_prev = 2*M

        r_next = (2*M) * (1 + math.e**((r_star - r_prev) / (2*M)))

        i = 0

        while abs((r_next - r_prev)/r_prev) > tol and i < 50: 

            i = i + 1

            r_prev = r_next

            r_next = 2*M*(1 + math.e**((r_star - r_prev) / (2*M)))

        return r_next 


    else: 
        #r initial needs to be whatever r is for -10M = r_star, and same for r_final 

        r_initial = 2.1*M
        r_final = 10*M

        sol = optimize.root_scalar(lambda r: r_star_to_r_condition(r_star, r, M), bracket=[r_initial, r_final], method = 'brentq')

        return sol.root 


# In[ ]:




