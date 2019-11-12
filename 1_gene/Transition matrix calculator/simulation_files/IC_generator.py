# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 18:42:26 2018

@author: John
"""
import numpy as np
from model_parameters import kp, G, ds, dp, N, k, c

#-------------------------

def hill(x):
    
    hill_top = 0
    hill_bot = 0
    for i in range(0, N + 1):
        hill_top = hill_top + k[i]*c[i]*(x**i)
        hill_bot = hill_bot + c[i]*(x**i)
    
    hill_func = (hill_top/hill_bot) 

    
    result = hill_func 
    return result

def f(p):
    
    result = G*hill(p) - ((ds*dp)/(kp))*p
    return result


# Classical Runge-Kutta RK4 solver.
def RK4(x_old, step_size):
    j1 = f(x_old)
    j2 = f(x_old + (step_size/2)*j1)  
    j3 = f(x_old + (step_size/2)*j2)
    j4 = f(x_old + (step_size)*j3)
    
    x_new = x_old + (step_size/6)*(j1 + 2*j2 + 2*j3 + j4)
     
    return x_new


# Distance metric. Trivial in 1D, but useful to define in higher dim.
def dist(x1, x2):
    t1 = (x1 - x2)**2
    
    dist = np.sqrt(t1)
    return dist


def get_IC(test_x, step_size, stop_tol):

    flag = 0
    while flag == 0:
        sol = RK4(test_x, step_size) # Runge-Kutta method
        
        if test_x == sol:
            flag = 1
        elif (dist(test_x, sol) < stop_tol):
            flag = 1
        else:
            test_x = sol
    x0 = test_x
    
    return x0
    


def generate_IC(tol, step_size, stop_tol):

    pss_lowerbound = 0.01
   # pss_lowerbound = (k[0]*kp*G)/(ds*dp)            # lower bound for p_ss
    pss_upperbound = (k[N]*kp*G)/(ds*dp)           # upper bound for p_ss
    
    pss_low = get_IC(pss_lowerbound, step_size, stop_tol)   # find closest steady state to lower bound
    pss_high = get_IC(pss_upperbound, step_size, stop_tol)  # find closest steady state to upper bound
    
    
    difference = np.abs(pss_low - pss_high)
    if difference < tol:
        # Only 1 steady state. Difference between pss_low and pss_high probably due to numerical error.
        num_steadystates = 1        
        protein_IC = [(pss_low + pss_high)/2]
    else:
        # Probably 2 steady states.
        num_steadystates = 2  
        protein_IC = [pss_low, pss_high]
        
        
    print(pss_low)
    print(pss_high)

    return num_steadystates, protein_IC