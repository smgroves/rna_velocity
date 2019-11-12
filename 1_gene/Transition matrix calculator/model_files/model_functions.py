# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:17:00 2018

@author: John
"""

# This file contains functions related to the model's SDE, including the functions governing the
# deterministic and noise dynamics. 

# These functions depend on the parameters found in eqn_parameters.py . 

import numpy as np
from model_parameters import kp, G, ds, dp, N, k, c, b_r, beta



# ==============================================================
    

# Hill function from protein SDE related to protein translation.    
def hill(x):
    
    hill_top = 0
    hill_bot = 0
    for i in range(0, N + 1):
        hill_top = hill_top + k[i]*c[i]*(x**i)
        hill_bot = hill_bot + c[i]*(x**i)
    
    hill_func = (hill_top/hill_bot) 

    
    result = hill_func 
    return result


# Hill function related to gene-protein binding from the noise part of the protein SDE.
def hill_bind(x):
    hill_top = 0
    hill_bot = 0
    for i in range(1,N+1):
        hill_top = hill_top + b_r[i]*c[i]*(x**i)
    
    for i in range(0,N+1):
        hill_bot = hill_bot + c[i]*(x**i)
        
    hill_func = 2*G*(hill_top/hill_bot)
    
    result = hill_func
    return result

# ==================================================

# Deterministic part of unspliced SDE.
def f_u(x, params):
    u = x[0]
    p = x[2]
    
    result = G*hill(p) - beta*u
    return result


# Deterministic part of spliced SDE.
def f_s(x, params):
    u = x[0]
    s = x[1]
    
    result = beta*u - ds*s
    return result


# Deterministic part of protein SDE.
def f_p(x, params):
    s = x[1]
    p = x[2]
    
    result = kp*s - dp*p
    return result

# ==================================================

# Noise magnitude for transcription reactions. 
def g_transcription(x, params):
    p = x[2]
    
    result = np.sqrt(G*hill(p))
    return result

# Noise magnitude for splicing reaction.
def g_splice(x, params):
    u = x[0]
    
    result = np.sqrt(beta*u)
    return result

# Noise magnitude for spliced degradation.
def g_sdeg(x, params):
    s = x[1]
    
    result = np.sqrt(ds*s)
    return result

# Noise magnitude for protein SDE.
def g_protein(x, params):
    s = x[1]
    p = x[2]
    
    result = np.sqrt( kp*s + dp*p + hill_bind(p) )
    return result

# ====================================================================
 


# Additive noise functions
    
def g_add(x, params):
    # Additive noise
    return params[0]


def g_prime_add(x, params):
    return 0

def g_add_dummy(x, params):
    return 1






# Multiplicative noise functions
    
def g_mult(x, params):
    return params[1]*x + params[2]

def g_prime_mult(x, params):
    return params[1]

def g_mult_dummy(x, params):
    return x