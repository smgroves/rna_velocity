# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:17:00 2018

@author: John
"""

# This file contains functions related to the model's SDE, including the functions governing the
# deterministic and noise dynamics. 

# These functions depend on the parameters found in eqn_parameters.py . 

import numpy as np
from model_parameters import G, N, c

from model_parameters import b_x, k_A_x, k_I_x
from model_parameters import kp_x, ds_x, dp_x, beta_x

from model_parameters import b_y, k_A_y, k_I_y
from model_parameters import kp_y, ds_y, dp_y, beta_y

from model_parameters import sig_x, sig_y



# ==============================================================
    

# Hill function from protein SDE related to protein translation.    
#def hill(x):
#    
#    hill_top = 0
#    hill_bot = 0
#    for i in range(0, N + 1):
#        hill_top = hill_top + k[i]*c[i]*(x**i)
#        hill_bot = hill_bot + c[i]*(x**i)
#    
#    hill_func = (hill_top/hill_bot) 
#
#    
#    result = hill_func 
#    return result



def hill_plain(x, const):
    top = x**N
    bot = const + x**N
    result = top/bot
    
    return result


# Hill function related to gene-protein binding from the noise part of the protein SDE.
#def hill_bind(x):
#    hill_top = 0
#    hill_bot = 0
#    for i in range(1,N+1):
#        hill_top = hill_top + b_r[i]*c[i]*(x**i)
#    
#    for i in range(0,N+1):
#        hill_bot = hill_bot + c[i]*(x**i)
#        
#    hill_func = 2*G*(hill_top/hill_bot)
#    
#    result = hill_func
#    return result

# ==================================================

# Deterministic part of unspliced SDE.
def f_u_x(x, y, params):
    u_x = x[0]
    p_x = x[2]
    p_y = y[2]
    
    result = b_x  - k_I_x*G*hill_plain(p_y, c**N) + k_A_x*G*hill_plain(p_x, c**N) - beta_x*u_x
    return result

def f_u_y(x, y, params):
    u_y = y[0]
    p_y = y[2]
    p_x = x[2]
    
    result = b_y - k_I_y*G*hill_plain(p_x, c**N) + k_A_y*G*hill_plain(p_y, c**N) - beta_y*u_y
    return result


# Deterministic part of spliced SDE.
def f_s_x(x, y, params):
    u_x = x[0]
    s_x = x[1]
    
    result = beta_x*u_x - ds_x*s_x
    return result

def f_s_y(x, y, params):
    u_y = y[0]
    s_y = y[1]
    
    result = beta_y*u_y - ds_y*s_y
    return result


# Deterministic part of protein SDE.
def f_p_x(x, y, params):
    s_x = x[1]
    p_x = x[2]
    
    result = kp_x*s_x - dp_x*p_x
    return result

def f_p_y(x, y, params):
    s_y = y[1]
    p_y = y[2]
    
    result = kp_y*s_y - dp_y*p_y
    return result

# ==================================================

# Noise magnitude for transcription reactions. 
def g_tran_x(x, y, params):
    p_x = x[2]
    p_y = y[2]  
    
    result = np.sqrt( b_x - k_I_x*G*hill_plain(p_y, c**N) + k_A_x*G*hill_plain(p_x, c**N) )
    return result

def g_tran_y(x, y, params):
    p_x = x[2]
    p_y = y[2]  
    
    result = np.sqrt( b_y - k_I_y*G*hill_plain(p_x, c**N) + k_A_y*G*hill_plain(p_y, c**N) )
    return result

# Noise magnitude for splicing reaction.
def g_splice_x(x, y, params):
    u_x = x[0]
    
    result = np.sqrt(beta_x*u_x)
    return result

def g_splice_y(x, y, params):
    u_y = y[0]
    
    result = np.sqrt(beta_y*u_y)
    return result


# Noise magnitude for spliced degradation.
def g_sdeg_x(x, y, params):
    s_x = x[1]
    
    result = np.sqrt(ds_x*s_x)
    return result

def g_sdeg_y(x, y, params):
    s_y = y[1]
    
    result = np.sqrt(ds_y*s_y)
    return result

# Noise magnitude for protein SDE.
def g_protein_x(x, y, params):
    s_x = x[1]
    p_x = x[2]
    
    result = np.sqrt( kp_x*s_x + dp_x*p_x + sig_x)
    return result

def g_protein_y(x, y, params):
    s_y = y[1]
    p_y = y[2]
    
    result = np.sqrt( kp_y*s_y + dp_y*p_y + sig_y)
    return result

# ====================================================================
 


# Additive noise functions
    
#def g_add(x, params):
#    # Additive noise
#    return params[0]
#
#
#def g_prime_add(x, params):
#    return 0
#
#def g_add_dummy(x, params):
#    return 1
#
#
#
#
#
#
## Multiplicative noise functions
#    
#def g_mult(x, params):
#    return params[1]*x + params[2]
#
#def g_prime_mult(x, params):
#    return params[1]
#
#def g_mult_dummy(x, params):
#    return x