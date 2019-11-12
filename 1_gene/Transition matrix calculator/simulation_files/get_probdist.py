# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 09:37:09 2019

@author: John
"""
import numpy as np
from scipy.stats import kde
from model_functions import f_u, f_s, f_p
from model_functions import g_transcription, g_splice, g_sdeg, g_protein


#def euler_maruyama(x, f, g, step_size, params):
#    #D = len(x)
#    D = 3
#    x_new = -100*np.ones(D)
#    
#    for i in range(0, D):      
#        while x_new[i] < 0:
#            r = np.random.normal(0,1)
#            x_new[i] = x[i] + f[i](x, params)*step_size + g[i](x, params)*np.sqrt(step_size)*r
#   
#    return x_new  


def u_step(x, step_size, params):
    u = x[0]
      
    u_new = -1
    while u_new < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_transcription(x, params)*r1 - g_splice(x, params)*np.sqrt(step_size)*r2)
        u_new = u + f_u(x, params)*step_size + noise
   
    return u_new  


def s_step(x, step_size, params):
    s = x[1]   

    s_new = -1
    while s_new < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_splice(x, params)*r1 + g_sdeg(x, params)*np.sqrt(step_size)*r2)
        s_new = s + f_s(x, params)*step_size + noise
   
    return s_new  


def p_step(x, step_size, params):
    p = x[2]
      
    p_new = -1
    while p_new < 0:
        r = np.random.normal(0,1)
        noise = np.sqrt(step_size)*g_protein(x, params)*r
        p_new = p + f_p(x, params)*step_size + noise
   
    return p_new  


def euler_step(x, step_size, params):
    u_new = u_step(x, step_size, params)
    s_new = s_step(x, step_size, params)
    p_new = p_step(x, step_size, params)
    
    return [u_new, s_new, p_new]

#=====================================================

    
def one_traj(x0, num_timesteps, step_size, params):
    D = 3
    xf = np.zeros((D, num_timesteps + 1))
    
    xf[:,0] = x0
    for i in range(1, num_timesteps + 1):
        xf[:,i] = euler_step(xf[:,i-1], step_size, params) 
    return xf


def N_timesteps(x0, num_timesteps, step_size, params):
    x_next = x0
    for i in range(0, num_timesteps):
        x_next = euler_step(x_next, step_size, params) 
    
    return x_next


def N_traj(N, x0, num_timesteps, step_size, params):
    # Initialize xf arrays.
    D = 3
    xf = np.zeros((D, N))
    
    # Do a bunch of simulations.
    for j in range(0, N):
        xf[:,j] = N_timesteps(x0, num_timesteps, step_size, params)
            
    return xf    


def get_prob_dist(N, x0, num_timesteps, step_size, params, num_KDE_pts): 
    # Get N trajectories.
    xf = N_traj(N, x0, num_timesteps, step_size, params)
    
    # Get bounds
    #xmin, xmax = my_min[0], my_max[0]
    #ymin, ymax = my_min[1], my_max[1]
    
    # Generate KDE from xf points.      
#    X, Y = np.mgrid[xmin:xmax:np.complex(0, num_KDE_pts), ymin:ymax:np.complex(0, num_KDE_pts)]
#    positions = np.vstack([X.ravel(), Y.ravel()])
#    values = np.vstack([xf[0,:], xf[1,:]])             
#    kernel = stats.gaussian_kde(values)
#    KDE = np.reshape(kernel(positions).T, X.shape)
    
    #k = kde.gaussian_kde(data.T)
    my_min =[xf[0,:].min(), xf[1,:].min()]
    my_max = [xf[0,:].max(), xf[1,:].max()]
    
    X, Y = np.mgrid[my_min[0]:my_max[0]:num_KDE_pts*1j, my_min[1]:my_max[1]:num_KDE_pts*1j]
    values = np.vstack([xf[0,:], xf[1,:]])             
    k = kde.gaussian_kde(values)
    KDE = k(np.vstack([X.flatten(), Y.flatten()]))
    


    
    return KDE, X, Y, my_min, my_max, xf