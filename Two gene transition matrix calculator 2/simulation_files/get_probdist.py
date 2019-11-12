# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 09:37:09 2019

@author: John
"""
import numpy as np
from scipy.stats import kde
from model_functions import f_u_x, f_s_x, f_p_x
from model_functions import f_u_y, f_s_y, f_p_y
from model_functions import g_tran_x, g_splice_x, g_sdeg_x, g_protein_x
from model_functions import g_tran_y, g_splice_y, g_sdeg_y, g_protein_y


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


def ux_step(x, y, step_size, params):
    u_x = x[0]
      
    u_new_x = -1
    while u_new_x < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_tran_x(x, y, params)*r1 - g_splice_x(x, y, params)*np.sqrt(step_size)*r2)
        u_new_x = u_x + f_u_x(x, y, params)*step_size + noise
   
    return u_new_x  

def uy_step(x, y, step_size, params):
    u_y = y[0]
      
    u_new_y = -1
    while u_new_y < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_tran_y(x, y, params)*r1 - g_splice_y(x, y, params)*np.sqrt(step_size)*r2)
        u_new_y = u_y + f_u_y(x, y, params)*step_size + noise
   
    return u_new_y  


def sx_step(x, y, step_size, params):
    s_x = x[1]   

    s_new_x = -1
    while s_new_x < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_splice_x(x, y, params)*r1 + g_sdeg_x(x, y, params)*np.sqrt(step_size)*r2)
        s_new_x = s_x + f_s_x(x, y, params)*step_size + noise
   
    return s_new_x  

def sy_step(x, y, step_size, params):
    s_y = y[1]   

    s_new_y = -1
    while s_new_y < 0:
        r1 = np.random.normal(0,1)
        r2 = np.random.normal(0,1)
        noise = np.sqrt(step_size)*(g_splice_y(x, y, params)*r1 + g_sdeg_y(x, y, params)*np.sqrt(step_size)*r2)
        s_new_y = s_y + f_s_y(x, y, params)*step_size + noise
   
    return s_new_y  


def px_step(x, y, step_size, params):
    p_x = x[2]
      
    p_new_x = -1
    while p_new_x < 0:
        r = np.random.normal(0,1)
        noise = np.sqrt(step_size)*g_protein_x(x, y, params)*r
        p_new_x = p_x + f_p_x(x, y, params)*step_size + noise
   
    return p_new_x  

def py_step(x, y, step_size, params):
    p_y = y[2]
      
    p_new_y = -1
    while p_new_y < 0:
        r = np.random.normal(0,1)
        noise = np.sqrt(step_size)*g_protein_y(x, y, params)*r
        p_new_y = p_y + f_p_y(x, y, params)*step_size + noise
   
    return p_new_y  


def euler_step(x, y, step_size, params):
    u_new_x = ux_step(x, y, step_size, params)
    s_new_x = sx_step(x, y, step_size, params)
    p_new_x = px_step(x, y, step_size, params)

    u_new_y = uy_step(x, y, step_size, params)
    s_new_y = sy_step(x, y, step_size, params)
    p_new_y = py_step(x, y, step_size, params)
    
    x_new = [u_new_x, s_new_x, p_new_x]
    y_new = [u_new_y, s_new_y, p_new_y]
    
    return x_new, y_new

#=====================================================

    
def one_traj(x0, y0, num_timesteps, step_size, params):
    D = 3
    xf = np.zeros((D, num_timesteps + 1))
    yf = np.zeros((D, num_timesteps + 1))
    
    xf[:,0] = x0
    yf[:,0] = y0
    for i in range(1, num_timesteps + 1):
        xf[:,i], yf[:,i] = euler_step(xf[:,i-1], yf[:,i-1], step_size, params) 
    return xf, yf


def N_timesteps(x0, y0, num_timesteps, step_size, params):
    x_next = x0
    y_next = y0
    for i in range(0, num_timesteps):
        x_next, y_next = euler_step(x_next, y_next, step_size, params) 
    
    return x_next, y_next


def N_traj(N, x0, y0, num_timesteps, step_size, params):
    # Initialize xf arrays.
    D = 3
    xf = np.zeros((D, N))
    yf = np.zeros((D, N))
    
    # Do a bunch of simulations.
    for j in range(0, N):
        xf[:,j], yf[:,j] = N_timesteps(x0, y0, num_timesteps, step_size, params)
            
    return xf, yf    


#def get_prob_dist(N, x0, y0, num_timesteps, step_size, params, num_KDE_pts): 
#    # Get N trajectories.
#    xf, yf = N_traj(N, x0, y0, num_timesteps, step_size, params)
#    
#    # Get bounds
#    #xmin, xmax = my_min[0], my_max[0]
#    #ymin, ymax = my_min[1], my_max[1]
#    
#    # Generate KDE from xf points.      
##    X, Y = np.mgrid[xmin:xmax:np.complex(0, num_KDE_pts), ymin:ymax:np.complex(0, num_KDE_pts)]
##    positions = np.vstack([X.ravel(), Y.ravel()])
##    values = np.vstack([xf[0,:], xf[1,:]])             
##    kernel = stats.gaussian_kde(values)
##    KDE = np.reshape(kernel(positions).T, X.shape)
#    
#    #k = kde.gaussian_kde(data.T)
#    my_min =[xf[0,:].min(), xf[1,:].min()]
#    my_max = [xf[0,:].max(), xf[1,:].max()]
#    
#    X, Y = np.mgrid[my_min[0]:my_max[0]:num_KDE_pts*1j, my_min[1]:my_max[1]:num_KDE_pts*1j]
#    values = np.vstack([xf[0,:], xf[1,:]])             
#    k = kde.gaussian_kde(values)
#    KDE = k(np.vstack([X.flatten(), Y.flatten()]))
#    
#
#
#    
#    return KDE, X, Y, my_min, my_max, xf