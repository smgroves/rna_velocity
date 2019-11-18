# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:28:25 2019

@author: johnv
"""

import numpy as np
from simulation_parameters import N_runs_ptrans
from model_parameters import kp_x, dp_x, ds_x, beta_x
from model_parameters import kp_y, dp_y, ds_y, beta_y
from get_probdist import N_traj


from scipy.spatial import cKDTree


def get_state_list(s_x, s_y):
    #print(s)
    sx_flat = s_x.flatten()                          # Flatten.
    sy_flat = s_y.flatten() 
    
    s_size = sx_flat.size
    s_together = np.zeros((s_size, 2))
    s_together[:,0], s_together[:,1] = sx_flat, sy_flat
    
    print("s together shape ", s_together.shape)
    state_list = np.unique(s_together, axis=0)             # Remove duplicates
    #print("uniques shape ", s_uniques.shape)
    #state_list = np.sort(s_uniques)          # Sort from low to high to get final list of states 
    print(state_list)
    print(state_list.shape)

    num_states = int(state_list.size/2)
    print("Num states is ", num_states)
    
    voronoi_kdtree = cKDTree(state_list)
    
    return state_list, num_states, voronoi_kdtree



def get_p_trans(cur_state, num_states, voronoi_kdtree, step_size, timesteps_per_measurement):
    p_trans = np.zeros(num_states)
    
    s_x, s_y = cur_state
    u_x, u_y = (ds_x/beta_x)*s_x , (ds_y/beta_y)*s_y
    p_x, p_y = (kp_x/dp_x)*s_x , (kp_y/dp_y)*s_y
    
    x0 = [u_x, s_x, p_x]
    y0 = [u_y, s_y, p_y]
        
    params = []
    
    xf, yf = N_traj(N_runs_ptrans, x0, y0, timesteps_per_measurement, step_size, params)
    
    s_runs = np.zeros((N_runs_ptrans, 2))
    s_runs[:,0], s_runs[:,1] = xf[0,:], xf[1,:]
    
    
    test_point_dist, test_point_regions = voronoi_kdtree.query(s_runs, k=1)

    # test_point_regions Now holds an array of shape (n_test, 1) with the indices 
    # of the points in voronoi_points closest to each of your test_points.

    counts = np.zeros(num_states)

    for i in range(0, N_runs_ptrans):
        my_index = test_point_regions[i]
        counts[my_index] = counts[my_index] + 1

        
    p_trans = counts/N_runs_ptrans
    

    return p_trans


def get_t_matrix(state_list, num_states, voronoi_kdtree, step_size, timesteps_per_measurement):
    
    t_matrix = np.zeros((num_states, num_states))
    
    # For each state in state list, compute P_trans to every other state.
    for i in range(0, num_states):
        cur_state = state_list[i,:]
        t_matrix[i,:] = get_p_trans(cur_state, num_states, voronoi_kdtree, step_size, timesteps_per_measurement)
    
    
    
    
    
    return t_matrix