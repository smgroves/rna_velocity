# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:28:25 2019

@author: johnv
"""

import numpy as np
from simulation_parameters import N_runs_ptrans, timesteps_per_splice, step_size_tmatrix
from model_parameters import kp, dp, ds, beta
from get_probdist import N_traj




def get_state_list(s):
    #print(s)
    s_flat = s.flatten()                          # Flatten.
    s_flat_uniques = np.unique(s_flat)            # Remove duplicates
    state_list = np.sort(s_flat_uniques)          # Sort from low to high to get final list of states 
    #print(state_list)

    num_states = state_list.size
    #print("Num states is ", num_states)
    
    state_boundaries = np.zeros(num_states + 1)
    for i in range(1, num_states):
        state_boundaries[i] = (state_list[i] + state_list[i-1])/2

    state_boundaries[num_states] = 100000

    #print(state_boundaries)
    
    return state_list, num_states, state_boundaries



def get_p_trans(cur_state, num_states, state_boundaries):
    p_trans = np.zeros(num_states)
    
    s = cur_state
    u = (ds/beta)*s
    p = (kp/dp)*s
    x0 = [u, s, p]
        
    params = []
    
    xf = N_traj(N_runs_ptrans, x0, timesteps_per_splice, step_size_tmatrix, params)
    sf = np.sort(xf[1,:])
    
    j = 0
    for i in range(0, N_runs_ptrans):
        while (sf[i] > state_boundaries[j]):
            j = j + 1
        p_trans[j - 1] = p_trans[j - 1] + 1
        
    p_trans = p_trans/N_runs_ptrans
    

    return p_trans


def get_t_matrix(state_list, num_states, state_boundaries):
    
    t_matrix = np.zeros((num_states, num_states))
    
    # For each state in state list, compute P_trans to every other state.
    for i in range(0, num_states):
        cur_state = state_list[i]
        t_matrix[i,:] = get_p_trans(cur_state, num_states, state_boundaries)
    
    return t_matrix