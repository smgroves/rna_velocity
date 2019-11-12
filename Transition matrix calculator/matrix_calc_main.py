# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:36:09 2018

@author: John
"""

# Timer
import time as ti
start_time = ti.time()


import sys
sys.path.insert(1, 'model_files')
sys.path.insert(2, 'simulation_files')

from get_fake_data import generate_fake_data
from get_t_matrix import get_state_list, get_t_matrix
#from IC_generator import generate_IC


# ============================================================================================

calc_t_matrix = 'no'    #options: 'yes' and 'no'

# =============================================================================
# tol = 0.1
# step_size = 0.001
# stop_tol = 0.000001
# num_attractors, protein_IC = generate_IC(tol, step_size, stop_tol)
# =============================================================================
#print(protein_IC)

N_cells = 20
num_timepoints = 2                   # number of time points, INCLUDING initial time point
time_between_measurements =  1        # currently arbitrary units


# Generate fake data
print("Generating fake data for ", N_cells, " cells with ", num_timepoints, " time points per cell...")
u, s, step_size, timesteps_per_measurement = generate_fake_data(N_cells, num_timepoints, time_between_measurements)
print("Fake data successfully generated!")


# Generate 'states'.
state_list, num_states, state_boundaries = get_state_list(s)


if calc_t_matrix == 'yes':
    print("Calculating transition matrix...")
    
    
    # Compute P_trans starting at each state for fixed time. 
    t_matrix = get_t_matrix(state_list, num_states, state_boundaries)




    print("Final transition matrix")
    print(t_matrix)
    print("T matrix size", t_matrix.size)


# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))



