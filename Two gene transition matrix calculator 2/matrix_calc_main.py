# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:36:09 2018

@author: John
"""

# Timer
import time as ti
start_time = ti.time()

import numpy as np

import sys
sys.path.insert(1, 'model_files')
sys.path.insert(2, 'simulation_files')

from get_fake_data import generate_fake_data
from get_t_matrix import get_state_list, get_t_matrix
from network_utils import save_data
#from generate_IC import generate_IC
outdir= '/Users/sarahmaddox/Documents/workspace/rna_velocity/output/'
import random
# ============================================================================================
seed = random.randint(0, 1000000000)
calc_t_matrix = 'no'    #options: 'yes' and 'no'


#num_zeros, protein_IC = generate_IC()
#print(num_zeros)
#print(protein_IC)

N_cells = 20
num_timepoints = 2                   # number of time points, INCLUDING initial time point
time_between_measurements =  1        # currently arbitrary units


# Generate fake data
print("Generating fake data for ", N_cells, " cells with ", num_timepoints, " time points per cell...")
u_x, u_y, s_x, s_y, step_size, timesteps_per_measurement = generate_fake_data(N_cells, num_timepoints, time_between_measurements)
print("Fake data successfully generated!")

#print("s_x is ", s_x.flatten())
#print("s_y is ", s_y.flatten())

#Save fake data
save_data([u_x,u_y],[s_x,s_y],outdir,N_cells,num_timepoints, model = "two_gene",seed = seed)



# -----------------------------------------------------------

if calc_t_matrix == 'yes':
    print("Calculating transition matrix...")
    # Generate 'states'.
    state_list, num_states, voronoi_kdtree = get_state_list(s_x, s_y)
    
    
    # Compute P_trans starting at each state for fixed time. 
    t_matrix = get_t_matrix(state_list, num_states, voronoi_kdtree, step_size, timesteps_per_measurement)
    t_matrix2 = get_t_matrix(state_list, num_states, voronoi_kdtree, step_size, timesteps_per_measurement)
    
    
    
    print("Final transition matrix")
    print(t_matrix)
    print(t_matrix2)
    print("T matrix size", t_matrix.size)
    
    m_norm = np.linalg.norm(t_matrix - t_matrix2)
    print("Matrix norm :", m_norm)
    w, v = np.linalg.eig(t_matrix)
    print(w)
    
    w2, v2 = np.linalg.eig(t_matrix2)
    print(w2)


# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))



