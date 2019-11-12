# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 19:33:56 2019

@author: johnv
"""

import numpy as np

from get_probdist import one_traj
from simulation_parameters import x0, params, step_size_naive

def generate_fake_data(N_cells, num_timepoints, time_between_measurements):
    
    # Determine number of time steps and step size.
    #total_time = (num_timepoints - 1)*time_between_measurements
    timesteps_per_measurement = int(np.ceil(time_between_measurements/step_size_naive))
    step_size = time_between_measurements/timesteps_per_measurement
    
    num_timesteps = timesteps_per_measurement*(num_timepoints - 1)
    
    
    
    u_measured = np.zeros((N_cells, num_timepoints))
    s_measured = np.zeros((N_cells, num_timepoints))
    
    for i in range(0, N_cells):
        x = one_traj(x0, num_timesteps, step_size, params)
        for j in range(0, num_timesteps + 1):
            if j%timesteps_per_measurement == 0:
                my_j = int(j/timesteps_per_measurement)
                u_measured[i,my_j] = x[0,j]
                s_measured[i,my_j] = x[1,j]
    

    return u_measured, s_measured, step_size, timesteps_per_measurement