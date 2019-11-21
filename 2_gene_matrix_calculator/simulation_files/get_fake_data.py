# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 19:33:56 2019

@author: johnv
"""

import numpy as np

from get_probdist import one_traj
from simulation_parameters import x0, y0, params, step_size_naive


def generate_fake_data(N_cells, num_timepoints, time_between_measurements, rand_start = False):
    # for random starting points
    if rand_start == True:
        x0 = np.random.rand(N_cells)
        y0 = np.random.rand(N_cells)

    # Determine number of time steps and step size.
    timesteps_per_measurement = int(np.ceil(time_between_measurements/step_size_naive))
    step_size = time_between_measurements/timesteps_per_measurement
    
    num_timesteps = timesteps_per_measurement*(num_timepoints - 1)
      
    
    ux_measured = np.zeros((N_cells, num_timepoints))
    sx_measured = np.zeros((N_cells, num_timepoints))

    uy_measured = np.zeros((N_cells, num_timepoints))
    sy_measured = np.zeros((N_cells, num_timepoints))
    
    for i in range(0, N_cells):
        if rand_start == True:
            x, y = one_traj(x0[i]*5, y0[i]*5, num_timesteps, step_size, params)
        else:
            x, y = one_traj(x0, y0, num_timesteps, step_size, params)
        for j in range(0, num_timesteps + 1):
            if j%timesteps_per_measurement == 0:
                my_j = int(j/timesteps_per_measurement)
                
                ux_measured[i,my_j] = x[0,j]
                sx_measured[i,my_j] = x[1,j]
                
                uy_measured[i,my_j] = y[0,j]
                sy_measured[i,my_j] = y[1,j]
    

    return ux_measured, uy_measured, sx_measured, sy_measured, step_size, timesteps_per_measurement