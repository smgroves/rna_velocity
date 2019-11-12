# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 19:21:02 2019

@author: johnv
"""
import numpy as np
from model_parameters import beta, kp, dp, ds
params = []

p0 = 11.80201080612191         # coordinates of two attractors
#p0 = 46.858926564535096

s0 = (kp/dp)*p0
u0 = (ds/beta)*s0

#u0 = 10
#s0 = 10
#p0 = 10
x0 = [u0, s0, p0]

t0 = 0
tf = 80

step_size_naive = 0.1
num_timesteps = int(np.ceil(tf/step_size_naive))
step_size = tf/num_timesteps

dt_splice = 1/beta
timesteps_per_splice = 100
step_size_tmatrix = dt_splice/timesteps_per_splice
    
t = np.linspace(t0, tf, num_timesteps + 1)

N_runs_ptrans = 100