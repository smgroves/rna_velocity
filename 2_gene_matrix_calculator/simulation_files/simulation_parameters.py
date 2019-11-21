# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 19:21:02 2019

@author: johnv
"""
import numpy as np

from model_files.model_parameters import ds_x, dp_x, kp_x, beta_x
from model_files.model_parameters import ds_y, dp_y, kp_y, beta_y

from choose_model_here import p_x0, p_y0

params = []



u_x0, u_y0 = ((ds_x*dp_x)/(kp_x*beta_x))*p_x0, ((ds_y*dp_y)/(kp_y*beta_y))*p_y0
s_x0, s_y0 = (beta_x/ds_x)*u_x0, (beta_y/ds_y)*u_y0

#u0 = 10
#s0 = 10
#p0 = 10

#Default: starting at a steady state
# x0 = [u_x0, s_x0, p_x0]
# y0 = [u_y0, s_y0, p_y0]

# Updated: starting at a different state
x0 = [.1*u_x0, 2*s_x0, 0*p_x0]
y0 = [.1*u_y0, 2*s_y0, 0*p_y0]

t0 = 0
tf = 80

step_size_naive = 0.1
num_timesteps = int(np.ceil(tf/step_size_naive))
step_size = tf/num_timesteps

t = np.linspace(t0, tf, num_timesteps + 1)

N_runs_ptrans = 2000