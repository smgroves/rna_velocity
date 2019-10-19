# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:26:43 2019

@author: johnv
"""

from __future__ import print_function
from pysb.simulator.bng import BngSimulator
from firstorder_gene_us import model

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------


# Get u-s fit line.
beta = 2
d_s = 1

u_min = 0
u_max = 17
num_u = 100
u_line = np.linspace(u_min, u_max, num_u)
s_line = np.zeros(num_u)
for i in range(0, num_u):
    s_line[i] = (beta/d_s)*u_line[i]
# ----------------------------------------------------------   

# Define observed times.
tmin = 0
tmax = 50
num_t = 100
t = np.linspace(tmin, tmax, num_t)


# Simulate model.
simulator = BngSimulator(model, tspan=t)
simresult = simulator.run(method='ssa')  # ssa: discrete stochastic simulations
yout = simresult.all


# Plot results.
plt.plot(t, yout['o_g_0'], label="$g_0$")
plt.plot(t, yout['o_g_1'], label="$g_1$")
plt.plot(t, yout['o_u'], label='$u$')
plt.plot(t, yout['o_s'], label="$s$")
plt.plot(t, yout['o_p'], label="$p$")



plt.xlabel('$t$')
plt.ylabel('molecule number')
plt.title('Molecule number vs time')

plt.legend()
plt.show()


plt.scatter(yout['o_u'], yout['o_s'])
plt.plot(u_line, s_line, linestyle='--', linewidth=3, alpha=0.6, color='black')
plt.xlabel('$u$')
plt.ylabel('$s$')
plt.title('Spliced vs unspliced')
plt.show()