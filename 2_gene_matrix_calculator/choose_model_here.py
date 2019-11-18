# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 22:18:22 2019

@author: johnv
"""

model = 'two'
# options: 'two', 'three', 'four'

c = 0.5

if model=='two':
    b = (6/4)*c
    kA = (1/4)*c
    kI = (1)*c
    p_x0, p_y0 = 0.7494341170382606, 0.35890619161456505
elif model == 'three':
    b = (1/2)*c
    kA = (2.1/2)*c
    kI = (0.8/2)*c 
    p_x0, p_y0 = 0.6145058082866004, 0.11227898726758452
elif model == 'four':
    b = (1/4)*c
    kA = (7/4)*c
    kI = (1/10)*c 
    p_x0, p_y0 = 0.933402929122054, 0.07935865405488161