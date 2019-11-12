# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:11:05 2018

@author: John
"""

# This file contains all model parameters. Change them here.

import numpy as np



# Equation parameters.
G = 1                      # total number of available genes/gene sites
N = 4                      # order of feedback (how many times protein can bind to gene) 
k = np.zeros(N+1)       # transcription rate array
k[0] = 10
k[N] = 50

beta = 1                   # splicing rate

kp = 1                     # protein translation rate

du = 1                     # unspliced degradation rate
ds = 1                     # spliced degradation rate
dp = 1                     # protein degradation rate 


c = np.zeros(N+1)              # Hill function coefficient array
c[0] = 1
c[N] = (1/25.32)**N    
#c[N] = (1/25.5)**N    
#c[N] = 1  


b_r_val = 1              # Backward gene-protein binding rate
b_r = np.zeros(N+1)        
for i in range(1,N+1): 
    b_r[i] = b_r_val       # For simplicity, assume all backward gene-protein binding rates are equal. 