# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:11:05 2018

@author: John
"""

# This file contains all model parameters. Change them here.



from choose_model_here import c, b, kA, kI


# Equation parameters.
G = 1
N = 4

b_x, b_y = b , b 
k_A_x, k_A_y = kA, kA 
k_I_x, k_I_y = kI, kI
#k_A_x = 0.2
#k_b_y = 0.42
#k_A_y = 0.23
#k = np.zeros((N+1,N+1))
#k[0,0] = k_b
#k[N,0] = k_b + k_A
#k[N,N] = k_A

#c = 0.5
#d_N0_x = 1
#d_0N_x = 1
#d_N0_y = 1
#d_0N_y = 1
#d = np.zeros((N+1, N+1))
#d[0,0] = 1
#d[N,0] = d_N0
#d[0,N] = d_0N
#d[N,N] = d_N0*d_0N

beta_x = 1                   # splicing rate
beta_y = 1

kp_x = 1                     # protein translation rate
kp_y = 1  

du_x = 1                     # unspliced degradation rate
ds_x = 1                     # spliced degradation rate
dp_x = 1                     # protein degradation rate 
du_y = 1                     
ds_y = 1                     
dp_y = 1                     

sig_x, sig_y = 1, 1

#G_x = 1                      # total number of available genes/gene sites
#G_y = 1
#
#N_x = 4                      # order of feedback (how many times protein can bind to gene) 
#N_y = 4
#k = np.zeros(N+1)       # transcription rate array
#k[0] = 10
#k[N] = 50
#
#beta_x = 1                   # splicing rate
#beta_y = 1
#
#kp_x = 1                     # protein translation rate
#kp_y = 1  
#
#du_x = 1                     # unspliced degradation rate
#ds_x = 1                     # spliced degradation rate
#dp_x = 1                     # protein degradation rate 
#du_y = 1                     
#ds_y = 1                     
#dp_y = 1                     
#
#c = np.zeros(N+1)              # Hill function coefficient array
#c[0] = 1
#c[N] = (1/25.32)**N    
#c[N] = (1/25.5)**N    
#c[N] = 1  


#b_r_val = 1              # Backward gene-protein binding rate
#b_r = np.zeros(N+1)        
#for i in range(1,N+1): 
#    b_r[i] = b_r_val       # For simplicity, assume all backward gene-protein binding rates are equal. 