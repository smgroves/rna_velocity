import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
from velocity_tools.visual_utils import *
from matplotlib import cm
## Initialize variables for model
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'

# m = "one_gene" # model type
# tp = 6 # Number of timepoints
# seed = 621008493 # seed used to generate fake data
# N_cells = 10
# multifactor = 1
#
# m = "one_gene" # model type
# tp = 6 # Number of timepoints
# seed = 50119327 # seed used to generate fake data
# N_cells = 100
# multifactor = 1

# m = "two_gene" # model type
# tp = 6 # Number of timepoints
# seed = 93922959 # seed used to generate fake data
# N_cells = 20
# multifactor = 1
#
# m = "two_gene" # model type
# tp = 20 # Number of timepoints
# seed = 335410403 # seed used to generate fake data
# N_cells = 20
# multifactor = 1
#
m = "two_gene" # model type
tp = 4 # Number of timepoints
seed = 761358175 # seed used to generate fake data
N_cells = 500
multifactor = 1

################################################################################################################
''' 0. Import data '''
##Convert data from model to loom file and read loom as VelocytoLoom
data2loom(indir,model = m,timepoints=tp,seed = seed, N_cells = N_cells, multfactor=multifactor)
vlm = import_loom(op.join(indir,f"{m}/{seed}.loom"))

# If you want to do the calculation on only one timepoint, use filter_cells to choose that timepoint.
# vlm.filter_cells(vlm.ca['SampleID'] == '1')

################################################################################################################
'''1. Normalize data and perform PCA: input .U, .S --> .U_norm, .S_norm, .U_sz, .S_sz, .pca, .pcs
 (using .S_norm which is log_transformed, unlike .S_sz)'''
## Plot_fractions from velocyto_utils, which is an edit of the function originally from Velocyto to make it usable
# for fake data

# print("Plotting U/S fractions...")
# plt.figure()
# plot_fractions(vlm)
# plt.show()

vlm._normalize_S()
vlm._normalize_U()

PCA(vlm, n_components=2, pick_cutoff=False)

################################################################################################################
'''2. knn imputation and gamma fit'''
#If knn == False, the below function just sets Sx= S and Ux = U for future functions
# .Ux, .Sx, .Ux_sz, .Sx_sz set equal to .U, .S, .U_sz, .S_sz
vlm = knn_impute(vlm, n_comps=2, knn = False, normalized=False)
# Uses fake-imputed, normalized data .Sx_sz, .Us_sz to fit gamma, q, and R2

## TODO: if fit_offset = True and you use "unnormalized" data, the fit is very incorrect (the slope is usually
# negative). I haven't figured out why this would be yet, so I've set fit_offset = False. This uses a simply linear
# regression from scipy on Sx,Ux or Sx_sz, Ux_sz to fit gamma.
# Since weighted is set to False, this fit uses all of the datapoints to fit gamma. The default is weighted = True,
# weight= minmax_diag, which only uses the upper and lower percentiles to calculate gamma (assumed to be at steady
# state) [found in velocyto.estimation module]

vlm.fit_gammas(use_imputed_data=True, use_size_norm=False, weighted=False, fit_offset=False)
################################################################################################################
'''3. Velocity estimation'''
# Predicting U with which_S = "Sx_sz" sets which_U_for)pred = "Sx_sz"  --> .velocity = Ux_sz - Upred

vlm.predict_U(which_S="Sx")
vlm.calculate_velocity() # Only Sx_sz or Sx are implemented in velocyto

# print(colored("...saving hdf5 file as vlm", "blue"))
# vlm.to_hdf5(f"{seed}")

vlm.set_clusters(vlm.ca['SampleID'], colormap=cm.get_cmap('tab20b'))

################################################################################################################
'''4. Extrapolation at time t'''
# Creates delta_S = The variation in gene expression and Sx_sz_t = extrapolated data at time t
delta_t=.1

vlm.calculate_shift(assumption = 'constant_velocity') # Only Sx_sz or Sx are implemented in velocyto
vlm.extrapolate_cell_at_t(delta_t=delta_t) # Only Sx_sz or Sx are implemented in velocyto

################################################################################################################
'''5. Plotting data'''
# plot_2_gene_arrows(vlm, which_S="Sx")
vel = pd.DataFrame(vlm.velocity, index = vlm.ra['Gene'], columns = vlm.ca['CellID'])
vel.to_csv(f"{indir}/{m}/vel_{seed}_deltat_{delta_t}.csv")
# plot_phase_portraits(vlm, genes=['gene_0', 'gene_1'], which_S='Sx')

#attractor states are written as [[x1, x2], [y1,y2] for two attractors
plot_distr_over_time(vlm, attractors=[0.7494341170382606, 0.35890619161456505]) #for 2 genes

################################################################################################################
'''6. Use velocity to calculate transition probabilities and transition matrix'''
# estimate_transition_prob:
# Use cosine correlation to estimate transition probabilities for every cells to its embedding neighborhood

# estimate_transition_prob_smg(vlm, hidim = 'Sx', embed='pcs', transform = 'log', n_neighbors=100, sampled_fraction=1)
# vlm.calculate_embedding_shift(sigma_corr=0.05, expression_scaling=False)
# vlm.calculate_grid_arrows(smooth = 0.8, steps =(40,40), n_neighbors=100)
# vlm.prepare_markov()
# print(vlm.tr)