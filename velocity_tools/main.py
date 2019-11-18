import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
from velocity_tools.visual_utils import *

## Initialize variables for model
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'
m = "two_gene" # model type
tp = 2 # Number of timepoints
seed = 50307154 # seed used to generate fake data
N_cells = 20

##Convert data from model to loom file and read loom as VelocytoLoom
data2loom(indir,model = m,timepoints=tp,seed = seed, N_cells = N_cells, multfactor=100)
vlm = import_loom(op.join(indir,f"{m}/{seed}.loom"))

## Plot_fractions from velocyto_utils, which is an edit of the function originally from Velocyto to make it usable
# for fake data

print("Plotting U/S fractions...")
plt.figure()
plot_fractions(vlm)
plt.show()

################################################################################################################
# 1. Normalize data and perform PCA: input .U, .S --> .U_norm, .S_norm, .U_sz, .S_sz, .pca (using .S_norm which is
# log_transformed, unlike .S_sz)

vlm.normalize()
PCA(vlm, n_components=2, pick_cutoff=False)

################################################################################################################
# 2. knn imputation and gamma fit
#If knn == False, the below function just sets Sx= S and Ux = U for future functions
# .U, .S, .Ux_sz, .Sx_sz set equal to .U, .S, .U_sz, .S_sz
vlm = knn_impute(vlm, n_comps=2, knn = False)

# Uses fake-imputed, normalized data .Sx_sz, .Us_sz to fit gamma, q, and R2
vlm.fit_gammas(use_imputed_data=True, use_size_norm=True, weighted=False)

################################################################################################################
# 4. Velocity estimation

vlm.set_clusters(vlm.ca['SampleID'])

#The below code predicts U, calculates velocity, calculates shift, and extrapolates cell at time t
# Calculate_velocity will only use Sx or Sx_sz (not S)
vlm = calc_velocity_and_shift(vlm, delta_t = 1, which_s="Sx_sz")
plot_original_arrows(vlm)

################################################################################################################
#5. Extrapolation at time t
