import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
from velocity_tools.visual_utils import *
from matplotlib import cm
## Initialize variables for model
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'
m = "two_gene" # model type
tp = 6 # Number of timepoints
seed = 93922959 # seed used to generate fake data
N_cells = 20

################################################################################################################
# 0. Import data

##Convert data from model to loom file and read loom as VelocytoLoom
data2loom(indir,model = m,timepoints=tp,seed = seed, N_cells = N_cells, multfactor=100)
vlm = import_loom(op.join(indir,f"{m}/{seed}.loom"))

################################################################################################################
# 1. Normalize data and perform PCA: input .U, .S --> .U_norm, .S_norm, .U_sz, .S_sz, .pca (using .S_norm which is
# log_transformed, unlike .S_sz)

## Plot_fractions from velocyto_utils, which is an edit of the function originally from Velocyto to make it usable
# for fake data

print("Plotting U/S fractions...")
plt.figure()
plot_fractions(vlm)
plt.show()

vlm.normalize("S", size=True, log=True)
vlm.normalize("U", size=True,  log=True)

PCA(vlm, n_components=2, pick_cutoff=False)

################################################################################################################
# 2. knn imputation and gamma fit
#If knn == False, the below function just sets Sx= S and Ux = U for future functions
# .U, .S, .Ux_sz, .Sx_sz set equal to .U, .S, .U_sz, .S_sz
vlm = knn_impute(vlm, n_comps=2, knn = False, normalized=False)

# Uses fake-imputed, normalized data .Sx_sz, .Us_sz to fit gamma, q, and R2
vlm.fit_gammas(use_imputed_data=True, use_size_norm=False, weighted=False)

################################################################################################################
# 4. Velocity estimation
# Predicting U with which_S = "Sx_sz" sets which_U_for)pred = "Sx_sz"  --> .velocity = Ux_sz - Upred

vlm.predict_U(which_S="Sx")
vlm.calculate_velocity() # Only Sx_sz or Sx are implemented in velocyto

# print(colored("...saving hdf5 file as vlm", "blue"))
# vlm.to_hdf5(f"{seed}")

vlm.set_clusters(vlm.ca['SampleID'], colormap=cm.get_cmap('Set2'))

################################################################################################################
#5. Extrapolation at time t
# Creates delta_S = The variation in gene expression and Sx_sz_t = extrapolated data at time t
delta_t=.1

vlm.calculate_shift(assumption = 'constant_velocity') # Only Sx_sz or Sx are implemented in velocyto
vlm.extrapolate_cell_at_t(delta_t=delta_t) # Only Sx_sz or Sx are implemented in velocyto

################################################################################################################
# 6. Plotting data
plot_2_gene_arrows(vlm, which_S="Sx")

vel = pd.DataFrame(vlm.velocity, index = vlm.ra['Gene'], columns = vlm.ca['CellID'])
vel.to_csv(f"{indir}/{m}/vel_{seed}_deltat_{delta_t}.csv")