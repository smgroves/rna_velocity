import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
from velocity_tools.visual_utils import *
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'

#Convert data from model to loom file and read loom as VelocytoLoom
m = "two_gene" # model type
tp = 2 # Number of timepoints
seed = 50307154 # seed used to generate fake data

data2loom(indir,model = m,timepoints=tp,seed = seed, N_cells = 20, multfactor=100)
vlm = import_loom(op.join(indir,f"{m}/{seed}.loom"))

print("Plotting U/S fractions...")
plt.figure()
plot_fractions(vlm)
plt.show()

vlm = norm_vlm(vlm)
vlm = gamma_fit(vlm, n_components=2, pick_cutoff=False, knn = False, balanced=False)

vlm.Sx = vlm.S
vlm.Ux = vlm.U
vlm.set_clusters(vlm.ca['SampleID'])

vlm = calc_velocity_and_shift(vlm, delta_t = 1)
print(vlm.Sx)
print(vlm.Sx_t)

# vlm.estimate_transition_prob(hidim = 'Sx', embed = 'pca', transform = None)
# vlm.calculate_embedding_shift()
# vlm.calculate_grid_arrows(embed='pca')
plot_original_arrows(vlm)
#plot portrait figure: U vs. S for a specific gene with gamma fit overlaid
plt.figure()
plt.scatter(x = vlm.velocity[0], y = vlm.velocity[1])
plt.show()
