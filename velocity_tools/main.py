import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'

#Convert data from model to loom file and read loom as VelocytoLoom
m = "two_gene"
tp = 2
seed =50307154
# seed = 32663012
# data2loom(indir,model = m,timepoints=tp,seed = seed)
vlm = import_loom(op.join(indir,f"{m}/{seed}.loom"))

print("Plotting U/S fractions...")
plt.figure()
plot_fractions(vlm)
plt.show()