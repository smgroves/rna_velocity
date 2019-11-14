import numpy as np
import loompy
import os.path as op
from velocity_tools.velocity_utils import *
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'
m = "one_gene"
tp = 0
seed =32663012
data2loom(indir,model = m,timepoint=tp,seed = seed)

vlm = import_loom(op.join(indir,f"{m}/time_{tp}/{seed}.loom"))
print(vlm)