import numpy as np
import loompy

from velocity_tools.velocity_utils import read_data
indir = '/Users/sarahmaddox/Documents/workspace/rna_velocity/output'

read_data(indir,model = 'one_gene',timepoint=0,seed = 32663012)
