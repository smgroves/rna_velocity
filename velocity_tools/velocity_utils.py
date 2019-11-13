import numpy as np
import velocyto as vc
import loompy
import os.path as op
import pandas as pd

filename = "test.loom"
matrix = np.arange(10000).reshape(100,100)
row_attrs = { "SomeRowAttr": np.arange(100) }
col_attrs = { "SomeColAttr": np.arange(100) }
loompy.create(filename, matrix, row_attrs, col_attrs)
ds = loompy.connect("test.loom")


def read_data(indir, model, timepoint, N_cells = 20, delim = ','):
    #Indir is directory path to output folder
    # model is one_gene, two_gene, etc.
    #timepoint is a single number, starting at 0 for first timepoint
    read_dir = f"{indir}/{model}/time_{timepoint}"
    print("Data should have rows = genes, columns = cells, for a single timepoint.")
    data = pd.read_csv(op.join(read_dir, f's_{N_cells}_cells.csv'), index_col=0, header=0, delimiter=delim)