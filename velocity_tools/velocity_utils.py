import numpy as np
import velocyto as vc
import loompy
import os.path as op
import pandas as pd


def read_data(indir, model, timepoint, seed = 0, N_cells = 20, delim = ','):
    #Indir is directory path to output folder
    # model is one_gene, two_gene, etc.
    #timepoint is a single number, starting at 0 for first timepoint
    read_dir = f"{indir}/{model}/time_{str(timepoint)}"
    print("Data should have rows = genes, columns = cells, for a single timepoint.")
    data = pd.read_csv(op.join(read_dir, f's_{N_cells}_cells_{seed}.csv'), index_col=0, header=0, delimiter=delim)
    filename = f"{seed}.loom"
    matrix = np.array(data)
    row_attrs = {"gene": data.index.values}
    col_attrs = {"cellID": data.columns.values}
    loompy.create(filename, matrix, row_attrs, col_attrs)
