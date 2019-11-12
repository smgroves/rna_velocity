import os.path as op
import os
import numpy as np

def save_data(u,s,dir_prefix, N_cells, N_timepoints, deltat, model):
    try:
        os.mkdir(op.join(dir_prefix, f'{model}'))
    except FileExistsError:
        pass
    outfile = open(op.join(dir_prefix, f"{model}/u_{N_cells}cells_{N_timepoints}timepoints_{deltat}timestep.csv"),
                   "w+")
    data = np.array(u)
    outfile.write(str(data))
    outfile.close()
    outfile = open(op.join(dir_prefix, f"{model}/s_{N_cells}cells_{N_timepoints}timepoints_{deltat}timestep.csv"),
                   "w+")
    data = np.array(s)
    outfile.write(str(data))
    outfile.close()