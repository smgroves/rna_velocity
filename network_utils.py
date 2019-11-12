import os.path as op
import os
import numpy as np
import pandas as pd

def save_data(u_list,s_list,dir_prefix, N_cells, N_timepoints, model):
    for tp in range(N_timepoints):
        try:
            os.mkdir(op.join(dir_prefix, f'{model}/time_{tp}/'))
        except FileExistsError:
            pass


        data = pd.DataFrame([u[:,tp] for u in u_list], columns = [f"cell_{str(i)}" for i in range(N_cells)])
        data.index = [f"gene_{str(i)}" for i in range(len(u_list))]
        print(data)

        data.to_csv(op.join(dir_prefix, f'{model}/time_{tp}/u_{N_cells}_cells.csv'))

        data = pd.DataFrame([s[:,tp] for s in s_list], columns = [f"cell_{str(i)}" for i in range(N_cells)])
        data.index = [f"gene_{str(i)}" for i in range(len(u_list))]

        data.to_csv(op.join(dir_prefix,f'{model}/time_{tp}/s_{N_cells}_cells.csv'))
