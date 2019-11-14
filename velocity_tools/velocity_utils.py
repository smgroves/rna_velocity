import numpy as np
import velocyto as vcy
import loompy
import os.path as op
import pandas as pd


def data2loom(indir, model, timepoints = 2, seed = 0, N_cells = 20, delim = ','):
    #Indir is directory path to output folder
    # model is one_gene, two_gene, etc.
    #timepoint is a single number, starting at 0 for first timepoint
    datas = pd.DataFrame()
    datau = pd.DataFrame()
    time = []
    for timepoint in range(timepoints):
        read_dir = f"{indir}/{model}/time_{str(timepoint)}"
        print("Data should have rows = genes, columns = cells, for a single timepoint.")
        tmp_datas = pd.read_csv(op.join(read_dir, f's_{N_cells}_cells_{seed}.csv'), index_col=0, header=0,
                              delimiter=delim)
        tmp_datau = pd.read_csv(op.join(read_dir, f'u_{N_cells}_cells_{seed}.csv'), index_col=0, header=0,
                              delimiter=delim)
        if timepoint == 0:
            datas = tmp_datas
            datau = tmp_datau
        else:
            datas = datas.merge(tmp_datas, on = datas.index)
            datau = datau.merge(tmp_datau, on = datau.index)
        for i in [str(timepoint)]*N_cells:
            time.append(i)
    datas.index = datas['key_0']
    datas = datas.drop('key_0',axis = 1)
    datau.index = datau['key_0']
    datau = datau.drop('key_0', axis = 1)
    total = np.array(datas)+np.array(datau)
    filename = f"{seed}.loom"
    ra = {"Gene": datas.index.values}
    ca = {"CellID": datas.columns.values}
    out_dir = f"{indir}/{model}/"

    try:
        ds = loompy.create(op.join(read_dir,filename), total, ra, ca)

        for layer_name, matrix in zip(['S','U'], [datas,datau]):
            ds.set_layer(name=layer_name, matrix=matrix, dtype='float32')
        ds.attrs["velocyto.__version__"] = vcy.__version__
        ds.attrs["velocyto.logic"] = 'Default'
        ds.close()

    except AttributeError:
        # If user is using loompy2
        # NOTE maybe this is not super efficient if the type and order are already correct
        layers = {'spliced':np.array(datas), 'unspliced': np.array(datau), 'ambiguous': np.zeros(np.array(datau).shape)}
        tmp_layers = {"": total.astype("float32", order="C", copy=False)}

        tmp_layers.update({layer_name: layers[layer_name].astype('uint32', order="C", copy=False) for
                           layer_name in layers.keys()})
        loompy.create(filename=op.join(out_dir,filename), layers=tmp_layers, row_attrs=ra, col_attrs=ca,
                      file_attrs={"velocyto.__version__": vcy.__version__, "velocyto.logic": 'Default'})

def import_loom(loom):
    print('Importing loom file...')
    vlm = vcy.VelocytoLoom(loom)
    print("Column attributes:", vlm.ca)
    print("Row attributes:", vlm.ra['Gene'])
    print("Spliced shape:", vlm.S.shape, "Unspliced shape:", vlm.U.shape)
    return vlm