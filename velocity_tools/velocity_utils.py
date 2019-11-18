import numpy as np
import velocyto as vcy
import loompy
import os.path as op
import pandas as pd
import matplotlib.pyplot as plt
from termcolor import colored
import os
import igraph
from sklearn.neighbors import NearestNeighbors
import rpy2.robjects as robj
from rpy2.robjects.packages import importr

def data2loom(indir, model, timepoints = 2, seed = 0, N_cells = 20, delim = ',', multfactor = 1):
    '''
    :param indir: directory where "output" folder is
    :param model: one_gene, two_gene
    :param timepoints:  number of timepoints
            Default = 2
    :param seed: random seed used to generate fake data (in the name of the fake data files)
            Default = 0
    :param N_cells: number of fake cells for which data was generated
            Default = 20
    :param delim: delimiter in fake data file
            Default = ','
    :param multfactor: multiplicative factor for CME models, which output concentration, to generate "counts" data
            Default = 100
    :return: Nothing, but generates a file {seed}.loom in the model folder with the characteristics
            layers: X (spliced + unspliced counts), spliced, unspliced, ambiguous (currently 0s matrix)
            ra: Genes
            ca: CellID and SampleID (timepoint)
    '''

    # model is one_gene, two_gene, etc.
    #timepoint is a single number, starting at 0 for first timepoint
    datas = pd.DataFrame()
    datau = pd.DataFrame()
    time = []
    for timepoint in range(timepoints):
        read_dir = f"{indir}/{model}/time_{str(timepoint)}"
        print("Data should have rows = genes and columns = cells.")
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
            datas.index = datas['key_0']
            datas = datas.drop('key_0', axis=1)
            datau.index = datau['key_0']
            datau = datau.drop('key_0', axis=1)
        for i in [str(timepoint)]*N_cells:
            time.append(i)

    datas = datas*multfactor
    datau = datau*multfactor
    total = np.array(datas)+np.array(datau)
    filename = f"{seed}.loom"
    ra = {"Gene": datas.index.values}
    ca = {"CellID": datas.columns.values, 'SampleID': np.array(time)}
    out_dir = f"{indir}/{model}/"

    try:
        ds = loompy.create(op.join(read_dir,filename), total, ra, ca)

        for layer_name, matrix in zip(['spliced','unspliced'], [datas,datau]):
            ds.set_layer(name=layer_name, matrix=matrix, dtype='float32')
        ds.attrs["velocyto.__version__"] = vcy.__version__
        ds.attrs["velocyto.logic"] = 'Default'
        ds.attrs['SampleID'] = time
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
    '''
    :param loom: name of loom file to import with full path specified
    :return: vlm, a VelocytoLoom object
    '''
    print('Importing loom file...')
    vlm = vcy.VelocytoLoom(loom)
    print("Column attributes:", vlm.ca.keys())
    print("Row attributes:", vlm.ra['Gene'])
    print("Spliced shape:", vlm.S.shape, "Unspliced shape:", vlm.U.shape)
    return vlm

def plot_fractions(vlm, save2file: str = None):
    """Plots a barplot showing the abundance of spliced/unspliced molecules in the dataset
    Arguments
    ---------
    save2file: str (default: None)
        If not None specifies the file path to which plots get saved
    Returns
    -------
    Nothing, it plots a barplot
    """
    plt.figure(figsize=(4, 4))
    try:
        chips, chip_ix = np.unique(self.ca["SampleID"], return_inverse=1)
    except KeyError:
        chips, chip_ix = np.unique([i.split(":")[0] for i in self.ca["CellID"]], return_inverse=1)
    n = len(chips)
    for i in np.unique(chip_ix):
        tot_mol_cell_submatrixes = [X[:, chip_ix == i].sum(0) for X in [self.S, self.A, self.U]]
        total = np.sum(tot_mol_cell_submatrixes, 0)
        _mean = [np.mean(j / total) for j in tot_mol_cell_submatrixes]
        _std = [np.std(j / total) for j in tot_mol_cell_submatrixes]

        plt.ylabel("Fraction")
        plt.bar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, 0.5 / (n * 1.05), label=chips[i])
        plt.errorbar(np.linspace(-0.2, 0.2, n)[i] + np.arange(3), _mean, _std, c="k", fmt="none", lw=1, capsize=2)

        # Hide the right and top spines
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('bottom')
        plt.gca().spines['left'].set_bounds(0, 0.8)
        plt.legend()

    plt.xticks(np.arange(3), ["spliced", "ambiguous", "unspliced"])
    plt.tight_layout()
    if save2file:
        plt.savefig(save2file, bbox_inches="tight")

def PCA(vlm, n_components, pick_cutoff = True):
    print('...Performing PCA')
    vlm.perform_PCA(n_components=n_components)
    if pick_cutoff == True:
        plt.figure()
        plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
        try:
            n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
        except IndexError: n_comps = n_components
        plt.axvline(n_comps, c="k")
        plt.show()
        # variable: change the 0.002 below to reflect where to cut off the pca; 0.002 is the minimum difference
        # you want when adding another PC dimension; lower = more dimensions; higher = fewer dimensions = less noise
        co = input("Input cutoff for change in cumulative explained variance; default is 0.002")
        n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>float(co)))[0][0]
        plt.figure()
        plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
        plt.axvline(n_comps, c="k")
        plt.show()
    return vlm

def knn_impute(vlm, n_comps = 100, knn = True, b_sight = 100, k = 50, balanced = True, normalized = True):
    print("Fitting gamma (degradation rate) parameter (~ u/s at steady state)...")
    if knn == True:
        print('...performing kNN imputation (runs slowly)')
        if balanced == True:
            vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=b_sight, b_maxl=1500, n_jobs=16)
        else:
            vlm.knn_imputation(n_pca_dims=n_comps, k = k)
    else:
        vlm.Sx = vlm.S
        vlm.Ux = vlm.U
        if normalized == True:
            vlm.Sx_sz = vlm.S_sz
            vlm.Ux_sz = vlm.U_sz
    return vlm

def fit_transform_clust(vlm, fit_type = "UMAP", clust = False, n_clust = 4):
    seed = 2
    print("Fitting and transforming the data...")
    if fit_type == "UMAP":
        import umap
        um = umap.UMAP()
        vlm.ts = um.fit_transform(vlm.pcs[:, :20])
        print(f"...embedding by {fit_type} saved as vlm.ts")
    elif fit_type == "TSNE":
        from sklearn.manifold import TSNE
        bh_tsne = TSNE()
        vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :20])
        print(f"...embedding by {fit_type} saved as vlm.ts")
    if clust == True:
        print("...clustering based on embedding saved under vlm.ca['ClusterName]")
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=n_clust, random_state=0).fit(vlm.ts)
        vlm.ca['ClusterName'] = km.labels_
        vlm.set_clusters(vlm.ca["ClusterName"])

    return vlm

def calc_velocity_and_shift(vlm, name = None, assumption = "constant_velocity", which_s = 'Sx', delta_t = 1):

    print("Calculating velocity...")
    vlm.predict_U(which_S=which_s)
    vlm.calculate_velocity()
    vlm.calculate_shift(assumption=assumption)
    print(f"Extrapolating cell at t with assumption {assumption}...")
    vlm.extrapolate_cell_at_t(delta_t=delta_t)
    if name == None:
        print(colored("...saving hdf5 file as vlm","blue"))
        vlm.to_hdf5("vlm")
    else:
        print(colored(f"...saving hdf5 file as {name}", "blue"))
        vlm.to_hdf5(name)
    return vlm

def estimate_trans_prob(vlm, embed = "ts"):
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
    print(f"Estimating transition probabilities with embed = {embed}...")
    # Try running outside kernel (run in ipython/terminal)
    # import os
    # os.environ['KMP_DUPLICATE_LIB_OK']='True'
    vlm.estimate_transition_prob(hidim="Sx_sz", embed=embed, transform="sqrt", psc=1, n_neighbors=2000, knn_random=True,
                                 sampled_fraction=0.5)
    # vlm.to_hdf5('vlm_prob_TK02')
    return vlm

def calc_embedding_shift(vlm):

    print("Calculating shift...")
    vlm.calculate_embedding_shift(sigma_corr=0.05, expression_scaling=True)
    print("Calculating grid arrows...")
    vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
    plt.figure()
    plt.title("Flow_norm_magnitude")
    plt.hist(vlm.flow_norm_magnitude)
    plt.show()
    return vlm


def save_two_timepoints(vlm, name = "velocyto_data"):

    print("Saving t0 as {}...")
    t = "t0"
    pd.DataFrame(vlm.Sx_sz).to_csv(f'{name}_{t}.csv')
    print("Saving t1 as {}...")
    t = "t1"
    pd.DataFrame(vlm.Sx_sz_t).to_csv(f'{name}_{t}.csv')


def array_to_rmatrix(X):
    nr, nc = X.shape
    xvec = robj.FloatVector(X.transpose().reshape((X.size)))
    xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
    return xr


def principal_curve(X, pca=True):
    """
    input : numpy.array
    returns:
    Result::Object
        Methods:
        projections - the matrix of the projectiond
        ixsort - the order ot the points (as in argsort)
        arclength - the lenght of the arc from the beginning to the point
    """
    # convert array to R matrix
    xr = array_to_rmatrix(X)

    if pca:
        # perform pca
        t = robj.r.prcomp(xr)
        # determine dimensionality reduction
        usedcomp = max(sum(np.array(t[t.names.index('sdev')]) > 1.1), 4)
        usedcomp = min([usedcomp, sum(np.array(t[t.names.index('sdev')]) > 0.25), X.shape[0]])
        Xpc = np.array(t[t.names.index('x')])[:, :usedcomp]
        # convert array to R matrix
        xr = array_to_rmatrix(Xpc)

    # import the correct namespace
    princurve = importr("princurve", on_conflict="warn")

    # call the function
    fit1 = princurve.principal_curve(xr)

    # extract the outputs
    class Results:
        pass

    results = Results()
    results.projections = np.array(fit1[0])
    results.ixsort = np.array(fit1[1]) - 1  # R is 1 indexed
    results.arclength = np.array(fit1[2])
    results.dist = np.array(fit1[3])

    if pca:
        results.PCs = np.array(xr)  # only the used components

    return results

def pc_clusters(vlm):
    nn = NearestNeighbors(50)
    nn.fit(vlm.pcs[:,:4])
    knn_pca = nn.kneighbors_graph(mode='distance')
    knn_pca = knn_pca.tocoo()
    G = igraph.Graph(list(zip(knn_pca.row, knn_pca.col)), directed=False, edge_attrs={'weight': knn_pca.data})
    VxCl = G.community_multilevel(return_levels=False, weights="weight")
    labels = np.array(VxCl.membership)

    from numpy_groupies import aggregate_np


    # Make a principal curve for a subset of the data based on cluster
    # pc_obj = principal_curve(vlm.pcs[vlm.cluster_labels == b"8",:4], False)

    pc_obj = principal_curve(vlm.pcs[:,:8], False)
    pc_obj.arclength = np.max(pc_obj.arclength) - pc_obj.arclength
    labels = np.argsort(np.argsort(aggregate_np(labels, pc_obj.arclength, func=np.median)))[labels]
    manual_annotation = {str(i):[i] for i in labels}
    annotation_dict = {v:k for k, values in manual_annotation.items() for v in values }
    clusters = np.array([annotation_dict[i] for i in labels])
    colors20 = np.vstack((plt.cm.tab20b(np.linspace(0., 1, 20))[::2], plt.cm.tab20c(np.linspace(0, 1, 20))[1::2]))
    vlm.set_clusters(clusters, cluster_colors_dict={k:colors20[v[0] % 20,:] for k,v in manual_annotation.items()})
    vlm.ca['Clusters'] = [int(i) for i in clusters]

# must be run after calculating velocity and shift!!
def average_vel_per_cluster(vlm, embed = True):
    if embed == True:
        magnitude_embed = []
        for i in range(len(vlm.delta_embedding[:, 1])):
            magnitude_embed.append(np.linalg.norm(vlm.delta_embedding[i, 1:10]))
    else:
        magnitude_vec = []
        for i in range(len(vlm.delta_S[0])):
            magnitude_vec.append(np.linalg.norm(vlm.delta_S[:, i]))
    # calculate average velocity for each cluster
    ave_vel = []
    for q in list(len(vlm.embedding[:, 0] + 1)):
        ind = [i for i, x in enumerate(vlm.cluster_labels == str(q)) if x]
        ave_vel.append(np.average([magnitude_embed[x] for x in ind]))

    col_ave_vel = []
    for i in range(len(vlm.embedding[:, 0] + 1)):
        cl = int(vlm.cluster_labels[i])
        col_ave_vel.append(ave_vel[cl])
    #plot
    plt.figure(None, (9, 9))
    ax = plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
                     c=col_ave_vel, cmap='RdBu', alpha=0.5, s=20, lw=0,
                     edgecolor='', rasterized=True)

    for i in range(max(vlm.ca["Clusters"]) + 1):
        pc_m = np.median(vlm.pcs[[x == str(i) for x in vlm.cluster_labels], :], 0)
        plt.text(pc_m[0], pc_m[1], str(vlm.cluster_labels[[x == str(i) for x in vlm.cluster_labels]][0]),
                 fontsize=13, bbox={"facecolor": "w", "alpha": 0.6})
    plt.axis("off")
    cbar = plt.colorbar(ax)
    cbar.ax.set_yticks(ave_vel, fontsize=8)
    plt.show()
    return col_ave_vel