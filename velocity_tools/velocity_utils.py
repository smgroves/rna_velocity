import numpy as np
import velocyto as vcy
import loompy
import os.path as op
import pandas as pd
import matplotlib.pyplot as plt
from termcolor import colored
import os
import igraph
import scipy as sc
from sklearn.neighbors import NearestNeighbors
import rpy2.robjects as robj
from rpy2.robjects.packages import importr
import random
import logging
from numba import jit
from velocyto.estimation import colDeltaCor, colDeltaCorSqrt, colDeltaCorLog10, colDeltaCorpartial, colDeltaCorSqrtpartial, \
    colDeltaCorLog10partial


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
        # print("Data should have rows = genes and columns = cells.")
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

        tmp_layers.update({layer_name: layers[layer_name].astype('float32', order="C", copy=False) for
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

def plot_fractions(self, save2file: str = None):
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

def PCA(vlm, n_components, pick_cutoff = True, which_S = 'S_norm'):
    '''

    :param vlm: VelocytoLoom object
    :param n_components: number of components for the PCA (initially)
    :param pick_cutoff: if True, user will be prompted in the console to provide a cutoff of the difference
    in cumulative explained variance necessary to keep PCs (i.e., this specifies how many PCs to keep in future steps)
    :return: vlm object with .pcs attribute, and plots cumulative explained variance if pick_cutoff == True
    '''
    print('...Performing PCA')
    vlm.perform_PCA(n_components=n_components,which = which_S)
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
    '''

    :param vlm: Velocyto object
    :param n_comps: number of components to use for the imputation
    :param knn: if True, perform imputation. Otherwise, add attributes for future compatibility.
    :param b_sight: If balanced == True, sight of each cell.
    :param k: number of nearest neighbors to calculate
    :param balanced: if True, run bKNN, otherwise run kNN
    :param normalized: True if a normalization step has previously been completed to generate .S_sz and U_sz. Only
    necessary if kNN == False.
    :return:
    '''
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


def estimate_transition_prob_smg(self, hidim= "Sx_sz", embed= "ts", transform = "sqrt",
                             ndims = None, n_sight = None, psc = None,
                             knn_random= True, sampled_fraction = 0.3,
                             sampling_probs = (0.5, 0.1), max_dist_embed =  None,
                             n_jobs= 4, threads = None, calculate_randomized = True,
                             random_seed = 15071990, **kwargs):
    """Use correlation to estimate transition probabilities for every cells to its embedding neighborhood

    Arguments
    ---------
    hidim: str, default="Sx_sz"
        The name of the attribute containing the high dimensional space. It will be retrieved as getattr(self, hidim)
        The updated vector at time t is assumed to be getattr(self, hidim + "_t")
        Appending .T to the string will transpose the matrix (useful in case we want to use S or Sx)
    embed: str, default="ts"
        The name of the attribute containing the embedding. It will be retrieved as getattr(self, embed)
    transform: str, default="sqrt"
        The transformation that is applies on the high dimensional space.
        If None the raw data will be used
    ndims: int, default=None
        The number of dimensions of the high dimensional space to work with. If None all will be considered
        It makes sense only when using principal components
    n_sight: int, default=None (also n_neighbors)
        The number of neighbors to take into account when performing the projection
    psc: float, default=None
        pseudocount added in variance normalizing transform
        If None, 1 would be used for log, 0 otherwise
    knn_random: bool, default=True
        whether to random sample the neighborhoods to speedup calculation
    sampling_probs: Tuple, default=(0.5, 1)
    max_dist_embed: float, default=None
        CURRENTLY NOT USED
        The maximum distance allowed
        If None it will be set to 0.25 * average_distance_two_points_taken_at_random
    n_jobs: int, default=4
        number of jobs to calculate knn
        this only applies to the knn search, for the more time consuming correlation computation see threads
    threads: int, default=None
        The threads will be used for the actual correlation computation by default half of the total.
    calculate_randomized: bool, default=True
        Calculate the transition probabilities with randomized residuals.
        This can be plotted downstream as a negative control and can be used to adjust the visualization scale of the velocity field.
    random_seed: int, default=15071990
        Random seed to make knn_random mode reproducible

    Returns
    -------
    """
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'  # needed so that python doesn't crash
    print(f"Estimating transition probabilities with embed = {embed}...")
    # Try running outside kernel (run in ipython/terminal)
    # import os
    # os.environ['KMP_DUPLICATE_LIB_OK']='True'
    random.seed(random_seed)
    self.which_hidim = hidim

    if "n_neighbors" in kwargs:
        n_neighbors = kwargs.pop("n_neighbors")
        if len(kwargs) > 0:
            logging.warning(f"keyword arguments were passed but could not be interpreted {kwargs}")
    else:
        n_neighbors = None

    if n_sight is None and n_neighbors is None:
        n_neighbors = int(self.S.shape[1] / 5)

    if (n_sight is not None) and (n_neighbors is not None) and n_neighbors != n_sight:
        raise ValueError(
            "n_sight and n_neighbors are different names for the same parameter, they cannot be set differently")

    if n_sight is not None and n_neighbors is None:
        n_neighbors = n_sight

    if psc is None:
        if transform == "log" or transform == "logratio":
            psc = 1.
        elif transform == "sqrt":
            psc = 1e-10  # for numerical stablity
        else:  # transform == "linear":
            psc = 0

    if knn_random:
        np.random.seed(random_seed)
        self.corr_calc = "knn_random"
        if "pcs" in hidim:  # sic
            hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
            hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
        else:
            if ndims is not None:
                raise ValueError(
                    f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
            hi_dim = getattr(self, hidim)  # [:, :ndims]
            hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
            if calculate_randomized:
                self.delta_S_rndm = np.copy(self.delta_S)
                permute_rows_nsign(self.delta_S_rndm)
                hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm

        embedding = getattr(self, embed)
        self.embedding = embedding
        logging.debug("Calculate KNN in the embedding space")
        nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
        nn.fit(embedding)  # NOTE should support knn in high dimensions
        self.embedding_knn = nn.kneighbors_graph(mode="connectivity")

        # Pick random neighbours and prune the rest
        neigh_ixs = self.embedding_knn.indices.reshape((-1, n_neighbors + 1))
        p = np.linspace(sampling_probs[0], sampling_probs[1], neigh_ixs.shape[1])
        p = p / p.sum()

        # There was a problem of API consistency because the random.choice can pick the diagonal value (or not)
        # resulting self.corrcoeff with different number of nonzero entry per row.
        # Not updated yet not to break previous analyses
        # Fix is substituting below `neigh_ixs.shape[1]` with `np.arange(1,neigh_ixs.shape[1]-1)`
        # I change it here since I am doing some breaking changes
        sampling_ixs = np.stack((np.random.choice(neigh_ixs.shape[1],
                                                  size=(int(sampled_fraction * (n_neighbors + 1)),),
                                                  replace=False,
                                                  p=p) for i in range(neigh_ixs.shape[0])), 0)
        self.sampling_ixs = sampling_ixs
        neigh_ixs = neigh_ixs[np.arange(neigh_ixs.shape[0])[:, None], sampling_ixs]
        nonzero = neigh_ixs.shape[0] * neigh_ixs.shape[1]
        self.embedding_knn = sc.sparse.csr_matrix((np.ones(nonzero),
                                                neigh_ixs.ravel(),
                                                np.arange(0, nonzero + 1, neigh_ixs.shape[1])),
                                               shape=(neigh_ixs.shape[0],
                                                      neigh_ixs.shape[0]))

        logging.debug(f"Correlation Calculation '{self.corr_calc}'")
        if transform == "log":
            delta_hi_dim = hi_dim_t - hi_dim
            self.corrcoef = colDeltaCorLog10partial(hi_dim,
                                                    np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim),
                                                    neigh_ixs, threads=threads, psc=psc)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                self.corrcoef_random = colDeltaCorLog10partial(hi_dim,
                                                               np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(
                                                                   delta_hi_dim_rndm), neigh_ixs, threads=threads,
                                                               psc=psc)
        elif transform == "logratio":
            log2hidim = np.log2(hi_dim + psc)
            delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
            self.corrcoef = colDeltaCorpartial(log2hidim, delta_hi_dim, neigh_ixs, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + psc) - log2hidim
                self.corrcoef_random = colDeltaCorpartial(log2hidim, delta_hi_dim_rndm, neigh_ixs, threads=threads)
        elif transform == "linear":
            self.corrcoef = colDeltaCorpartial(hi_dim, hi_dim_t - hi_dim, neigh_ixs, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                self.corrcoef_random = colDeltaCorpartial(hi_dim, hi_dim_t_rndm - hi_dim, neigh_ixs,
                                                          threads=threads)
        elif transform == "sqrt":
            delta_hi_dim = hi_dim_t - hi_dim
            self.corrcoef = colDeltaCorSqrtpartial(hi_dim.astype(np.float64).copy(order='C'),
                                                   (np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(
                                                       delta_hi_dim)).copy(order='C').astype(np.float64),
                                                   neigh_ixs, threads=threads, psc=psc)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                self.corrcoef_random = colDeltaCorSqrtpartial(hi_dim,
                                                              np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(
                                                                  delta_hi_dim_rndm), neigh_ixs, threads=threads,
                                                              psc=psc)
        else:
            raise NotImplementedError(f"transform={transform} is not a valid parameter")
        np.fill_diagonal(self.corrcoef, 0)
        if np.any(np.isnan(self.corrcoef)):
            self.corrcoef[np.isnan(self.corrcoef)] = 1
            logging.warning(
                "Nans encountered in corrcoef and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
        if calculate_randomized:
            np.fill_diagonal(self.corrcoef_random, 0)
            if np.any(np.isnan(self.corrcoef_random)):
                self.corrcoef_random[np.isnan(self.corrcoef_random)] = 1
                logging.warning(
                    "Nans encountered in corrcoef_random and corrected to 1s. If not identical cells were present it is probably a small isolated cluster converging after imputation.")
        logging.debug(f"Done Correlation Calculation")
    else:
        self.corr_calc = "full"
        if "pcs" in hidim:  # sic
            hi_dim = np.array(getattr(self, hidim).T[:, :ndims], order="C")
            hi_dim_t = np.array(getattr(self, hidim + "_t").T[:, :ndims], order="C")
        else:
            if ndims is not None:
                raise ValueError(
                    f"ndims was set to {ndims} but hidim != 'pcs'. Set ndims = None for hidim='{hidim}'")
            hi_dim = getattr(self, hidim)  # [:, :ndims]
            hi_dim_t = hi_dim + self.used_delta_t * self.delta_S  # [:, :ndims] [:, :ndims]
            if calculate_randomized:
                self.delta_S_rndm = np.copy(self.delta_S)
                permute_rows_nsign(self.delta_S_rndm)
                hi_dim_t_rndm = hi_dim + self.used_delta_t * self.delta_S_rndm

        embedding = getattr(self, embed)
        self.embedding = embedding
        logging.debug("Calculate KNN in the embedding space")
        nn = NearestNeighbors(n_neighbors=n_neighbors + 1, n_jobs=n_jobs)
        nn.fit(embedding)
        self.embedding_knn = nn.kneighbors_graph(mode="connectivity")

        logging.debug("Correlation Calculation 'full'")
        if transform == "log":
            delta_hi_dim = hi_dim_t - hi_dim
            self.corrcoef = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim),
                                             threads=threads, psc=psc)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                self.corrcoef_random = colDeltaCorLog10(hi_dim, np.log10(np.abs(delta_hi_dim_rndm) + psc) * np.sign(
                    delta_hi_dim_rndm), threads=threads, psc=psc)
        elif transform == "logratio":
            log2hidim = np.log2(hi_dim + psc)
            delta_hi_dim = np.log2(np.abs(hi_dim_t) + psc) - log2hidim
            self.corrcoef = colDeltaCor(log2hidim, delta_hi_dim, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = np.log2(np.abs(hi_dim_t_rndm) + 1) - log2hidim
                self.corrcoef_random = colDeltaCor(log2hidim, delta_hi_dim_rndm, threads=threads)
        elif transform == "linear":
            self.corrcoef = colDeltaCor(hi_dim, hi_dim_t - hi_dim, threads=threads)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                self.corrcoef_random = colDeltaCor(hi_dim, hi_dim_t_rndm - hi_dim, threads=threads, psc=psc)
        elif transform == "sqrt":
            delta_hi_dim = hi_dim_t - hi_dim
            self.corrcoef = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim) + psc) * np.sign(delta_hi_dim),
                                            threads=threads, psc=psc)
            if calculate_randomized:
                logging.debug(f"Correlation Calculation for negative control")
                delta_hi_dim_rndm = hi_dim_t_rndm - hi_dim
                self.corrcoef_random = colDeltaCorSqrt(hi_dim, np.sqrt(np.abs(delta_hi_dim_rndm) + psc) * np.sign(
                    delta_hi_dim_rndm), threads=threads, psc=psc)
        else:
            raise NotImplementedError(f"transform={transform} is not a valid parameter")
        np.fill_diagonal(self.corrcoef, 0)
        if calculate_randomized:
            np.fill_diagonal(self.corrcoef_random, 0)
    return self

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

@jit(nopython=True)
def permute_rows_nsign(A: np.ndarray) -> None:
    """Permute in place the entries and randomly switch the sign for each row of a matrix independently.
    """
    plmi = np.array([+1, -1])
    for i in range(A.shape[0]):
        np.random.shuffle(A[i, :])
        A[i, :] = A[i, :] * np.random.choice(plmi, size=A.shape[1])