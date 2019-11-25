import numpy as np
import matplotlib.pyplot as plt
import velocyto as vcy
from velocyto.analysis import scatter_viz
import seaborn as sns
import logging
import matplotlib.patches as mpatches
from matplotlib import cm
## Plot in original dimensions (for low-dimensional data)
def plot_2_gene_arrows(vlm, gene_0 = 0,gene_1 = 1, which_S = 'Sx_sz', plot_ss = False, separate = False):
    if separate == False:
        if which_S == 'Sx_sz':
            arrowprops = dict(
                arrowstyle="->")
            x = vlm.Sx_sz[gene_0]
            y = vlm.Sx_sz[gene_1]
            fig, ax = plt.subplots()
            scatter = ax.scatter(x,y, c=vlm.colorandum, label=vlm.colorandum)
            legend1 = ax.legend(*scatter.legend_elements(),
                                loc="upper right", title="Time points")
            ax.add_artist(legend1)
            for i in range(len(x)):
                ax.annotate('',(vlm.Sx_sz_t[gene_0][i],vlm.Sx_sz_t[gene_1][i]), (x[i],y[i]), arrowprops = arrowprops)
            plt.xlabel(f"Gene {gene_0} Counts")
            plt.xlabel(f"Gene {gene_1} Counts")
            plt.title("Two Gene Model Spliced Counts")

        elif which_S == 'Sx':
            arrowprops = dict(
                arrowstyle="->")
            x = vlm.Sx[gene_0]
            y = vlm.Sx[gene_1]
            fig, ax = plt.subplots()
            scatter = ax.scatter(x,y, c = vlm.colorandum, label = vlm.colorandum)

            for i in range(len(x)):
                ax.annotate('',(vlm.Sx_t[gene_0][i],vlm.Sx_t[gene_1][i]), (x[i],y[i]), arrowprops = arrowprops)
                plt.tight_layout()
            plt.xlabel(f"Gene {gene_0} Counts")
            plt.xlabel(f"Gene {gene_1} Counts")
            plt.title("Two Gene Model Spliced Counts")
    else:
        tps = vlm.ca['SampleID']
        for tp in (np.unique(tps)):
            arrowprops = dict(
                arrowstyle="->")
            x = vlm.Sx[gene_0][vlm.ca['SampleID'] == tp]
            y = vlm.Sx[gene_1][vlm.ca['SampleID'] == tp]
            x_t = vlm.Sx_t[gene_0][vlm.ca['SampleID'] == tp]
            y_t = vlm.Sx_t[gene_1][vlm.ca['SampleID'] == tp]
            fig, ax = plt.subplots()
            ax.scatter(x, y)
            for i in range(len(x)):
                ax.annotate('', (x_t[i], y_t[i]), (x[i], y[i]), arrowprops=arrowprops)
            plt.xlabel(f"Gene {gene_0} Counts")
            plt.xlabel(f"Gene {gene_1} Counts")
            plt.xlim([-1, np.max(vlm.Sx[gene_0])])
            plt.ylim([-1, np.max(vlm.Sx[gene_1])])
    plt.show()

def plot_1_gene_arrows(vlm, gene_0 = 0, which_S = 'Sx_sz'):
    if which_S == 'Sx_sz':
        arrowprops = dict(
            arrowstyle="->")
        y = vlm.Sx_sz[gene_0]
        x = [int(i) for i in vlm.ca['SampleID']]
        fig, ax = plt.subplots()
        scatter = ax.scatter(x,y, c=vlm.colorandum, label=vlm.colorandum)
        legend1 = ax.legend(*scatter.legend_elements(),
                            loc="upper right", title="Time points")
        ax.add_artist(legend1)
        for i in range(len(x)):
            ax.annotate('',(x[i]+1,vlm.Sx_sz_t[gene_0][i]), (x[i],y[i]), arrowprops = arrowprops)
        plt.xlabel(f"Gene {gene_0} Counts")
        plt.xlabel(f"Time Points")
        plt.title("One Gene Model Spliced Counts")

        plt.show()
    elif which_S == 'Sx':
        arrowprops = dict(
            arrowstyle="->")
        y = vlm.Sx[gene_0]
        x = [int(i) for i in vlm.ca['SampleID']]
        fig, ax = plt.subplots()
        scatter = ax.scatter(x,y, c = vlm.colorandum, label = vlm.colorandum)

        for i in range(len(x)):
            ax.annotate('',(x[i]+1,vlm.Sx_t[gene_0][i]), (x[i],y[i]), arrowprops = arrowprops)
            plt.tight_layout()
        plt.xlabel(f"Gene {gene_0} Counts")
        plt.xlabel(f"Time Points")
        plt.title("One Gene Model Spliced Counts")
        plt.show()

def hist_last_tp(vlm, gene_0, tp = None):
        if tp == None: last_tp = vlm.Sx[gene_0]
        else:
            last_tp = vlm.Sx[gene_0][vlm.ca['SampleID'] == tp]
        plt.figure()
        plt.hist(x = last_tp, bins = 50)
        plt.show()

def plot_distr_over_time(vlm, attractors, gene_0 = 0, gene_1 = 1):
    tps = vlm.ca['SampleID']
    sqrtn = int(np.ceil(np.sqrt(len(np.unique(tps)))))
    fig = plt.figure(figsize=(20,20))
    gs = plt.GridSpec(sqrtn, int(np.ceil(len(np.unique(tps)) / sqrtn)))
    for i, tp in enumerate(np.unique(tps)):
        x_sub = vlm.Sx[gene_0][vlm.ca['SampleID'] == tp]
        y_sub = vlm.Sx[gene_1][vlm.ca['SampleID'] == tp]
        plt.subplot(gs[i])
        sns.kdeplot(data=x_sub, data2=y_sub, shade=True, shade_lowest=False, cmap = 'Blues')
        plt.scatter(x = x_sub, y = y_sub, c="w", s=30, linewidth=1, marker="+")
        plt.scatter(x = attractors[0], y = attractors[1], c="r", s=50, linewidth=3, marker="+")
        plt.xlim([-1, np.max(vlm.Sx[gene_0])])
        plt.ylim([-1, np.max(vlm.Sx[gene_1])])
    plt.show()

def plot_arrows(vlm,filename = "grid_arrows",type = "field", quiver_scale = 1.5, marker_overlay = None):

    if type == "field":
        print("Plotting grid arrows...")
        plt.figure(None, (20, 10))
        plot_grid_arrows_smg(vlm, quiver_scale=quiver_scale,
                             scatter_kwargs_dict={"alpha": 0.35, "lw": 0.35, "edgecolor": "0.4", "s": 20,
                                                  "rasterized": True}, min_mass=1.5, angles='xy', scale_units='xy',
                             headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,scale_type='absolute',
                             plot_random=False)

        plt.show()
        plt.savefig(f'{filename}.pdf')
    if type == "field_mask":
        print("Plotting grid arrows...")
        plt.figure(None, (20, 10))
        vlm.plot_grid_arrows(quiver_scale=quiver_scale,
                             scatter_kwargs_dict={"alpha": 0.35, "lw": 0.35, "edgecolor": "0.4", "s": 20,
                                                  "rasterized": True, "cmap":plt.cm.get_cmap('viridis'), "c":marker_overlay},
                             min_mass=1.5, angles='xy', scale_units='xy',
                             headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                             plot_random=True)
        # plt.show()
        plt.savefig(f'{filename}.pdf')
    if type == "embedded":
        print("Plotting embedded arrows...")
        plot_arrows_embedding_smg(vlm,choice = vlm.S.shape[1], quiver_scale=1.4, scale_type='absolute',
                                  #color_arrow=marker_overlay,
                                  plot_scatter=False, new_fig=False, plot_random=False)
                                  # scatter_kwargs_dict = {"cmap":plt.cm.get_cmap("viridis")})
        plt.show()
        # plt.savefig(f'{filename}.pdf')

def plot_phase_portraits(vlm, genes, which_S = 'Sx_sz'):
        """Plot spliced-unspliced scatterplots resembling phase portraits

        Arguments
        ---------
        genes: List[str]
            A list of gene symbols.
        """
        n = len(genes)
        sqrtn = int(np.ceil(np.sqrt(n)))
        gs = plt.GridSpec(sqrtn, int(np.ceil(n / sqrtn)))
        for i, gn in enumerate(genes):
            _plot_phase_portrait_smg(vlm,gn, gs[i], which_S)
        plt.show()

#######################################################################################################################
# Changes to original plotting functions
def _plot_phase_portrait_smg(vlm, gene, gs_i, which_S='Sx_sz'):
    """Plot spliced-unspliced scatterplot resembling phase portrait
    """
    if gene is None:
        plt.subplot(111)
    else:
        plt.subplot(gs_i)
    ix = np.where(vlm.ra["Gene"] == gene)[0][0]
    if which_S == 'Sx_sz':
        scatter_viz(vlm.Sx_sz[ix, :], vlm.Ux_sz[ix, :], c=vlm.colorandum, s=5, alpha=0.4)
        plt.title(gene)
        xnew = np.linspace(0, vlm.Sx_sz[ix, :].max())
        plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
    elif which_S == 'Sx':
        scatter_viz(vlm.Sx[ix, :], vlm.Ux[ix, :], c=vlm.colorandum, s=5, alpha=0.4)
        plt.title(gene)
        xnew = np.linspace(0, vlm.Sx[ix, :].max())
        plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")

def plot_arrows_embedding_smg(self, choice = "auto", quiver_scale = "auto",color_arrow = 'cluster',
                          scale_type = "relative", new_fig = False, plot_random = True,
                          scatter_kwargs = {}, plot_scatter = False, **quiver_kwargs):
    """Plots velocity on the embedding cell-wise

    Arguments
    ---------
    choice: int, default = "auto"
        the number of cells to randomly pick to plot the arrows (To avoid overcrowding)
    quiver_scale: float, default="auto"
        Rescaling factor applied to the arrow field to enhance visibility
        If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
        If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
        NOTE: Despite a similar option than plot_grid_arrows, here there is no strong motivation to calculate the scale relative to the randomized control
        This is because the randomized doesn't have to have smaller velocity cell-wise, there might be for example
        scattered cells that will have strong velocity but they will, correctly just average out when calculating the average velocity field.
    scale_type: str, default="relative"
        How to interpret `quiver_scale`:
        If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
        (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
        If "absolute" the value will be passed to the matplotlib quiver function
    plot_scatter: bool, default = False
        whether to plot the points
    scatter_kwargs: Dict
        A dictionary containing all the keywords arguments to pass to matplotlib scatter
        by default the following are passed: c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3. But they can be overridden.
    color_arrow: str, default = "cluster"
        the color of the arrows, if "cluster" the arrows are colored the same as the cluster
    new_fig: bool, default=False
        whether to create a new figure
    plot_random: bool, default=True
        whether to plot the randomized control next to the plot
    **quiver_kwargs: dict
        keyword arguments to pass to quiver
        By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.

    Returns
    -------
    Nothing, just plots the tsne with arrows
    """
    if choice == "auto":
        choice = int(self.S.shape[1] / 3)
        logging.warning(
            f"Only {choice} arrows will be shown to avoid overcrowding, you can choose the exact number setting the `choice` argument")
    _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
    _scatter_kwargs = dict(c="0.8", alpha=0.4, s=10, edgecolor=(0, 0, 0, 1), lw=0.3)
    _scatter_kwargs.update(scatter_kwargs)
    if new_fig:
        if plot_random and hasattr(self, "delta_embedding_random"):
            plt.figure(figsize=(22, 12))
        else:
            plt.figure(figsize=(14, 14))

    ix_choice = np.random.choice(self.embedding.shape[0], size=choice, replace=False)

    # Determine quiver scale
    if scale_type == "relative":
        if hasattr(self, "delta_embedding_random"):
            plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0),
                                        2)  # Diagonal of the plot
            arrows_scale = np.percentile(np.linalg.norm(self.delta_embedding_random, 2, 1),
                                         80)  # Tipical length of an arrow
            if quiver_scale == "auto":
                quiver_scale = arrows_scale / (plot_scale * 0.005)
            else:
                quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.005)
        else:
            raise ValueError("""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
            Please run estimate_transition_prob or set `scale_type` to `absolute`""")
    else:
        logging.warning(
            "The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

    if color_arrow == "cluster":
        colorandum = self.colorandum[ix_choice, :]
    else:
        colorandum = color_arrow

    _quiver_kwargs.update({"color": colorandum})
    _quiver_kwargs.update(quiver_kwargs)

    if plot_random and hasattr(self, "delta_embedding_random"):
        plt.subplot(122)
        plt.title("Randomized")
        if plot_scatter:
            plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)
        plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
                   self.delta_embedding_random[ix_choice, 0], self.delta_embedding_random[ix_choice, 1],
                   scale=quiver_scale, **_quiver_kwargs)
        plt.axis("off")
        plt.subplot(121)
        plt.title("Data")

    if plot_scatter:
        plt.scatter(self.embedding[:, 0], self.embedding[:, 1], **_scatter_kwargs)

    plt.quiver(self.embedding[ix_choice, 0], self.embedding[ix_choice, 1],
               self.delta_embedding[ix_choice, 0], self.delta_embedding[ix_choice, 1],
               scale=quiver_scale, **_quiver_kwargs)
    patches = []
    for i in range(len(np.unique(self.ca['SampleID']))):
        patches.append(mpatches.Patch(color=cm.get_cmap('Set2')(i), label=f'Time {i+1}'))
    plt.legend(handles=patches)

def plot_grid_arrows_smg(self, quiver_scale = "auto", scale_type = "relative",
                     min_mass = 1, min_magnitude = None,
                     scatter_kwargs_dict = None, plot_dots= False, plot_random = False,
                     **quiver_kwargs):
    """Plots vector field averaging velocity vectors on a grid

    Arguments
    ---------
    quiver_scale: float, default="auto"
        Rescaling factor applied to the arrow field to enhance visibility
        If "auto" the scale is selected using the randomized (negative) control (even if `plot_random`=False)
        If a float is provided the interpretation of the value depends on the parameter `scale_type`, see below.
        NOTE: In the situation where "auto" is resulting in very small or big velocities, pass a float to this parameter
        The float will be interpreted as a scaling, importantly both the data and the control will be scaled
        in this way you can rescale the velocity arbitrarily without the risk of observing just an overfit of the noise
    scale_type: str, default="relative"
        How to interpret `quiver_scale`:
        If "relative" (default) the value will be used as a scaling factor and multiplied by the value from "auto"
        (it follows that quiver_scale="auto" is equivalent to quiver_scale=1)
        If "absolute" the value will be passed to the matplotlib quiver function (not recommended if you are not sure what this implies)
    min_mass: float, default=1
        the minimum density around a grid point for it to be considered and plotted
    min_magnitude: float, default=None
        the minimum magnitude of the velocity for it to be considered and plotted
    scatter_kwargs_dict: dict, default=None
        a dictionary of keyword arguments to pass to scatter
        by default the following are passed: s=20, zorder=-1, alpha=0.2, lw=0, c=self.colorandum. But they can be overridden.
    plot_dots: bool, default=True
        whether to plot dots in correspondence of all low velocity grid points
    plot_random: bool, default=True
        whether to plot the randomized control next to the plot
    **quiver_kwargs: dict
        keyword arguments to pass to quiver
        By default the following are passed angles='xy', scale_units='xy', minlength=1.5. But they can be overridden.
    """
    # plt.figure(figsize=(10, 10))
    _quiver_kwargs = {"angles": 'xy', "scale_units": 'xy', "minlength": 1.5}
    _quiver_kwargs.update(quiver_kwargs)

    scatter_dict = {"s": 20, "zorder": -1, "alpha": 0.2, "lw": 0, "c": self.colorandum}
    if scatter_kwargs_dict is not None:
        scatter_dict.update(scatter_kwargs_dict)

    # Determine quiver scale
    if scale_type == "relative":
        if hasattr(self, "flow_rndm"):
            plot_scale = np.linalg.norm(np.max(self.flow_grid, 0) - np.min(self.flow_grid, 0),
                                        2)  # Diagonal of the plot
            arrows_scale = np.percentile(np.linalg.norm(self.flow_rndm[self.total_p_mass >= min_mass, :], 2, 1),
                                         90)  # Tipical lenght of an arrow
            if quiver_scale == "auto":
                quiver_scale = arrows_scale / (plot_scale * 0.0025)
            else:
                quiver_scale = quiver_scale * arrows_scale / (plot_scale * 0.0025)
        else:
            raise ValueError(""""`scale_type` was set to 'relative' but the randomized control was not computed when running estimate_transition_prob
            Please run estimate_transition_prob or set `scale_type` to `absolute`""")
    else:
        logging.warning(
            "The arrow scale was set to be 'absolute' make sure you know how to properly interpret the plots")

    mass_filter = self.total_p_mass < min_mass
    if min_magnitude is None:
        XY, UV = np.copy(self.flow_grid), np.copy(self.flow)
        if not plot_dots:
            UV = UV[~mass_filter, :]
            XY = XY[~mass_filter, :]
        else:
            UV[mass_filter, :] = 0
    else:
        XY, UV = np.copy(self.flow_grid), np.copy(self.flow_norm)
        if not plot_dots:
            UV = UV[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
            XY = XY[~(mass_filter | (self.flow_norm_magnitude < min_magnitude)), :]
        else:
            UV[mass_filter | (self.flow_norm_magnitude < min_magnitude), :] = 0

    if plot_random:
        if min_magnitude is None:
            XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_rndm)
            if not plot_dots:
                UV_rndm = UV_rndm[~mass_filter, :]
                XY = XY[~mass_filter, :]
            else:
                UV_rndm[mass_filter, :] = 0
        else:
            XY, UV_rndm = np.copy(self.flow_grid), np.copy(self.flow_norm_rndm)
            if not plot_dots:
                UV_rndm = UV_rndm[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
                XY = XY[~(mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude)), :]
            else:
                UV_rndm[mass_filter | (self.flow_norm_magnitude_rndm < min_magnitude), :] = 0

        plt.subplot(122)
        plt.title("Randomized")
        plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
        plt.quiver(XY[:, 0], XY[:, 1], UV_rndm[:, 0], UV_rndm[:, 1],
                   scale=quiver_scale, zorder=20000, **_quiver_kwargs)
        plt.axis("off")
        plt.subplot(121)
        plt.title("Data")

    plt.scatter(self.flow_embedding[:, 0], self.flow_embedding[:, 1], **scatter_dict)
    plt.quiver(XY[:, 0], XY[:, 1], UV[:, 0], UV[:, 1],
               scale=quiver_scale, zorder=20000, **_quiver_kwargs)
#######################################################################################################################
# Plotting utility functions

def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')


def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end)) + 1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""] * len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)


def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end)) + 1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""] * len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)

def gaussian_kernel(X, mu = 0, sigma=1):
    return np.exp(-(X - mu)**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)