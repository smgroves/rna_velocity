import numpy as np
import matplotlib.pyplot as plt
import velocyto as vcy
from velocyto.analysis import scatter_viz

## Plot in original dimensions (for low-dimensional data)
def plot_2_gene_arrows(vlm, gene_0 = 0,gene_1 = 1, which_S = 'Sx_sz', plot_ss = False):
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

def plot_arrows(vlm,filename = "grid_arrows",type = "field", quiver_scale = 1.5, marker_overlay = None):

    if type == "field":
        print("Plotting grid arrows...")
        plt.figure(None, (20, 10))
        vlm.plot_grid_arrows(quiver_scale=quiver_scale,
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
        vlm.plot_arrows_embedding(quiver_scale=1.4, scale_type='relative', #color_arrow=marker_overlay,
                                  plot_scatter=False, new_fig=False, plot_random=True)
                                  # scatter_kwargs_dict = {"cmap":plt.cm.get_cmap("viridis")})
        plt.savefig(f'{filename}.pdf')



def _plot_phase_portrait_smg(vlm, gene, gs_i, which_S = 'Sx_sz'):
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

def plot_by_cluster(vlm, h =4, l = 3, cluster_list = (["NON-NE", "NE", "NE-V1", "NE-V2"]), type = 'int'):

    plt.figure(None, (14, 14))
    gs = plt.GridSpec(h, l)
    for i, t in enumerate(range(max(vlm.ca["Clusters"]) + 1)):
        quiver_scale = 2.5
        if type == 'int':
            q = vlm.cluster_labels == str(t)
            treat_index = [i for i, x in enumerate(q) if x]
        if type == 'str':
            treat_index = np.where(vlm.ca['Clusters'] == str(t))[0]
        treat_index = np.random.choice(treat_index, size=len(treat_index), replace=False)

        plt.subplot(gs[i])
        plt.title(t, size='xx-large')

        plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1], c="k", alpha=0.2, s=2, edgecolor="")

        ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0] / 1.), replace=False)
        quiver_kwargs = dict(headaxislength=20, headlength=20, headwidth=20, linewidths=0.25, width=0.00045,
                             edgecolors="k",
                             color=vlm.colorandum[treat_index], alpha=0.7)
        plt.quiver(vlm.embedding[treat_index, 0], vlm.embedding[treat_index, 1],
                   vlm.delta_embedding[treat_index, 0], vlm.delta_embedding[treat_index, 1],
                   scale=quiver_scale, **quiver_kwargs)

    plt.tight_layout()
    plt.tight_layout()
    plt.savefig("clusterplot.pdf")


# vlm must have Sx_sz attribute and cluster_labels that are colored before running!
def plot_multigenes(glist, vlm):
    plt.figure(None, (17, 18), dpi=300)
    gs = plt.GridSpec(6, 6)
    for i, gn in enumerate(glist):
        print(gn)
        ax = plt.subplot(gs[i * 3])
        try:
            ix = np.where(vlm.ra["Gene"] == gn)[0][0]
            print(ix)
        except:
            continue
        vcy.scatter_viz(vlm.Sx_sz[ix, :], vlm.Ux_sz[ix, :], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
        plt.title(gn)
        xnew = np.linspace(0, vlm.Sx[ix, :].max())
        plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
        plt.ylim(0, np.max(vlm.Ux_sz[ix, :]) * 1.02)
        plt.xlim(0, np.max(vlm.Sx_sz[ix, :]) * 1.02)
        minimal_yticks(0, np.max(vlm.Ux_sz[ix, :]) * 1.02)
        minimal_xticks(0, np.max(vlm.Sx_sz[ix, :]) * 1.02)
        despline()

        vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i * 3 + 1], s=3, rasterized=True, which_tsne='embedding')

        vlm.plot_expression_as_color(gene_name=gn, gs=gs[i * 3 + 2], s=3, rasterized=True, which_tsne='embedding')

    plt.tight_layout()
    plt.show()


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