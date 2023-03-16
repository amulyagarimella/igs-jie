import ipywidgets as wg
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import igraph
import math
import copy

from mpl_toolkits.mplot3d import Axes3D
from IPython.display import Image
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


from jie.jie.demo import mock_experiment, error_types, polymer_model, polymer_skip, find_polymer

from jie.jie.aligner import (log_bond,
                         find_chr,
                         find_all_chr,
                         edge_weights,
                         boundary_init)

from jie.jie.utilities import (cartesian_esqsum, 
                           cartesian_sqdiff, 
                           cartesian_diff, 
                           check_lp_wgaps,
                           find_loci_dist)

from collections import Counter
import statistics

# make mock hyb based on SPATIAL distance
def join_pos (cell_chr, l_p_bp=150, nm_per_bp = .0004, pixel_dist=100.):
    l_p = l_p_bp * nm_per_bp / pixel_dist
    threshold = l_p*10
    new = cell_chr.copy(deep=True).sort_values(by="pos",ascending=True).reset_index()
    for i in range(len(new)-1):
        p1 = new.iloc[i,1:4]
        p2 = new.iloc[i+1,1:4]
        if math.dist(p1,p2) <= threshold:
            #m = (new.pos[i+1] + new.pos[i])/2
            new.iloc[i+1,1:4] = p1
            print("found close match")
            i+=1
        
    for i in range(len(new)-1):
        if new.pos[i+1] - new.pos[i] < l_p_bp*100:
            new.pos[i+1] = new.pos[i]
            i+=1
    
    return new
            
# 1500 kb together make same 
def find_binsize_genedist (cell_chr):
    # find bin size - median base pair interval between genomic loci
    # divide by chr
    ind = sorted(list(set(cell_chr.pos)))
    if len(ind) <= 1:
        raise ValueError
    # diffs = []
    diffs = [j-i for i, j in zip(ind[:-1], ind[1:])]
    binsize = statistics.median(diffs)
    genedist = [0]
    for i in range(len(diffs)):
        genedist.append(diffs[i] + genedist[i])
    genedist.append(ind[-1])
    return binsize, genedist

def realign_chr (embryo_id=1, cell_id=1, lines=True):
    matplotlib.rcParams.update({'font.size': 5})
    # divide by cell id
    all_cells_file = 'coord_table_57_embryos_pub_210423.csv'
    all_cells = pd.read_csv(all_cells_file)
    cell_data = all_cells[(all_cells.embryo_id==embryo_id)&(all_cells.cell_id==cell_id)]

    # make data processing edits - see walkthru 00
    cols_to_use = ['x_um_abs', 'y_um_abs', 'z_um_abs', 'pos']
    new_cols = ['x_hat', 'y_hat', 'z_hat', 'pos']
    # gene_dist = sorted(list(set(cell_pts_input.hyb.tolist())))
    chrs = sorted(list(Counter(cell_data.chr).keys()))
    fig = plt.figure()
    for chr in chrs:
        print(f"chr {chr}")
        chr_data = cell_data[cell_data.chr == chr]
        chr_data = chr_data.loc[:,cols_to_use]
        chr_data.columns = new_cols
        chr_data.loc[:, 'sig_x'] = 1
        chr_data.loc[:, 'sig_y'] = 1
        chr_data.loc[:, 'sig_z'] = 1 
        # chr_data = chr_data.sort_values(by="hyb",ascending=True)
        try:
            new = join_pos(chr_data)
            bin_size, gene_dist = find_binsize_genedist(new)
            pos_unique = list(Counter(new.pos).keys())
            # new['hyb'] = new['pos']
            new['hyb'] = [pos_unique.index(p) for p in new.pos]
            print(f"{len(pos_unique)} mock hybridization points")
            new.drop(columns=["pos"], inplace=True)
            res = find_all_chr(new,
                        gene_dist=gene_dist,
                        bin_size=bin_size,
                        nm_per_bp = .0004,
                        pixel_dist=100.,
                        num_skip=10)
            print(f"{len(res)} fiber(s)")
            print(res)
            # generate 3d viz
            colors = ['r','b','g','c','m','y']
            ax = fig.add_subplot(4,5,chr, projection='3d')
            if len(res) == 0:
                pass
            for i in range(len(res)):
                r = res[i]
                orig_data = chr_data.loc[list(r['index'])].sort_values(by="pos",ascending=True)
                ax.scatter(xs=orig_data['x_hat'], ys=orig_data['y_hat'], zs=orig_data['z_hat'], c=colors[i])
                if lines:
                    ax.plot(xs=orig_data['x_hat'], ys=orig_data['y_hat'], zs=orig_data['z_hat'], c=colors[i])
            # get data fw no match
            clustered = pd.concat(res)
            print(chr_data)
            unclustered = chr_data.loc[[x for x in chr_data.index if x not in list(clustered["index"])]].sort_values(by="pos",ascending=True)
            ax.scatter(xs=unclustered['x_hat'], ys=unclustered['y_hat'], zs=unclustered['z_hat'], c=colors[len(res)+1])
            ax.set_title(f"chr {chr}")
            print("\n")
        except Exception as e:
            print(f"wasn't able to find_chr for {chr}\n")
            print(e)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.suptitle(f"embryo {embryo_id}, cell {cell_id}")
    plt.tight_layout()
    plt.show()
    plt.clf()

realign_chr(embryo_id=51)

"""
- extract prob fn two points 
"""