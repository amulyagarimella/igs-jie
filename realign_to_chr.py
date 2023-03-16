import ipywidgets as wg
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import igraph
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

def find_binsize_genedist (cell_chr):
    # find bin size - median base pair interval between genomic loci
    # divide by chr
    ind = cell_chr.hyb.sort_values(ascending=True).tolist()
    if len(ind) > 1:
        diffs = [ind[i+1] - ind[i] for i in range(len(ind)-1)]
        # print(diffs)
        binsize = statistics.median(diffs)
        genedist = [0]
        for i in range(len(diffs)):
            genedist.append(diffs[i] + genedist[i])
        genedist.append(ind[-1])
        return binsize, genedist
    else:
        raise ValueError(ind)


def realign_chr (cell_id=1):
    # divide by cell id
    all_cells_file = 'Table_S1_pgp1_data_table.csv'
    all_cells = pd.read_csv(all_cells_file)
    cell_data = all_cells[all_cells.cell_id==cell_id]

    # make data processing edits - see walkthru 00
    cols_to_use = ['x_um', 'y_um', 'z_um', 'hg38_pos']
    new_cols = ['x_hat', 'y_hat', 'z_hat', 'hyb']
    """
    cell_pts_input = cell_data.copy(deep=True).loc[:, cols_to_use]
    cell_pts_input.columns = new_cols
    cell_pts_input.loc[:, 'sig_x'] = 1
    cell_pts_input.loc[:, 'sig_y'] = 1
    cell_pts_input.loc[:, 'sig_z'] = 1 
    cell_pts_input = cell_pts_input.sort_values(by="hyb",ascending=True)
    """
    # gene_dist = sorted(list(set(cell_pts_input.hyb.tolist())))

    # feed cell into find_all_chr
    for chr in sorted(list(Counter(cell_data.hg38_chr).keys())):
        print(f"chr {chr}")
        chr_data = cell_data[cell_data.hg38_chr == chr]
        #print(chr_data.head())
        chr_data = chr_data.loc[:,cols_to_use]
        chr_data.columns = new_cols
        chr_data.loc[:, 'sig_x'] = 1
        chr_data.loc[:, 'sig_y'] = 1
        chr_data.loc[:, 'sig_z'] = 1 
        # chr_data = chr_data.sort_values(by="hyb",ascending=True)

        try:
            bin_size, gene_dist = find_binsize_genedist(chr_data)
            chr_data['hyb'] = range(len(chr_data))
            res = find_all_chr(chr_data,
                        gene_dist=gene_dist,
                        bin_size=bin_size)
            print(f"{len(res)} fiber(s)")
            print(res)
        except:
            print(f"wasn't able to find_chr for {chr}")

realign_chr()

"""
- extract prob fn two points 
"""