import ipywidgets as wg
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import igraph
import copy

from mpl_toolkits.mplot3d import Axes3D
from IPython.display import Image

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

def find_bin_size (cell_data):
    # find bin size - median base pair interval between genomic loci
    # divide by chr
    loci_ind_diffs = []
    chr_diffs = []
    for i in list(Counter(cell_data.hg38_chr).keys()):
        cell_chr = cell_data[cell_data.hg38_chr == i]
        ind = cell_chr.hg38_ind.sort_values(ascending=False).tolist()
        if len(ind) > 1:
            chr_diffs = [ind[i] - ind[i+1] for i in range(len(ind)-1)]
            loci_ind_diffs = loci_ind_diffs + chr_diffs
        else:
            continue
    return chr_diffs, statistics.median(loci_ind_diffs)

def find_gene_dist ():
    # find gene dist - reference genomic distances between locis imaged on given chr

    pass

def realign_chr (cell_id=1):
    # divide by cell id
    all_cells_file = 'Table_S1_pgp1_data_table.csv'
    all_cells = pd.read_csv(all_cells_file)
    cell_data = all_cells[all_cells.cell_id==cell_id]

    # TODO - idk how to do gene dist
    _, bin_size = find_bin_size(cell_data)

    # make data processing edits - see walkthru 00
    cols_to_use = ['x_um', 'y_um', 'z_um', 'hg38_ind']
    new_cols = ['x_hat', 'y_hat', 'z_hat', 'hyb']
    cell_pts_input = cell_data.copy(deep=True).loc[:, cols_to_use]
    cell_pts_input.columns = new_cols
    cell_pts_input.loc[:, 'sig_x'] = 1
    cell_pts_input.loc[:, 'sig_y'] = 1
    cell_pts_input.loc[:, 'sig_z'] = 1 
    cell_pts_input = cell_pts_input.sort_values(by="hyb",ascending=True)
    gene_dist = sorted(list(set(cell_pts_input.hyb.tolist())))

    # feed cell into find_all_chr
    
    
    """ValueError: Spatial distance estimated from genomic intervals separating
    loci (contour length) must be >> l_p (pixel dist). Double check the input
    persistence length (l_p_bp) and reference genome intervals (gene_dist)."""
    # if not np.all([elem/l_p >= 10 for elem in gdintervals]):
    return find_all_chr(cell_pts_input,
                        gene_dist=gene_dist,
                        bin_size=bin_size,
                        nm_per_bp = 500, #PROBABLY WRONG
                        pixel_dist=1000)
    # profit


print(realign_chr())