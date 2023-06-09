import pickle
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, inconsistent
import matplotlib.pyplot as plt
import itertools
import copy
import seaborn as sns
import os
from matplotlib.colors import ListedColormap


# cloncalcell_cell_allpairs_ppca_20.pkl

cell_cell_file = 'clonal/cell_cell_allpairs_ppca_newID_20.pkl'
with open(cell_cell_file, 'rb') as f:
    cell_cell_20 = pickle.load(f)[:, :, 0]

cell_chr_fiber_file = "clonal/E14_clonal_24hr_coord_table_with_cluster.csv"
cell_chr_fiber_file_icloud = "clonal/.E14_clonal_24hr_coord_table_with_cluster.csv.icloud"
os.system(f"brctl download {cell_chr_fiber_file_icloud}")
cell_chr_fiber = pd.read_csv(cell_chr_fiber_file)
cell_chr_fiber['newCellID'] = cell_chr_fiber["fov"].astype(str) + cell_chr_fiber["cellID"].apply(lambda a : str(a).zfill(2))
cell_chr_fiber['newCellID'] = cell_chr_fiber['newCellID'].astype(int)

"""
with open('clonal/cell_cell_allpairs_ppca_newID_100.pkl', 'rb') as f:
    cell_cell_100 = pickle.load(f)[:, :, 0]
"""
clonal_info = pd.read_csv("clonal/E14_clonal_comb_rna_table.csv")
clonal_info['newCellID'] = clonal_info["fov"].astype(str) + clonal_info["cellID"].apply(lambda a : str(a).zfill(2))
clonal_info['newCellID'] = clonal_info['newCellID'].astype(int)
clonal_info_24hr = clonal_info.loc[clonal_info.expt == 'E14_clonal_24hr']
gfp_threshold = 5
gfp_expressing_idx = np.where(clonal_info_24hr.mEGFP>gfp_threshold)[0]

#non_gfp_expressing = set(clonal_info_24hr.loc[clonal_info_24hr.mEGFP<=gfp_threshold,'newCellID'])
a1 = np.tile(clonal_info_24hr.mEGFP<=gfp_threshold, (len(set(clonal_info_24hr.newCellID)),1))
a2 = a1.T
same_clone = a1 == a2

# determine accuracy
def acc (clusters):
    c1 = np.tile(clusters, (len(clusters),1))
    clust_clone = c1 == c1.T
    agree = same_clone == clust_clone
    acc = agree.sum()/agree.size
    return agree, acc

def cluster_fdr (cluster_idxs):
    agree = np.array([all([same_clone[i][j] for i,j in itertools.combinations(c_idx, 2)]) for c_idx in cluster_idxs])
    cluster_idxs = np.array(cluster_idxs, dtype=object)
    sisters = cluster_idxs[agree]
    false_sisters = cluster_idxs[~agree]
    return 1 - np.count_nonzero(agree)/np.size(agree), sisters, false_sisters

def all_cluster_pair_distances (cluster_idxs, distance):
    distances_allclust = []
    for c_idx in cluster_idxs:
        pairs = itertools.combinations(c_idx, 2)
        distances_allclust = distances_allclust + [[i,j,distance[i][j]] for i,j in pairs]
    return np.array(distances_allclust)


# hierarchical clustering
def hierarchical_cluster (cell_cell):
    linkage_cell_cell = linkage(cell_cell)
    # linkage_allfov_100 = linkage(cell_cell_100)
    inconsistencies = inconsistent(linkage_cell_cell)[:,3]
    distances = linkage_cell_cell[:,2]
    # decide whether to use inconsistency or distance here
    thresholds = np.unique(distances)
    threshold_fdr = np.vstack([thresholds, np.zeros(len(thresholds)), np.zeros(len(thresholds)), np.zeros(len(thresholds))])
    threshold_fdr[1,:] = np.nan
    threshold_fdr = np.array(threshold_fdr, dtype=object)

    for i in range(len(thresholds)):
        t = thresholds[i]
        clusters = fcluster(linkage_cell_cell, t=t, criterion='distance')
        cluster_set = set(clusters)
        # get indices of cells in each cluster
        cluster_idxs = [np.where(clusters == c)[0] for c in cluster_set]
        cluster_idxs_len = np.array([len(c_idx) for c_idx in cluster_idxs])
        if np.all(cluster_idxs_len == 1) or len(cluster_set) == 1:
            continue
        # FDR, sisters, false sisters
        fdr, sisters, false_sisters = cluster_fdr(cluster_idxs)
        distances_sisters = all_cluster_pair_distances(sisters, cell_cell)
        distances_false_sisters =  all_cluster_pair_distances(false_sisters, cell_cell)
        if len(distances_sisters) < 1:
            continue
        threshold_fdr[1,i] = fdr
        threshold_fdr[2,i] = distances_sisters
        threshold_fdr[3,i] = distances_false_sisters
        #print(all_cluster_pair_distances(cluster_idxs, cell_cell))
        
    min_fdr = min(threshold_fdr[1,np.where(threshold_fdr[1,:] > 0)[0]])
    min_fdr_full = threshold_fdr[:,np.where(threshold_fdr[1,:] == min_fdr)[0]]
    t = min_fdr_full[0]
    dendrogram(linkage_cell_cell)
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    label_colors = {idx : 'g' if int(idx.get_text()) in gfp_expressing_idx else 'k' for idx in xlbls}
    for lbl in xlbls:
        lbl.set_color(label_colors[lbl])
    plt.axhline(t, color='k', ls='--')
    plt.show()
    plt.clf()
    distances_sisters = min_fdr_full[2]
    # TODO
    print("distances sisters")
    print(distances_sisters)
    print(distances_sisters[0])
    #closest_idx = np.where(distances_sisters[:,3] == min(distances_sisters[:,3]))
    #print(closest_idx)
    #farthest_idx = np.where(distances_sisters[:,3] == max(distances_sisters[:,3]))
    two_cell_viz(cell_chr_fiber,distances_sisters[0][0][0], distances_sisters[0][0][1])
    distances_false_sisters = min_fdr_full[2]
    """if len(distances_false_sisters) > 1:
        closest_idx = np.where(distances_false_sisters[:,3] == min(distances_false_sisters[:,3]))
        print(closest_idx)
        farthest_idx = np.where(distances_false_sisters[:,3] == max(distances_false_sisters[:,3]))
        two_cell_viz(cell_chr_fiber,distances_false_sisters[closest_idx][0], distances_false_sisters[closest_idx][1])"""
        
    
    #print(f"min_cutoff: {min_fdr_t}")
    # get cluster_idxs
    # get all cluster_idxs that are sisters and not sisters
    # pick sister clusters with highest and lowest distances
    

def two_cell_viz (full_cell_data, cellid1, cellid2, x_colname='x',y_colname='y',z_colname='z',cellid_colname='newCellID',chrom_colname='chromID'):
    # TODO sharex and sharey but idt nec
    data1 = np.array(full_cell_data.loc[full_cell_data[cellid_colname] == cellid1, [x_colname, y_colname, chrom_colname]])
    data2 = np.array(full_cell_data.loc[full_cell_data[cellid_colname] == cellid2, [x_colname, y_colname, chrom_colname]])
    chrom = np.sort(np.intersect1d(data1[:,2], data2[:,2]))
    nsubplots = 2*(np.size(chrom) + 1)
    nrow = int(np.floor(np.sqrt(nsubplots)))
    ncol = nsubplots//nrow + 1
    # TODO make two belowlines more efficient?
    data1 = data1[np.isin(data1[:,2],chrom)]
    data2 = data2[np.isin(data2[:,2],chrom)]
    # fig, axs = plt.subplots(nrow, ncol)
    for i in range(nsubplots):
        ax = plt.subplot(nrow, ncol, i + 1)
        # shallow copy to save on space
        if i % 2 == 0:
            data = data1.copy()
        else:
            data = data2.copy()
        if i <= 1:
            sns.scatterplot(x=data[:,0],y=data[:,1],hue=data[:,2],ax=ax)
        else:  
            # cmap = ListedColormap(colors)
            # add a col describinb whether it is the chrom we want
            chrom_idx = int(np.floor((i-2)/2))
            sns.scatterplot(x=data[:,0],y=data[:,1],hue=data[:,2]==chrom[chrom_idx],ax=ax, size=0.5, palette = {True: sns.color_palette()[np.mod(chrom_idx, len(sns.color_palette()))], False: (0.71, 0.71, 0.71, 0.05)}, edgecolors='face')
        
        ax.get_legend().remove()
        plt.axis('off')
        # For even i use data1, odd i use data2
        # for first two i seaborn hue = chrom
        # for rest, figure out chr number and plot separately with rest grey OR do that in seaborn without plotting differently
    plt.show()


# TODO plot by fov
for fov in set(clonal_info_24hr["fov"]):
    ids = clonal_info_24hr.loc[clonal_info_24hr.fov==fov, "newCellID"]
    ids_in_total = np.isin(list(clonal_info_24hr.newCellID), list(ids))
    ids_idx = np.where(ids_in_total)[0]
    cell_cell_fov = cell_cell_20[ids_idx,:][:,ids_idx]
    hierarchical_cluster(cell_cell_fov)


# per fov?
