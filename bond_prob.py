from jie.jie.aligner import log_bond
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
from collections import Counter

# defaults
l_p_bp = 150.
nm_per_bp = .0004
pixel_dist = 100.

# todo why is bond prob so high?!?
def calc_bond_prob (p1, p2, corr_fac=nm_per_bp/pixel_dist, l_p_bp=l_p_bp):
    '''
        p1, p2: should be in (x,y,z, bp) format
        corr_fac: scale genomic dist (bp) into pixels (e.g nm_per_bp / pixel_dist)
        l_p_bp: persistence length in bp
    '''
    # persistence length
    l_p = l_p_bp * corr_fac
    # observed distance between points
    ideal = abs(p2[3] - p1[3])
    observed = math.dist(p1[0:3], p2[0:3])
    return log_bond(l_p, ideal, observed)

# create heatmap of bond probs
def pairwise_bond_prob (data):
    '''
    input:
        data: must have cols [x_um, y_um, z_um, hg38_pos]
    output:
        probs_matrix: matrix of probabilities (dataframe)
        probs_list: dataframe of pairwise probabilities with corresponding points
                    and extra info (genomic and euclidean distance)
    '''
    nrows = len(data.index)
    row_range = range(nrows)
    print(data.head())
    probs_matrix = pd.DataFrame(index=data.hg38_pos, columns=data.hg38_pos, dtype=float)
    out_cols = ['x_um_1', 'y_um_1', 'z_um_1', 'hg38_bp_1', 
                'x_um_2', 'y_um_2', 'z_um_2', 'hg38_bp_2', 
                'prob', 'dist_bp', 'dist_bp_abs', 'dist_euc']
    probs_w_coords = pd.DataFrame(columns=out_cols, dtype=float)

    for i in row_range:
        for j in range(i+1,nrows):
            p1 = data.iloc[i]
            p2 = data.iloc[j]
            prob = calc_bond_prob(tuple(p1), tuple(p2))
            probs_matrix.iloc[i,j] = prob
            probs_matrix.iloc[j,i] = prob
            dist_bp = p2[3] - p1[3]
            dist_euc = math.dist(p1[0:3], p2[0:3])
            out_data = list(p1) + list(p2) + [prob, dist_bp, abs(dist_bp), dist_euc]
            new_row = pd.DataFrame(out_data, index=out_cols).T 
            probs_w_coords = pd.concat([probs_w_coords, new_row],axis=0)

    return probs_matrix, probs_w_coords

def score_heatmap (probs_matrix):
    sns.heatmap(data=probs_matrix.sort_index(axis=0).sort_index(axis=1))
    plt.title('Log bond probabilities ')
    plt.tight_layout()
    plt.show()
    #plt.savefig('')
    plt.clf()

def score_hist (probs):
    sns.histplot(x=probs['prob'])
    plt.xlabel('log bond probabilities')
    plt.show()
    #plt.savefig('')
    plt.clf()

# score scatterplots
def score_v_x (data, x, xlabel=None):
    sns.scatterplot(data, x=x, y="prob")
    plt.ylabel('log bond probabilities')
    plt.title(f'Log bond probabilites vs. {xlabel}')
    if xlabel is not None:
        plt.xlabel(xlabel)
    plt.show()
    #plt.savefig('')
    plt.clf()

# A scatter plot of bond score versus genomic distance
def score_v_distance_bp (probs):
    score_v_x(data=probs, x='dist_bp_abs', xlabel="Absolute genomic distance")

# and another one of bond scores vs euclidian distance
def score_v_distance (probs):
    score_v_x(data=probs, x='dist_euc', xlabel="Euclidean distance")

'''Plot all of the chr1 reads in xyz space and for any given point, color 
all other points by the bond score with that point (maybe make the given point
bigger or a completely different color so we can find it)'''
def visualize_reads (probs, p):
    """
    input:
        probs:
        p: (x, y, z, bp)
    """
    # todo weird column stuff is not ideal
    # isolate pairs that include p - todo does this work??
    p_pairs_1 = probs.loc[(probs.x_um_1 == p[0])&
                        (probs.y_um_1 == p[1])&
                        (probs.z_um_1 == p[2])].drop(columns=['x_um_1', 'y_um_1', 'z_um_1', 'hg38_bp_1'])
    print(p_pairs_1)
    p_pairs_2 = probs.loc[(probs.x_um_2 == p[0])&
                        (probs.y_um_2 == p[1])&
                        (probs.z_um_2 == p[2])].drop(columns=['x_um_2', 'y_um_2', 'z_um_2', 'hg38_bp_2'])
    new_cols = ['x_um', 'y_um', 'z_um', 'hg35_bp', 'prob', 'dist_bp', 'dist_bp_abs', 'dist_euc']
    p_pairs_1.columns = new_cols
    p_pairs_2.columns = new_cols
    p_pairs = pd.concat([p_pairs_1, p_pairs_2], axis=0, ignore_index=True)
    plot = plt.axes(projection='3d')
    plot.scatter(xs=p_pairs['x_um'], ys=p_pairs['y_um'], zs=p_pairs['z_um'], c=p_pairs['prob'])
    plot.scatter(xs=p[0], ys=p[1], zs=p[2], s=[50],c="black",marker="X")
    plt.title(f'3D visualization of chromosome points\nColored by distance from ({p[0]:.2f},{p[1]:.2f},{p[2]:.2f})')
    plt.show()
    plt.clf

    # sanity check - color by euc/gen distance TODO


def plot_chr (cell_id=1):
    # divide by cell id
    all_cells_file = 'Table_S1_pgp1_data_table.csv'
    all_cells = pd.read_csv(all_cells_file)
    cell_data = all_cells[all_cells.cell_id==cell_id]

    # make data processing edits - see walkthru 00
    cols_to_use = ['x_um', 'y_um', 'z_um', 'hg38_pos']
    new_cols = ['x_hat', 'y_hat', 'z_hat', 'hyb']
    """cell_pts_input = cell_data.copy(deep=True).loc[:, cols_to_use]
    cell_pts_input.columns = new_cols
    cell_pts_input.loc[:, 'sig_x'] = 1
    cell_pts_input.loc[:, 'sig_y'] = 1
    cell_pts_input.loc[:, 'sig_z'] = 1 
    cell_pts_input = cell_pts_input.sort_values(by="hyb",ascending=True)"""
    # gene_dist = sorted(list(set(cell_pts_input.hyb.tolist())))

    # feed cell into find_all_chr
    chrs = sorted(list(Counter(cell_data.hg38_chr).keys()))
    for chr in range(1,2):
        chr_data = cell_data.loc[cell_data.hg38_chr == chr, cols_to_use].reset_index(drop=True)
        """chr_data.x_um = chr_data.x_um.multiply(1000)
        chr_data.y_um = chr_data.y_um.multiply(1000)
        chr_data.z_um = chr_data.z_um.multiply(1000)"""
        probs_matrix, probs = pairwise_bond_prob(chr_data)
        score_heatmap(probs_matrix)
        score_hist(probs)
        score_v_distance(probs)
        score_v_distance_bp(probs)
        print(chr_data)
        visualize_reads(probs, chr_data.iloc[0])
    
plot_chr()