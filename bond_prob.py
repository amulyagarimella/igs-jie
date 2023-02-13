from jie.jie.aligner import (log_bond)
from pandas import pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

# defaults
l_p_bp = 150.
nm_per_bp = .0004
pixel_dist = 100.

def calc_bond_prob (p1, p2, corr_fac=nm_per_bp/pixel_dist, l_p_bp=l_p_bp):
    '''
        p1, p2: should be in (x,y,z, bp) format
        corr_fac: scale genomic dist (bp) into pixels (e.g nm_per_bp / pixel_dist)
        l_p_bp: persistence length in bp
    '''
    # persistence length
    l_p = l_p_bp * corr_fac
    # ideal bond length is same as observed since only two points
    # observed distance between points
    observed = (p2[3] - p1[3]) * corr_fac
    return log_bond(l_p, observed, observed)

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
    probs_matrix = pd.DataFrame(index=row_range, columns=row_range)
    out_cols = ['x_um_1', 'y_um_1', 'z_um_1', 'hg38_bp_1', 
                'x_um_2', 'y_um_2', 'z_um_2', 'hg38_bp_2', 
                'prob', 'dist_bp', 'dist_bp_abs', 'dist_euc']
    probs_w_coords = pd.DataFrame(index=range(nrows**2), columns=out_cols)

    for i in row_range:
        for j in range(i,nrows):
            p1 = data.iloc[i]
            p2 = data.iloc[j]
            prob = calc_bond_prob(tuple(p1), tuple(p2))
            probs_matrix[i][j] = prob
            probs_matrix[j][i] = prob

            dist_bp = p2[3] - p1[3] 
            dist_euc = math.dist(p1[0:3], p2[0:3])
            out_data = list(p1) + list(p2) + [prob, dist_bp, abs(dist_bp), dist_euc]
            new_row = pd.Series(out_data, index=out_cols)
            # todo ensure this works
            probs_w_coords = pd.concat([probs_w_coords, new_row])

    return probs_matrix, probs_w_coords

def score_heatmap (probs_matrix):
    sns.heatmap(data=probs_matrix)
    sns.xlabel('log bond probabilities')
    plt.show()
    #plt.savefig('')
    plt.clf()

def score_hist (probs):
    sns.histplot(x=probs)
    sns.xlabel('log bond probabilities')
    plt.show()
    #plt.savefig('')
    plt.clf()

# score scatterplots
def score_v_x (data, x, xlabel=None):
    sns.scatter(data, x=x, y="prob")
    sns.ylabel('log bond probabilities')
    if xlabel is not None:
        sns.xlabel(xlabel)
    plt.show()
    #plt.savefig('')
    plt.clf()

# A scatter plot of bond score versus genomic distance
def score_v_distance_bp (probs):
    score_v_x(data=probs, x="dist_bp_abs")

# and another one of bond scores vs euclidian distance
def score_v_distance (probs):
    score_v_x(data=probs, x="dist_euc")

'''TODO: Plot all of the chr1 reads in xyz space and for any given point, color 
all other points by the bond score with that point (maybe make the given point
bigger or a completely different color so we can find it)'''
def visualize_reads (probs, p1):
    pass