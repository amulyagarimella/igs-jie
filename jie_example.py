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
                         edge_weights,
                         boundary_init)

from jie.jie.utilities import (cartesian_esqsum, 
                           cartesian_sqdiff, 
                           cartesian_diff, 
                           check_lp_wgaps,
                           find_loci_dist)

# plot
fig = plt.figure(figsize = (18, 9))
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
# load data
df = pd.read_csv("../data/df_spalign_demo.csv")
df.columns = ['x_hat', 'y_hat', 'z_hat', 'Annot', 't', 'bp']
df.reset_index(drop = True, inplace = True)

# subset data
df_fiber = df[df['Annot'] == 'Fiber']
df_data = df[(df['Annot'] != 'Fiber') & (df['Annot'] != 'FN')].sort_values(by='t')

# add indexing
df_list = []
for hyb, group in df_data.groupby(by='t'):
    _df = copy.deepcopy(group)
    _df['i'] = list(range(_df.shape[0]))
    df_list.append(_df)
df_data = pd.concat(df_list)

df_data = df_data[['x_hat', 'y_hat', 'z_hat', 'Annot', 't', 'i', 'bp']]
    
# add localization err
df_data.loc[:, 'sig_x'] = 1
df_data.loc[:, 'sig_y'] = 1
df_data.loc[:, 'sig_z'] = 1  

# make all positions positive
xmin = np.min(df_data['x_hat'].values)
ymin = np.min(df_data['y_hat'].values)
zmin = np.min(df_data['z_hat'].values)

df_data.loc[:, 'x_hat'] = df_data['x_hat'].values - xmin
df_data.loc[:, 'y_hat'] = df_data['y_hat'].values - ymin
df_data.loc[:, 'z_hat'] = df_data['z_hat'].values - zmin

# translate fiber accordingly (for visualization)
df_fiber.loc[:, 'x_hat'] = df_fiber['x_hat'].values - xmin
df_fiber.loc[:, 'y_hat'] = df_fiber['y_hat'].values - ymin
df_fiber.loc[:, 'z_hat'] = df_fiber['z_hat'].values - zmin


# Draw imaged loci
cdict = {0:'Red', 1:'Green', 2:'Blue', 3:'Purple'}
ax1.scatter(df_data.loc[:, 'z_hat'],
           df_data.loc[:, 'y_hat'],
           df_data.loc[:, 'x_hat'],
           c=[cdict[elem] for elem in df_data.t], alpha = 1, marker='^', s = 100)

ax1.set_title("Input", fontsize = 20);

# Draw imaged loci
cdict = {0:'Red', 1:'Green', 2:'Blue', 3:'Purple'}
ax2.scatter(df_data.loc[:, 'z_hat'],
           df_data.loc[:, 'y_hat'],
           df_data.loc[:, 'x_hat'],
           c=[cdict[elem] for elem in df_data.t], alpha = 1, marker='^', s = 100)

# Draw fiber
ax2.plot(df_fiber.loc[:, 'z_hat'], df_fiber.loc[:, 'y_hat'], df_fiber.loc[:, 'x_hat'], c='Blue', alpha = 0.3, label = 'DNA (unobserved)')

# Draw true positive (true location)
ax2.scatter(df_data.loc[[0, 5000, 15000], 'z_hat'],
            df_data.loc[[0, 5000, 15000], 'y_hat'],
            df_data.loc[[0, 5000, 15000], 'x_hat'],
            alpha = 1, marker='o', s = 300, facecolors='none', edgecolors='r')

ax2.set_title("Output", fontsize = 20);

df_data[df_data['t'] <= 1]

# parms
nm_per_bp = 0.34
pixel_dist = 100
l_p_bp = 150
l_p = l_p_bp*nm_per_bp/pixel_dist
L = 5000

# calculate expected contour length
s = L*nm_per_bp/pixel_dist

# calculate Euclidean distance between observed spots
r_fp = np.linalg.norm(df_data.iloc[0, :3]- df_data.iloc[4, :3])
r_tp = np.linalg.norm(df_data.iloc[0, :3]- df_data.iloc[3, :3])

print("This is the bond probability between spot 1, hyb 0 - spot 1, hyb 1: ", 
      log_bond(l_p, s, r_fp))

print("This is the bond probability between spot 1, hyb 0 - spot 2, hyb 1: ",
      log_bond(l_p, s, r_tp))

cell_pts_input = copy.deepcopy(df_data)
cell_pts_input.rename(columns = {'t':'hyb'}, inplace = True)

mat, path, _ = find_chr(cell_pts_input = cell_pts_input, 
                        gene_dist = [0, 5000, 10000, 15000], 
                        bin_size = 5000,
                        nm_per_bp = 0.34,
                        num_skip = 3,
                        total_num_skip_frac = 1,
                        norm_skip_penalty = False,
                        stretch_factor = 1.01,
                        init_skip_frac = 0,
                        lim_init_skip = False,
                        lim_min_dist=False)

df_data.iloc[path]