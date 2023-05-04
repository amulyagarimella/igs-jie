import pickle
import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np


colors = ['r','g','g','c','m','y']
## selecting specific chr points (orig data)
available_chr = glob.glob('seqfishE14_demo01_res_chr*')
available_chr_names = [x.split('.')[0].split('_')[-1] for x in available_chr]
full_res = {}
sel_res = {}
with open('seqfishE14_demo01_chr_pts.pkl', 'rb') as f:
    chr_pts = pickle.load(f)
for c in available_chr:
    with open(c, 'rb') as f:
        full_res[c.split('.')[0].split('_')[-1]] = pickle.load(f)


def viz (fibers, full_data, ax, chr, pos_col, color, lines=True):
    if len(fibers) == 0:
        pass
    r = pd.concat(fibers).sort_values(by=pos_col,ascending=True)
    ax.scatter(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'], c=color)
    rest_of_cell = pd.concat([r, full_data], ignore_index=True)
    rest_of_cell = rest_of_cell.drop_duplicates(keep=False)
    print(rest_of_cell.head())
    ax.scatter(xs=rest_of_cell['x_hat'], ys=rest_of_cell['y_hat'], zs=rest_of_cell['z_hat'], c='0.5',alpha=.007)
    """for i in range(len(fibers)):
        r = fibers[i].sort_values(by=pos_col,ascending=True)
        if color is None:
            color = colors[i]
        print(color)
        ax.scatter(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'])
        if lines:
            ax.plot(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'])"""
    
    ax.set_title(f"{chr}")

def plot_cell_chr (cell, axes, chrs_to_plot = None, color = None):
    chosen_celltype = 'mesc'
    # include the other chr points
    ## selecting fibers
    # get fiber list for each chromosome
    if chrs_to_plot is None:
        chrs_to_plot = list(full_res.keys())
    for (c,fiber_list) in full_res.items():
        if c not in chrs_to_plot:
            continue
        for i in range(len(fiber_list)):
            if len(fiber_list[i]) < 1:
                continue
            # print(int(list(fiber_list[i][0]["finalcellID"])[0]))
            # print(list(fiber_list[i][0]["finalcellID"]))
            if int(list(fiber_list[i][0]["finalcellID"])[0]) == cell:
                sel_res[c] = fiber_list[i]
                break

    sel_pts = pd.concat(chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['finalcellID'] == cell)]['data'].tolist())
    for c in chrs_to_plot:
        idx_to_plot = chrs_to_plot.index(c) + 1
        # sel_chr_pts = chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['chr'] == c) & (chr_pts['finalcellID'] == cell)]['data'].tolist()[0]
        #  = int(c.replace("chr",""))
        #print(sel_res[c])
        try:
            viz(sel_res[c], sel_pts, axes[idx_to_plot], c, 'hyb', color=color)
        except:
            continue

fig = plt.figure()
axes = []
for i in range(19):
    axes.append(fig.add_subplot(4,5,i + 1, projection='3d'))

cells_to_plot = [1076, 1120]
for c in cells_to_plot:
    # chrs_to_plot = ['chr15', 'chr12']
    # print(c % 2)
    plot_cell_chr(c, axes, chrs_to_plot=None, color=colors[cells_to_plot.index(c)])
plt.tight_layout()
plt.show()



