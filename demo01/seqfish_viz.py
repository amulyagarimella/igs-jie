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


def viz (fibers, chr_data, full_data, chr, pos_col, ax, color, lines=True):
    if len(fibers) == 0:
        pass
    r = pd.concat(fibers).sort_values(by=pos_col,ascending=True)
    ax.scatter(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'], c=color)
    # if lines:
    #     ax.plot(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'])
    """for i in range(len(fibers)):
        r = fibers[i].sort_values(by=pos_col,ascending=True)
        if color is None:
            color = colors[i]
        print(color)
        ax.scatter(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'])
        if lines:
            ax.plot(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'])"""
    # get data fw no match
    # clustered = pd.concat(fibers)
    # unclustered = pd.concat([clustered, chr_data], ignore_index=True)
    # unclustered = unclustered.drop_duplicates(keep=False)
    # TODO adjust looks
    #if color is None:
        #color = colors[(len(fibers) + 1) % len(colors)]
    #ax.scatter(xs=unclustered['x_hat'], ys=unclustered['y_hat'], zs=unclustered['z_hat'], c=color,alpha=.1)
    # rest_of_cell = pd.concat([clustered, full_data], ignore_index=True)
    # rest_of_cell = rest_of_cell.drop_duplicates(keep=False)
    #ax.scatter(xs=rest_of_cell['x_hat'], ys=rest_of_cell['y_hat'], zs=rest_of_cell['z_hat'], c='0.5',alpha=.01)
    ax.set_title(f"{chr}")

def plot_cell_chr (cell, axes, chrs_to_plot = None, color = None):
    chosen_celltype = 'mesc'
    # fov = 0
    # sel_chr_pts = chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['chr'] == chosen_chrom)]['data'].tolist()
    sel_pts = pd.concat(chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['finalcellID'] == cell)]['data'].tolist())
    print(set(list(sel_pts.chr)))
    ## selecting fibers
    # get fiber list for each chromosome
    for (c,fiber_list) in full_res.items():
        if chrs_to_plot is not None:
            if c not in chrs_to_plot:
                continue
        for i in range(len(fiber_list)):
            if len(fiber_list[i]) > 0:
                print(c)
                # print(f"{int(list(fiber_list[i][0]['finalcellID'])[0])},{int(list(fiber_list[i][0]['FOV'])[0])}")
                if int(list(fiber_list[i][0]["finalcellID"])[0]) == cell:
                    sel_res[c] = fiber_list[i]
                    break
            else:
                break
    ## plot fibers for each chr in sel_res
    
    res_chr = list(sel_res.keys())
    if chrs_to_plot is None:
        size = int(np.ceil(np.sqrt(19)))
        chrs_to_plot = res_chr
    else:
        size = int(np.ceil(np.sqrt(len(chrs_to_plot))))

    for c in chrs_to_plot:
        print(c)
        if chrs_to_plot is None:
            try:
                idx_to_plot = int(c.replace("chr",""))
            except:
                continue
        else:
            idx_to_plot = chrs_to_plot.index(c) + 1
        sel_chr_pts = chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['chr'] == c) & (chr_pts['finalcellID'] == cell)]['data'].tolist()[0]
        #  = int(c.replace("chr",""))
        #print(sel_res[c])
        print(color)
        viz(sel_res[c], sel_chr_pts, sel_pts, c, 'hyb', axes[chrs_to_plot.index(c)], color=color)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1, projection='3d')
ax2 = fig.add_subplot(1,2,2, projection='3d')

for c in [106]:
    chrs_to_plot = ['chr7', 'chr12']
    print(c % 2)
    plot_cell_chr(c,(ax1, ax2), chrs_to_plot=chrs_to_plot, color=colors[c % 2])
plt.tight_layout()
plt.show()



