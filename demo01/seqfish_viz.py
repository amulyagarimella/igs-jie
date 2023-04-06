import pickle
import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np

def viz (fibers, chr_data, full_data, chr, fig, nrow, ncol, pos_col, idx, lines=True):
    colors = ['r','b','g','c','m','y']
    ax = fig.add_subplot(nrow,ncol,idx, projection='3d')
    if len(fibers) == 0:
        pass
    for i in range(len(fibers)):
        r = fibers[i].sort_values(by=pos_col,ascending=True)
        ax.scatter(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'], c=colors[i])
        if lines:
            ax.plot(xs=r['x_hat'], ys=r['y_hat'], zs=r['z_hat'], c=colors[i])
    # get data fw no match
    clustered = pd.concat(fibers)
    unclustered = pd.concat([clustered, chr_data], ignore_index=True)
    unclustered = unclustered.drop_duplicates(keep=False)
    # TODO adjust looks
    ax.scatter(xs=unclustered['x_hat'], ys=unclustered['y_hat'], zs=unclustered['z_hat'], c=colors[(len(fibers) + 1) % len(colors)],alpha=.1)
    rest_of_cell = pd.concat([clustered, full_data], ignore_index=True)
    rest_of_cell = rest_of_cell.drop_duplicates(keep=False)
    ax.scatter(xs=rest_of_cell['x_hat'], ys=rest_of_cell['y_hat'], zs=rest_of_cell['z_hat'], c='0.5',alpha=.01)
    ax.set_title(f"{chr}")


## selecting specific chr points (orig data)

# unpickle chr_pts
with open('seqfishE14_demo01_chr_pts.pkl', 'rb') as f:
    chr_pts = pickle.load(f)

chosen_celltype = 'mesc'
fov = 0
cell = 4
# sel_chr_pts = chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['chr'] == chosen_chrom)]['data'].tolist()
sel_pts = pd.concat(chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['FOV'] == fov) & (chr_pts['finalcellID'] == cell)]['data'].tolist())

## selecting fibers
available_chr = glob.glob('seqfishE14_demo01_res_chr*')
available_chr_names = [x.split('.')[0].split('_')[-1] for x in available_chr]
full_res = {}
sel_res = {}
print(available_chr_names)

for c in available_chr:
    with open(c, 'rb') as f:
        full_res[c.split('.')[0].split('_')[-1]] = pickle.load(f)

"""for (key,item) in full_res.items():
    print(item)"""

# get fiber list for each chromosome
for (c,fiber_list) in full_res.items():
    # get all fibers across all cells
    #print(fiber_list)
    #print(len(fiber_list))
    # select list of fibers that matches cell and FOV
    #print(fiber_list[0:5])
    for i in range(len(fiber_list)):
        if len(fiber_list[i]) > 0:
            print(c)
            print(f"{int(list(fiber_list[i][0]['finalcellID'])[0])},{int(list(fiber_list[i][0]['FOV'])[0])}")
            if int(list(fiber_list[i][0]["finalcellID"])[0]) == cell and int(list(fiber_list[i][0]["FOV"])[0]) == fov:
                sel_res[c] = fiber_list[i]
                break
        else:
            break
    

## plot fibers for each chr in sel_res
fig = plt.figure()
res_chr = list(sel_res.keys())
size = int(np.ceil(np.sqrt(len(res_chr))))

for c in res_chr:
    sel_chr_pts = chr_pts[(chr_pts['celltype'] == chosen_celltype) & (chr_pts['chr'] == c) & (chr_pts['FOV'] == fov) & (chr_pts['finalcellID'] == cell)]['data'].tolist()[0]
    #  = int(c.replace("chr",""))
    #print(sel_res[c])
    viz(sel_res[c], sel_chr_pts, sel_pts, c, fig,size,size,'hyb',res_chr.index(c) +1)
plt.tight_layout()
plt.show()
plt.clf()



