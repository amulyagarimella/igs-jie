import pickle
import glob
import pandas as pd
from collections import Counter
import math
import numpy as np
import scipy
from scipy.stats import pearsonr   
import matplotlib.pyplot as plt
import seaborn as sns
import time
np.seterr(all='print')


## compile fiber data across cells and chr
#cell_chr_fiber = pd.DataFrame(columns=['finalcellID', "FOV", "chrnum", "x_hat", "y_hat", "z_hat", "fiber"])
cell_chr_fiber = pd.read_csv("cell_chr_fiber_full.csv")
available_chr = glob.glob('seqfishE14_demo01_res_chr*')
available_chr = []
for chr_file in available_chr:
    print(chr_file)
    chr_fibers = pd.read_pickle(chr_file)
    for i in range(len(chr_fibers)):
        # for each cell
        all_cell_fibers = chr_fibers[i]
        # print(all_cell_fibers)
        for c in range(len(all_cell_fibers)):
            cell_fiber = all_cell_fibers[c].loc[:,['finalcellID','FOV','x_hat','y_hat','z_hat', 'chr']]
            cell_fiber['chrnum'] = cell_fiber['chr'].str.replace("chr","").astype(float)
            cell_fiber = cell_fiber.assign(fiber=str(c))
            # bin chr into chunks
            cell_chr_fiber = pd.concat([cell_chr_fiber, cell_fiber])
        # print(cell_chr_fiber)
    cell_chr_fiber.to_csv("cell_chr_fiber_full.csv")
    print("completed")

cell_chr_fiber = pd.read_csv("cell_chr_fiber_full.csv")

chrs = list(Counter(cell_chr_fiber['chrnum']).keys())
cells = list(Counter(cell_chr_fiber['finalcellID']).keys())

comps = int(len(chrs)*(len(chrs)-1)/2)
chr_vs_chr = np.empty([len(cells),comps])
## create #chr x #chr mean spatial distance matrix for all cells
for c in cells:
    # for chr i
    cvc = np.empty([comps])
    cvc[:] = np.nan
    idx = -1
    for i in range(len(chrs)):
        # for chr j != i
        for j in range(i+1,len(chrs)):
            idx+=1
            if j == i:
                continue
            dists = []
            fibers_i = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(c)) & (cell_chr_fiber.chrnum==int(chrs[i]))].groupby('fiber')
            fibers_j = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(c)) & (cell_chr_fiber.chrnum==int(chrs[j]))].groupby('fiber')
            for k in range(fibers_i.ngroups):
                for l in range(fibers_j.ngroups):
                    distmat = scipy.spatial.distance.cdist(fibers_i.nth(k)[['x_hat','y_hat','z_hat']],fibers_j.nth(l)[['x_hat','y_hat','z_hat']])
                    dists.append(np.nanmedian(distmat))
            cvc[idx] = np.nanmedian(dists)
    print(cvc)
    chr_vs_chr[cells.index(c),:] = cvc
print(chr_vs_chr)
print(chr_vs_chr[~np.isnan(chr_vs_chr).all(axis=1),:])
print(chr_vs_chr[~np.isnan(chr_vs_chr).any(axis=1),:])


cell_v_cell = np.empty([len(cells), len(cells), 2])

same_fov = []
diff_fov = [] 
for i in range(len(cells)):
    for j in range(i+1,len(cells)):
        try:
            x = chr_vs_chr[i]
            y = chr_vs_chr[j]
            nas = np.logical_or(np.isnan(x), np.isnan(y))
            r = pearsonr(x[~nas], y[~nas])[0]
            if len(x[~nas]) < len(chrs) or len(y[~nas]) < len(chrs):
                print("chr missing, skip")
                continue
            #print(f"({i},{j})")
        except:
            # print("not able to calculate correlation")
            continue
        if np.isnan(r):
            continue
        cell_v_cell[i][j][0] = r
        cell_v_cell[j][i][0] = r
        if list(cell_chr_fiber.loc[cell_chr_fiber.finalcellID==cells[i], "FOV"])[0] == list(cell_chr_fiber.loc[cell_chr_fiber.finalcellID==cells[j], "FOV"])[0]:
            same_fov.append(r)
            cell_v_cell[i][j][1] = 1
            cell_v_cell[j][i][1] = 1
        else:
            diff_fov.append(r)
            cell_v_cell[i][j][1] = 0
            cell_v_cell[j][i][1] = 0
"""
# create histograms of fov-based relationships
"""
print(f"# pairs with same fov: {len(same_fov)}")
print(f"# pairs with diff fov: {len(diff_fov)}")
with open('cell_v_cell_corr.pkl', 'wb') as f:
    pickle.dump(cell_v_cell, f)

ax = plt.gca()
sns.kdeplot(x=same_fov,color="r",ax=ax)
sns.kdeplot(x=diff_fov,color="b",ax=ax)
plt.xlim((-1,1))
plt.savefig("corr_memviz.png")
plt.savefig("corr_memviz.svg")
plt.show()

with open('same_fov_corr.pkl', 'wb') as f:
    pickle.dump(same_fov, f)

with open('diff_fov_corr.pkl', 'wb') as f:
    pickle.dump(diff_fov, f)
# plot for fov 3: TODO"""
"""

fov = 3;
fov_table = full_table(full_table.fov == fov,:);
 
figure;
scatter(fov_table.x,fov_table.y,10,colors(fov_table.chromID,:),'filled'); hold on;
axis equal;
 
fov_corr_mat = corr_mat(cell_list(:,1)==fov,cell_list(:,1)==fov);
num_fov_cells = sum(cell_list(:,1)==fov);
%figure; histogram(fov_corr_mat(fov_corr_mat>0))
corr_thresh = 0.4;
 
for i=1:num_fov_cells
    for j=1:num_fov_cells
        if i>j & fov_corr_mat(i,j) > corr_thresh
            i_loc = mean(fov_table{fov_table.cellID == i,3:4},1);
            j_loc = mean(fov_table{fov_table.cellID == j,3:4},1);
            plot([i_loc(1) j_loc(1)],[i_loc(2) j_loc(2)])
        end
    end
end
 

"""


# PCA the corr matrix: TODO