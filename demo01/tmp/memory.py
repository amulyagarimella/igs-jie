import pickle
import glob
import pandas as pd
from collections import Counter
import math
import numpy as np
from scipy.stats import pearsonr   
import matplotlib.pyplot as plt
import seaborn as sns
import time
np.seterr(all='print')


## compile fiber data across cells and chr

cell_chr_fiber = pd.DataFrame(columns=["cell", "FOV", "chrnum", "x_hat", "y_hat", "z_hat", "fiber"])

available_chr = glob.glob('seqfishE14_demo01_res_chr*')
# available_chr = [f'seqfishE14_demo01_res_chr{n}.pkl' for n in [1,2,3,10,11,12,13]]
for chr_file in available_chr:
    #with open(chr_file, 'rb') as f:
        # ModuleNotFoundError: No module named 'pandas.core.indexes.numeric'
    print(chr_file)
    chr_fibers = pd.read_pickle(chr_file)
    for i in range(len(chr_fibers)):
        # for each cell
        all_cell_fibers = chr_fibers[i]
        # print(all_cell_fibers)
        for c in range(len(all_cell_fibers)):
            cell_fibers = all_cell_fibers[c].loc[:,['finalcellID','FOV','hyb','x_hat','y_hat','z_hat', 'chr']]
            cell_fibers['chrnum'] = cell_fibers['chr'].str.replace("chr","").astype(float)
            # append, cell, fov, chrom, median pos
            #if (i >= 28):
                #print(cell_fibers)
                #print(cell_fibers.median(numeric_only=True))
            try:
                med_info = cell_fibers.sort_values(by='hyb').median(numeric_only=True)
            except:
                #print(cell_fibers)
                raise
            # assign fiber number
            # do chr num a better way TODO
            med_info['fiber'] = c
            # bin chr into chunks
            cell_chr_fiber = pd.concat([cell_chr_fiber,pd.DataFrame(med_info).T],axis=0)
        # print(cell_chr_fiber)
    cell_chr_fiber.to_csv("cell_chr_fiber_med.csv")
    print("completed")

cell_chr_fiber = pd.read_csv("cell_chr_fiber_med.csv")

## create #chr x #chr mean spatial distance matrix for all cells
cells = list(Counter(cell_chr_fiber['finalcellID']).keys())
chrs = list(Counter(cell_chr_fiber['chrnum']).keys())

print(cells)
print(chrs)
# initialize matrix
chr_vs_chr = []

#todo use set insted of counter

# for each cell
for c in cells:
    # for chr i
    cvc = np.empty([len(chrs), len(chrs)])
    cvc[:] = np.nan
    for i in range(len(chrs)):
        # for chr j != i
        for j in range(i+1,len(chrs)):
            if j == i:
                continue
            dists = []
            fibers_i = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(c)) & (cell_chr_fiber.chrnum==int(chrs[i]))]
            f_i_pts = fibers_i.loc[fibers_i.fiber==int(i),['x_hat', 'y_hat','z_hat']].values
            #print(f_i_pts)
            fibers_j = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(c)) & (cell_chr_fiber.chrnum==int(chrs[j]))]
            f_j_pts = fibers_j.loc[fibers_j.fiber==int(j),['x_hat', 'y_hat','z_hat']].values
            #print(f_j_pts)
            if len(f_i_pts) == 0 or len(f_j_pts) == 0:
                continue
            for f_i_pt in f_i_pts:
                for f_j_pt in f_j_pts:
                    d = math.dist(f_i_pt,f_j_pt)
                    dists.append(d)
                    #print(d)
            # take med
            # or, 2 fibers close => merge
            #print(dists)
            cvc[i][j] = np.nanmean(dists)
    chr_vs_chr.append(cvc)
    #print(cvc)
    #print(chr_vs_chr)

# create correlation matrix; # cells x # cells

cell_v_cell = np.empty([len(cells), len(cells), 2])

same_fov = []
diff_fov = [] 
for i in range(len(cells)):
    for j in range(i+1,len(cells)):
        try:
            x = chr_vs_chr[i].flatten()
            y = chr_vs_chr[j].flatten()
            nas = np.logical_or(np.isnan(x), np.isnan(y))
            r = pearsonr(x[~nas], y[~nas])[0]
            if len(x[~nas]) < len(chrs) or len(y[~nas]) < 3:
                print("chr missing, skip")
                continue
            #print(f"({i},{j})")
            if (r > .9):
                print(f"({i},{j})")
                print(x[~nas])
                print(y[~nas])
        except:
            # print("not able to calculate correlation")
            continue
        if np.isnan(r):
            continue
        print(r)
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

# create histograms of fov-based relationships

print(f"# pairs with same fov: {len(same_fov)}")
print(f"# pairs with diff fov: {len(diff_fov)}")


ax = plt.gca()
sns.kdeplot(x=same_fov,color="r",ax=ax)
sns.kdeplot(x=diff_fov,color="b",ax=ax)
plt.xlim((-1,1))
plt.savefig("allchr_memviz.png")
plt.savefig("allchr_memviz.svg")
plt.show()

with open('same_fov_corr.pkl', 'wb') as f:
    pickle.dump(same_fov, f)

with open('diff_fov_corr.pkl', 'wb') as f:
    pickle.dump(diff_fov, f)
# plot for fov 3: TODO
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