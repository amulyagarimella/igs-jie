import pickle
import glob
import pandas as pd
from collections import Counter
import math
import numpy as np
from scipy.stats.stats import pearsonr   
import matplotlib.pyplot as plt
import seaborn as sns

## compile fiber data across cells and chr

cell_chr_fiber = pd.DataFrame(columns=["cell", "FOV", "chrnum", "x_hat", "y_hat", "z_hat", "fiber"])

available_chr = glob.glob('seqfishE14_demo01_res_chr*')
for chr_file in available_chr[0:2]:
    with open(chr_file, 'rb') as f:
        chr_fibers = pickle.load(f)
        for i in range(len(chr_fibers)):
            # for each cell
            all_cell_fibers = chr_fibers[i]
            for c in range(len(all_cell_fibers)):
                cell_fibers = all_cell_fibers[c]
                cell_fibers['chrnum'] = cell_fibers['chr'].str.replace("chr","").astype(int)
                # append, cell, fov, chrom, median pos
                med_info = cell_fibers.loc[:,['cell','FOV','chrnum','x_hat','y_hat','z_hat']].median(numeric_only=None)
                # assign fiber number
                med_info['fiber'] = c
                cell_chr_fiber = pd.concat([cell_chr_fiber,pd.DataFrame(med_info).T],axis=0)
            # print(cell_chr_fiber)

## create #chr x #chr mean spatial distance matrix for all cells

cells = list(Counter(cell_chr_fiber['cell']).keys())
chrs = list(Counter(cell_chr_fiber['chrnum']).keys())

# initialize matrix
chr_vs_chr = []

# for each cell
for c in cells:
    # for chr i
    cvc = np.empty([len(chrs), len(chrs)])
    for i in range(len(chrs)):
        # for chr j != i
        for j in range(len(chrs)):
            if j != i:
                dists = []
                fibers_i = cell_chr_fiber.loc[(cell_chr_fiber.cell==int(c)) & (cell_chr_fiber.chrnum==int(chrs[i]))]
                f_i_pts = fibers_i.loc[fibers_i.fiber==int(i),['x_hat', 'y_hat','z_hat']].values
                fibers_j = cell_chr_fiber.loc[(cell_chr_fiber.cell==int(c)) & (cell_chr_fiber.chrnum==int(chrs[j]))]
                f_j_pts = fibers_j.loc[fibers_j.fiber==int(j),['x_hat', 'y_hat','z_hat']].values
                for f_i_pt in f_i_pts:
                    for f_j_pt in f_j_pts:
                        d = math.dist(f_i_pt,f_j_pt)
                        dists.append(d)
                cvc[i][j] = np.mean(dists)
            else:
                pass
    chr_vs_chr.append(cvc)

# create correlation matrix; # cells x # cells

cell_v_cell = np.empty([len(cells), len(cells), 2])

same_fov = []
diff_fov = [] 
for i in range(len(cells)):
    for j in range(i,len(cells)):
        try:
            r = pearsonr(chr_vs_chr[i].flatten(), chr_vs_chr[j].flatten())[0]
        except:
            pass
        cell_v_cell[i][j][0] = r
        cell_v_cell[j][i][0] = r
        if list(cell_chr_fiber.loc[cell_chr_fiber.cell==cells[i], "FOV"])[0] == list(cell_chr_fiber.loc[cell_chr_fiber.cell==cells[j], "FOV"])[0]:
            same_fov.append(r)
            cell_v_cell[i][j][1] = 1
            cell_v_cell[j][i][1] = 1
        else:
            diff_fov.append(r)
            cell_v_cell[i][j][1] = 0
            cell_v_cell[j][i][1] = 0

# create histograms of fov-based relationships: TODO

print(len(same_fov))

sns.displot(x=same_fov,color="r")
sns.displot(x=diff_fov,color="b")
plt.xlim((-1,1))
plt.savefig("chr1_chr2_memviz.png")
plt.savefig("chr1_chr2_memviz.svg")


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