import pickle
import glob
import pandas as pd
from collections import Counter
import math
import numpy as np
import scipy
from sklearn.decomposition import PCA
from ppca import PPCA
from scipy.stats import pearsonr   
import matplotlib.pyplot as plt
import seaborn as sns
import time
np.seterr(all='print')

def fiber_data (compression_fn=None):
    cell_chr_fiber = pd.DataFrame(columns=['finalcellID', "FOV", "chrnum", "x_hat", "y_hat", "z_hat", "fiber"])
    available_chr = glob.glob('seqfishE14_demo01_res_chr*')
    #cell_chr_fiber = pd.read_csv("cell_chr_fiber_full.csv")
    #available_chr = [f'seqfishE14_demo01_res_chr{n}.pkl' for n in [1,2,3,10,11,12,13]]
    for chr_file in available_chr:
        print(chr_file)
        chr_fibers = pd.read_pickle(chr_file)
        for i in range(len(chr_fibers)):
            # for each cell
            all_cell_fibers = chr_fibers[i]
            for c in range(len(all_cell_fibers)):
                cell_fiber = all_cell_fibers[c].loc[:,['finalcellID','FOV','x_hat','y_hat','z_hat', 'chr']]
                cell_fiber['chrnum'] = cell_fiber['chr'].str.replace("chr","").astype(float)
                cell_fiber = cell_fiber.assign(fiber=str(c))
                if compression_fn is not None:
                    cell_fiber = cell_fiber.apply(compression_fn,axis=0)
                cell_chr_fiber = pd.concat([cell_chr_fiber, cell_fiber])
        cell_chr_fiber.to_csv("cell_chr_fiber_full.csv")
        print("completed")
    return cell_chr_fiber

cell_chr_fiber = pd.read_csv("cell_chr_fiber_full.csv")

## compile fiber data across cells and chr

def fiber_dist (cell_chr_fiber, cell, chr1, chr2):
    dists = []
    fibers_i = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(cell)) & (cell_chr_fiber.chrnum==int(chr1))].groupby('fiber')
    fibers_j = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==int(cell)) & (cell_chr_fiber.chrnum==int(chr2))].groupby('fiber')
    for k in range(fibers_i.ngroups):
        for l in range(fibers_j.ngroups):
            distmat = scipy.spatial.distance.cdist(fibers_i.nth(k)[['x_hat','y_hat','z_hat']],fibers_j.nth(l)[['x_hat','y_hat','z_hat']])
            dists.append(np.nanmedian(distmat))
    return np.nanmedian(dists)

# TODO change/add diff option
def chr_chr (cell_chr_fiber, comps, cell, chrs):
    cvc = np.empty([comps])
    cvc[:] = np.nan
    idx = -1
    for i in range(len(chrs)):
        # for chr j != i
        for j in range(i+1,len(chrs)):
            idx+=1
            if j == i:
                continue
            cvc[idx] = fiber_dist(cell_chr_fiber, cell, chrs[i], chrs[j])
    return cvc
    
def chr_chr_all (cell_chr_fiber, pkl_name = 'chr_vs_chr'):
    chrs = list(Counter(cell_chr_fiber['chrnum']).keys())
    cells = list(Counter(cell_chr_fiber['finalcellID']).keys())
    comps = int(len(chrs)*(len(chrs)-1)/2)
    chr_vs_chr = np.empty([len(cells),comps])
    ## create #chr x #chr mean spatial distance matrix for all cells
    for c in cells:
        chr_vs_chr[cells.index(c),:] = chr_chr(cell_chr_fiber, comps, c, chrs)
    print(chr_vs_chr)
    print(chr_vs_chr[~np.isnan(chr_vs_chr).all(axis=1),:])
    print(chr_vs_chr[~np.isnan(chr_vs_chr).any(axis=1),:])
    with open(f'{pkl_name}.pkl', 'wb') as f:
        pickle.dump(chr_vs_chr, f)
    return chr_vs_chr

# future TODO ~ if you dropna for PCA, need to save which cell indices ur keeping
def ppca_transform (chr_vs_chr, ncomp=2):
    ppca = PPCA()
    ppca.fit(data=chr_vs_chr,d=ncomp,verbose=True)
    chr_vs_chr_transformed = ppca.transform()
    # TODO make graphs
    return chr_vs_chr_transformed

def cell_cell (cell_chr_fiber, chr_vs_chr, pkl_name = "cell_v_cell_corr", sim = lambda x,y : math.dist(x,y)):
    cells = list(Counter(cell_chr_fiber['finalcellID']).keys())
    cell_vs_cell = np.empty([len(cells), len(cells), 2])
    same_fov = []
    diff_fov = [] 
    for i in range(len(cells)):
        for j in range(i+1,len(cells)):
            try:
                sim = sim(chr_vs_chr[i], chr_vs_chr[j])
                print(sim)
            except:
                print("similarity was not able to be calculated")
                continue
            if np.isnan(sim):
                continue
    
            cell_vs_cell[i][j][0] = sim
            cell_vs_cell[j][i][0] = sim
            if list(cell_chr_fiber.loc[cell_chr_fiber.finalcellID==cells[i], "FOV"])[0] == list(cell_chr_fiber.loc[cell_chr_fiber.finalcellID==cells[j], "FOV"])[0]:
                same_fov.append(sim)
                cell_vs_cell[i][j][1] = 1
                cell_vs_cell[j][i][1] = 1
            else:
                diff_fov.append(sim)
                cell_vs_cell[i][j][1] = 0
                cell_vs_cell[j][i][1] = 0
    print(f"# pairs with same fov: {len(same_fov)}")
    print(f"# pairs with diff fov: {len(diff_fov)}")
    with open(f'{pkl_name}.pkl', 'wb') as f:
        pickle.dump(cell_vs_cell, f)
    with open('same_fov_ppca.pkl', 'wb') as f:
        pickle.dump(same_fov, f)
    with open('diff_fov_ppca.pkl', 'wb') as f:
        pickle.dump(diff_fov, f)
    return cell_vs_cell, same_fov, diff_fov

def fov_hist (same_fov, diff_fov, name):
    ax = plt.gca()
    sns.kdeplot(x=same_fov,color="r",ax=ax)
    sns.kdeplot(x=diff_fov,color="b",ax=ax)
    plt.savefig(f"{name}.png")
    plt.savefig(f"{name}.svg")
    plt.show()

def main ():
    pass

if __name__ == "__main__":
    main()


# todo

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