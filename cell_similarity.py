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

## compile fiber data across cells and chr

def allpairs_median_dists(fibers_i, fibers_j, x_colname,y_colname,z_colname):
    dists = []
    for k in range(fibers_i.ngroups):
        for l in range(fibers_j.ngroups):
            distmat = scipy.spatial.distance.cdist(fibers_i.nth(k)[[x_colname, y_colname, z_colname]],fibers_j.nth(l)[[x_colname, y_colname, z_colname]])
            dists.append(np.nanmedian(distmat))
    return dists

def first_second_points_dists(fibers_i, fibers_j, x_colname,y_colname,z_colname):
    dists = []
    for k in range(fibers_i.ngroups):
        for l in range(fibers_j.ngroups):
            i_pts = fibers_i.nth(k)[[x_colname, y_colname, z_colname]]
            j_pts = fibers_j.nth(l)[[x_colname, y_colname, z_colname]]
            ai = i_pts[0]
            aj = j_pts[0]
            bi = i_pts[-1]
            bj = j_pts[-1]
            dists.append(np.min(abs(ai-bj), abs(aj-bi)))
    return dists

def fiber_dist (cell_chr_fiber, cell, chr1, chr2,fiber_colname,x_colname,y_colname,z_colname, cellid_colname,chrom_colname, dists_fn = allpairs_median_dists, final_dist_fn = np.nanmedian):
    fibers_i = cell_chr_fiber.loc[(cell_chr_fiber[cellid_colname]==int(cell)) & (cell_chr_fiber[chrom_colname]==int(chr1))].groupby(fiber_colname)
    fibers_j = cell_chr_fiber.loc[(cell_chr_fiber[cellid_colname]==int(cell)) & (cell_chr_fiber[chrom_colname]==int(chr2))].groupby(fiber_colname)
    dists = dists_fn(fibers_i, fibers_j, x_colname,y_colname,z_colname)
    return final_dist_fn(dists)


def chr_chr (cell_chr_fiber, comps, cell, chrs,fiber_colname,x_colname,y_colname,z_colname, cellid_colname,chrom_colname):
    cvc = np.empty([comps])
    cvc[:] = np.nan
    idx = -1
    for i in range(len(chrs)):
        # for chr j != i
        for j in range(i+1,len(chrs)):
            idx+=1
            if j == i:
                continue
            cvc[idx] = fiber_dist(cell_chr_fiber, cell, chrs[i], chrs[j],fiber_colname,x_colname,y_colname,z_colname,cellid_colname,chrom_colname)
    return cvc
    
def chr_chr_all (cell_chr_fiber, pkl_name = 'chr_vs_chr_allpairs',fiber_colname='fiber',x_colname='x_hat',y_colname='y_hat',z_colname='z_hat',cellid_colname='finalcellID',chrom_colname='chrnum'):
    chrs = list(Counter(cell_chr_fiber[chrom_colname]).keys())
    cells = list(Counter(cell_chr_fiber[cellid_colname]).keys())
    comps = int(len(chrs)*(len(chrs)-1)/2)
    chr_vs_chr = np.empty([len(cells),comps])
    ## create #chr x #chr mean spatial distance matrix for all cells
    for c in cells:
        chr_vs_chr[cells.index(c),:] = chr_chr(cell_chr_fiber, comps, c, chrs,fiber_colname=fiber_colname,x_colname=x_colname,y_colname=y_colname,z_colname=z_colname,cellid_colname=cellid_colname,chrom_colname=chrom_colname)
    print(chr_vs_chr)
    print(chr_vs_chr[~np.isnan(chr_vs_chr).all(axis=1),:])
    print(chr_vs_chr[~np.isnan(chr_vs_chr).any(axis=1),:])
    with open(f'{pkl_name}.pkl', 'wb') as f:
        pickle.dump(chr_vs_chr, f)
    return chr_vs_chr

# input both as np arrays
def pcs_vs_orig (pcs, orig):
    # TODO orig var names
    for i in range(len(pcs)):
        for j in range(len(orig)):
            np.pearsonr(pcs[i],orig[j])

# future TODO ~ if you dropna for PCA, need to save which cell indices ur keeping
def ppca_transform (chr_vs_chr, ncomp=10,var_exp_name="var_exp_allpairs_ppca", pkl_name="chr_v_chr_allpairs_ppca_transformed"):
    ppca = PPCA()
    ppca.fit(data=chr_vs_chr,d=ncomp,verbose=True)
    chr_vs_chr_transformed = ppca.transform()
    variance_explained = list(ppca.var_exp)
    model_params = ppca.C
    cumulative_var_exp = [variance_explained[i]*100 for i in range(1,len(variance_explained))]
    sns.lineplot(x=range(1,len(cumulative_var_exp)+1), y=cumulative_var_exp)
    #xticks
    plt.xticks(range(1,len(cumulative_var_exp)+1))
    with open(f'{pkl_name}_{ncomp}_var_exp.pkl', 'wb') as f:
        pickle.dump(variance_explained, f)
    #with open(f'{pkl_name}_{ncomp}_components.pkl', 'wb') as f:
    #    pickle.dump(components, f)
    with open(f'{pkl_name}_{ncomp}_params.pkl', 'wb') as f:
        pickle.dump(model_params, f)
    plt.savefig(f"{var_exp_name}_{ncomp}.svg")
    with open(f'{pkl_name}_{ncomp}.pkl', 'wb') as f:
        pickle.dump(chr_vs_chr_transformed, f)
    plt.show()
    plt.clf()
    return chr_vs_chr_transformed

def cell_cell (cell_chr_fiber, chr_vs_chr, pkl_name = "cell_v_cell_allpairs_ppca", sim_fn = math.dist, cellid_colname='finalcellID',fov_colname='fov'):
    cells = set(cell_chr_fiber[cellid_colname])
    cell_vs_cell = np.empty([len(cells), len(cells), 2])
    same_fov = []
    diff_fov = [] 
    for i in range(len(cells)):
        for j in range(i+1,len(cells)):
            try:
                sim = sim_fn(chr_vs_chr[i], chr_vs_chr[j])
            except:
                print("similarity was not able to be calculated")
                continue
            if np.isnan(sim):
                continue
            cell_vs_cell[i][j][0] = sim
            cell_vs_cell[j][i][0] = sim
            if list(cell_chr_fiber.loc[cell_chr_fiber[cellid_colname]==cells[i], fov_colname])[0] == list(cell_chr_fiber.loc[cell_chr_fiber[cellid_colname]==cells[j], fov_colname])[0]:
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
    with open(f'{pkl_name}_samefov.pkl', 'wb') as f:
        pickle.dump(same_fov, f)
    with open(f'{pkl_name}_difffov.pkl', 'wb') as f:
        pickle.dump(diff_fov, f)
    return cell_vs_cell, same_fov, diff_fov

def fov_hist (same_fov, diff_fov, name):
    ax = plt.gca()
    sns.kdeplot(x=same_fov,color="r",ax=ax,label='same FOV')
    sns.kdeplot(x=diff_fov,color="b",ax=ax,label='diff FOV')
    #plt.xlim((-,20))
    plt.title('Distance between cells')
    plt.legend()
    plt.savefig(f"{name}.png")
    plt.savefig(f"{name}.svg")
    plt.show()

def main ():
    cell_chr_fiber = pd.read_csv("clonal/E14_clonal_24hr_coord_table_with_cluster.csv")
    cell_chr_fiber.sort_values(by="Start", inplace=True)
    cell_chr_fiber['newCellID'] =cell_chr_fiber["fov"].astype(str) + cell_chr_fiber["cellID"].apply(lambda a : str(a).zfill(2))
    cell_chr_fiber['newCellID'] = cell_chr_fiber['newCellID'].astype(int)
    # chr_vs_chr = chr_chr_all(cell_chr_fiber,pkl_name='clonal/chr_chr_allpairs_newID',fiber_colname='cluster',x_colname='x',y_colname='y',z_colname='z',cellid_colname='newCellID',chrom_colname='chromID')
    ncomp=50
    #chr_vs_chr_transformed = ppca_transform(chr_vs_chr,ncomp=ncomp,var_exp_name='clonal/var_exp_allpairs_ppca_newID',pkl_name='clonal/chr_chr_allpairs_ppca_newID')
    with open('clonal/chr_chr_allpairs_ppca_newID_100.pkl', 'rb') as f:
        chr_vs_chr_transformed = pickle.load(f)
    cell_vs_cell, same_fov, diff_fov = cell_cell(cell_chr_fiber, chr_vs_chr_transformed, pkl_name = f"clonal/cell_cell_allpairs_ppca_newID_{ncomp}", cellid_colname='newCellID')
    fov_hist([abs(x) for x in same_fov], [abs(x) for x in diff_fov], name=f"clonal/fov_comparison_ppca_newID_{ncomp}")

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