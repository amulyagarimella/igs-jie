import pickle

with open('seqfishE14_demo01_res_chrom_all.pkl', 'rb') as f:
    chr_fibers = pickle.load(f)


## assign fiber numbers for all cells

for 

## create #chr x #chr mean spatial distance matrix for all cells

# find mean 3-pt location for each fiber for each chr; store this

# for chr i
    # for chr j ≠ i
        # for each fiber of chr i
            # for all fibers of chr j
                # determine dist between two fibers; add to list
        # mean of list: distance between chr i <-> j






"""
corr_mat = zeros(num_cells,num_cells);
 
for i=1:num_cells
    for j=1:num_cells
        if i>j
            corr_mat(i,j) = corr(reshape(chr_dist_mat(i,:,:),[],1),reshape(chr_dist_mat(j,:,:),[],1));
        end
    end
end
 

"""

# create correlation matrix; # cells x # cells

"""
% creates a matrix where (i,j) = 1 indicates that cell i and j were from the same fov
same_fov = double(repmat(cell_list(:,1),1,num_cells) == repmat(cell_list(:,1),1,num_cells)')
 
same_fov_corr = corr_mat(find(same_fov));
same_fov_corr = same_fov_corr(same_fov_corr~=0);
 
diff_fov_corr = corr_mat(find(~same_fov));
diff_fov_corr = diff_fov_corr(diff_fov_corr~=0);
 
figure;
ksdensity(same_fov_corr); hold on;
ksdensity(diff_fov_corr); hold on;
 
figure;
histogram(same_fov_corr,-1:0.1:1); hold on;
histogram(diff_fov_corr,-1:0.1:1); hold on;
xlim([-1 1])

"""

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


# PCA the corr matrix 