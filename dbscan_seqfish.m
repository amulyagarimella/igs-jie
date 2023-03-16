%% read table
 
full_table = readtable('processed/seqfish/coord_tables/E14_rep1_coord_table.txt');
 
% convert pixel coordinates to microns
full_table.x = full_table.x*0.103;
full_table.y = full_table.y*0.103;
full_table.z = full_table.z*0.250;
 
%% get list of cells
 
cell_list = unique(full_table{:,1:2},'rows');
num_cells = size(cell_list,1);
num_chrs = 20;
 
%% visualize single cell
 
colors = distinguishable_colors(101);
 
idx = 4;
fov = cell_list(idx,1); 
cell = cell_list(idx,2); 
 
cell_table = full_table(full_table.fov == fov & full_table.cellID == cell,:);
 
figure;
scatter3(cell_table.x,cell_table.y,cell_table.z,10,colors(cell_table.chromID,:),'filled'); view(-45,30);
axis equal;
 
%% visualize clusters for single cell (colored by either chromosome, position, or dbscan cluster)
 
figure;
p = tight_subplot(3,7,[.01 .01],[.01 .01],[.01 .01]);
 
axes(p(1))
scatter3(cell_table.x,cell_table.y,cell_table.z,10,colors(cell_table.chromID,:),'filled')
view(-45,30); axis equal;
title(sprintf('fov %d, cell %d: %d reads',fov,cell,size(cell_table,1)))
 
for i=1:num_chrs
axes(p(i+1));
 
chr_table = cell_table(cell_table.chromID == i,:);
chr_table.cluster = dbscan(chr_table{:,3:5},0.75,5);
 
scatter3(cell_table.x,cell_table.y,cell_table.z,10,[.7 .7 .7],'filled','MarkerFaceAlpha',0.25); hold on;
%scatter3(chr_table.x,chr_table.y,chr_table.z,10,colors(chr_table.chromID,:),'filled');
%scatter3(chr_table.x,chr_table.y,chr_table.z,10,chr_table.Start,'filled');
scatter3(chr_table.x,chr_table.y,chr_table.z,10,cluster,'filled'); 
 
view(-45,30); axis equal;
title(sprintf('chr%d: %d reads',i,size(chr_table,1)))
 
end
 
%% dbscan clustering for all cells
 
full_table.cluster = zeros(size(full_table,1),1);
 
for cell_idx=1:num_cells; disp(cell_idx)
 
fov = cell_list(cell_idx,1);
cell = cell_list(cell_idx,2);
 
cell_table = full_table(full_table.fov == fov & full_table.cellID == cell,:);
cell_table.cluster = zeros(size(cell_table,1),1);
 
for chr=1:num_chrs
    chr_table = cell_table(cell_table.chromID == chr,:);
    if size(chr_table,1) <= 0
        continue
    end
    
    init_cluster = dbscan(chr_table{:,3:5},0.75,5);
    
    % reorder clusters by 
    
    [C ia ic] = unique(init_cluster(init_cluster>0));
    counts = accumarray(ic,1);
    [tmp sort_idx] = sort(counts,'descend');
    
    cluster = zeros(size(init_cluster));
    for i=1:max(init_cluster)
       cluster(init_cluster==sort_idx(i)) = i;
    end
    
    cell_table.cluster(cell_table.chromID == chr) = cluster;
    
end
 
full_table.cluster(full_table.fov == fov & full_table.cellID == cell) = cell_table.cluster;
 
end
 
%% create #chr x #chr mean spatial distance matrix for all cells
 
chr_dist_mat = zeros(num_cells,num_chrs,num_chrs);
 
for cell_idx=1:num_cells; disp(cell_idx)
 
fov = cell_list(cell_idx,1);
cell = cell_list(cell_idx,2);
    
cell_table = full_table(full_table.fov == fov & full_table.cellID == cell,:);
 
for chr_idx_1=1:num_chrs
    for chr_idx_2=1:num_chrs
        
        clust_mat = zeros(2,2);
        
        for clust_idx_1=1:2 
            for clust_idx_2=1:2
                mean_loc_1 = nanmean(cell_table{cell_table.chromID == chr_idx_1 & (cell_table.cluster == clust_idx_1),3:5},1);
                mean_loc_2 = nanmean(cell_table{cell_table.chromID == chr_idx_2 & (cell_table.cluster == clust_idx_2),3:5},1);
                clust_mat(clust_idx_1,clust_idx_2) = pdist2(mean_loc_1,mean_loc_2);
            end
        end
        
        %chr_mat(chr1,chr2) = mean(clust_mat,'all'); 
        chr_dist_mat(cell_idx,chr_idx_1,chr_idx_2) = mean(clust_mat(clust_mat>0),'all'); % removes zeros
        
    end
end
 
end
 
chr_dist_mat(isnan(chr_dist_mat)) = 0;
 
%%
 
corr_mat = zeros(num_cells,num_cells);
 
for i=1:num_cells
    for j=1:num_cells
        if i>j
            corr_mat(i,j) = corr(reshape(chr_dist_mat(i,:,:),[],1),reshape(chr_dist_mat(j,:,:),[],1));
        end
    end
end
 
%%
 
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
 
%%
 
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
 
%% save table with clusters
 
writetable(full_table,'processed/seqfish/coord_tables/E14_rep1_coord_table__with_cluster.txt');