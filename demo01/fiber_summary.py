import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pickle
import pandas as pd
from collections import Counter
import numpy as np

"""cell_chr_fiber = pd.DataFrame(columns=["cell", "FOV", "chrnum", "x_hat", "y_hat", "z_hat", "fiber"])

available_chr = glob.glob('seqfishE14_demo01_res_chr*')
for chr_file in available_chr:
    with open(chr_file, 'rb') as f:
        print(chr_file)
        chr_fibers = pd.read_pickle(chr_file)
        for i in range(len(chr_fibers)):
            # for each cell
            all_cell_fibers = chr_fibers[i]
            for c in range(len(all_cell_fibers)):
                cell_fibers = all_cell_fibers[c]
                cell_fibers['chrnum'] = cell_fibers['chr'].str.replace("chr","")
                # append, cell, fov, chrom, first pos
                med_info = cell_fibers.loc[:,['cell','FOV','chrnum','x_hat','y_hat','z_hat']].reset_index().iloc[0]
                # assign fiber number
                med_info['fiber'] = c
                # bin chr into chunks
                cell_chr_fiber = pd.concat([cell_chr_fiber,pd.DataFrame(med_info).T],axis=0)
            # print(cell_chr_fiber)
        print("completed")

cell_chr_fiber.to_csv("cell_chr_fiber_firstpos.csv")"""
cell_chr_fiber = pd.read_csv("cell_chr_fiber_med.csv")
fiber_num = []
cells = list(Counter(cell_chr_fiber['finalcellID']).keys())
chrs = list(Counter(cell_chr_fiber['chrnum']).keys())
print(chrs)
for cell in cells:
    for chr in chrs:
        fibers = cell_chr_fiber.loc[(cell_chr_fiber.finalcellID==cell) & (cell_chr_fiber.chrnum==int(chr))]
        fiber_num.append(len(list(Counter(fibers.fiber).keys())))

"""with open('fiber_nums.pkl', 'wb') as output_file:
    pickle.dump(fiber_num, output_file)"""

"""with open("fiber_nums.pkl", "rb") as input_file:
    fiber_num = pickle.load(input_file)"""

# print(np.mean(fiber_num))
sns.countplot(x=fiber_num)
# plt.yscale('log')
plt.show()