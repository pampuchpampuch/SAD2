import argparse
from scipy.sparse import csr_matrix
from scipy import stats
import anndata
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter

def plot_hists(adata,file_name,clip=False):

    X_shape = adata.shape
    n_zeros=X_shape[0]*X_shape[1]-len(adata.X.data)
    data_non_zero=adata.X.data
    raw_data_non_zero=adata.layers['counts'].data
    if clip:
        data_non_zero=data_non_zero[data_non_zero<clip[0]]
        raw_data_non_zero=raw_data_non_zero[raw_data_non_zero<clip[1]]

    data_with_zero=np.append(data_non_zero,[0]*n_zeros)  
    raw_data=np.append(raw_data_non_zero,[0]*n_zeros)
    

    fig, ((ax0, ax1),(ax2, ax3)) = plt.subplots(2,2,figsize=(10,5))

    # '''
    # name axes and titles
    # '''
    # '''
    # do the same for test and train?
    # '''
    ax0.hist(data_with_zero, bins=100)
    ax0.set_title("Processed data with zeros")
    ax1.hist(data_non_zero, bins=100)
    ax1.set_title("Processed data without zeros")
    ax2.hist(raw_data, bins=100)
    ax2.set_title("Raw data with zeros")
    ax3.hist(raw_data_non_zero, bins=100)
    ax3.set_title("Raw data without zeros")
    fig.tight_layout() 

    for ax in [ax0,ax1,ax2,ax3]:
        ax.set_xlabel("Read count")
        ax.set_ylabel("No. of samples")



    # t_data=anndata.AnnData(np.reshape(raw_data.copy(),adata.shape))
    # # print(t_data)

    # n = sc.pp.normalize_total(t_data, target_sum=1e4, inplace=False)['X']
    # l = sc.pp.log1p(t_data, copy=True)
    # # nl = sc.pp.log1p(n, copy=True)
    # s = sc.pp.scale(t_data, copy=True)
    # # ls = sc.pp.scale(l, copy=True)
    # # ns=sc.pp.scale(n, copy=True)
    # # nls = sc.pp.scale(nl, copy=True)
    # # print(t_data)
    # # print(adata.obs)

    # # v=t_data.flatten()

    # ax4.hist(n, bins=100)
    # ax5.hist(l, bins=100)
    # ax5.hist(s, bins=100)
    # # ax5.hist(ls, bins=100)
    # # ax5.hist(ns, bins=100)
    # # ax5.hist(ls, bins=100)
    # # ax5.hist(nls, bins=100)
    plt.savefig(file_name)
    plt.show()
    
    # plot_hists(adata)




parser = argparse.ArgumentParser()
parser.add_argument('data_file',
                    help='path to data file')
args = parser.parse_args()

adata = sc.read_h5ad(args.data_file)

# plot_hists(adata,"fig.svg",(10,10))

# adata = sc.read_h5ad(train_file)
n_zeros=adata.n_obs*adata.n_vars-len(adata.X.data)
raw_data_non_zero=adata.layers['counts'].data
raw_data=np.append(raw_data_non_zero,[0]*n_zeros)

d=anndata.AnnData(np.reshape(raw_data.copy(),adata.shape))


d = sc.pp.normalize_total(d, target_sum=1e4, inplace=False)['X']
# print('n)')
# d = sc.pp.log1p(d, copy=True)
# print('l')

d = sc.pp.scale(d, copy=True)
print('s')

# print(t_data)
# print(adata.obs)

# d=d.X
d=d.flatten()

plt.stairs(*np.histogram(d, bins=100),fill=True)
plt.savefig('ns')