import argparse
from scipy.sparse import csr_matrix
from scipy import stats
import anndata
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd

def plot_hists(adata,file_name,clip=10,zeros=False):
    
    data=adata.X.data
    raw_data=adata.layers['counts'].data

    if zeros:
        n_zeros=adata.n_obs*adata.n_vars-len(adata.X.data)
        data = np.append(data,[0]*n_zeros)
        raw_data = np.append(raw_data,[0]*n_zeros)
         
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2,figsize=(10,5))

    ax0.stairs(*np.histogram(data, bins=100), fill=True)
    ax0.set_title("Processed data")
    ax1.stairs(*np.histogram(data[data<clip], bins=100), fill=True)
    ax1.set_title("Clipped processed data")
    ax2.stairs(*np.histogram(raw_data, bins=100), fill=True)
    ax2.set_title("Raw data")
    ax3.stairs(*np.histogram(raw_data[raw_data<clip], bins=100), fill=True)
    ax3.set_title("Clipped raw data")

    for ax in [ax0,ax1,ax2,ax3]:
        ax.set_xlabel("Value")
        ax.set_ylabel("Frequency")
    
    fig.tight_layout()    
    plt.savefig(file_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('train_data',
                        help='path to train data file')
    parser.add_argument('test_data',
                        help='path to test data file')
    args = parser.parse_args()

    adata = sc.read_h5ad(args.train_data)
    adata_test = sc.read_h5ad(args.test_data)

    
    plot_hists(adata,'train_zeros.svg',5,zeros=True)
    plot_hists(adata_test,'test_zeros.svg',5,zeros=True)
    plot_hists(adata_test,'test.svg')
    plot_hists(adata,'train.svg')

    n_donors=len(pd.unique(adata.obs['DonorID']))
    n_cell_types=(len(pd.unique(adata.obs['cell_type'])))
    n_sites=(len(pd.unique(adata.obs['Site'])))

    print(f'Train dataset has {adata.n_obs} observations and {adata_test.n_vars} variables, test dataset has {adata_test.n_obs} observations and {adata_test.n_vars} variables.')
    print(f'No. of patients: {n_donors}, No. of cell types: {n_cell_types}, No. of labs: {n_sites}')


if __name__ == "__main__":
    main()
