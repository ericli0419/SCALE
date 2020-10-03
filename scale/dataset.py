#!/usr/bin/env python
"""
# Author: Yuzhe Li
# Created Time : Tue Sep 15 19:15:31 CST 2020
# File Name: dataset.py
# Description:
"""

import os
import numpy as np
import pandas as pd
import scipy
import time
from tqdm import tqdm

from torch.utils.data import Dataset
import anndata as ad
import scanpy as sc
from sklearn.preprocessing import maxabs_scale, MaxAbsScaler
from multiprocessing import Pool, cpu_count

from glob import glob

np.warnings.filterwarnings('ignore')



class SingleCellDataset(Dataset):
    """
    Dataset for dataloader
    """

    def __init__(self, adata,transforms=[]):
        self.adata = adata
        self.shape = adata.shape
        self.n_cells, self.n_peaks = adata.shape[0],adata.shape[1]
        self.data, self.peaks, self.barcode = adata.X, adata.var.index, adata.obs.index
        self.obs=adata.obs
        self.obsm=adata.obsm
        self.var=adata.var
        self.layers=adata.layers
        for transform in transforms:
            adata.X = transform(adata.X)
    
    def write(self,path,compression='gzip'):
        self.adata.write(path,compression='gzip')


    def __len__(self):
        return self.adata.X.shape[0]

    '''
    def __getitem__(self, index):
        x = self.adata.X[index].toarray().squeeze()
        return x, index
    '''
    
    def __getitem__(self, index):
        data = self.adata.X[index];
        if type(data) is not np.ndarray:
            data = data.toarray().squeeze()
        return data

    def info(self):
        print("\n===========================")
        print("Dataset Info")
        print('Cell number: {}\nPeak number: {}'.format(self.n_cells, self.n_peaks))
        print('===========================\n')


def filter_cells(
        adata,
        min_counts=None,
        min_peaks=0,
        max_counts=None,
        max_peaks=None,
        inplace=True,
        copy=False):
    """Filter cell outliers based on counts and numbers of genes expressed.

    For instance, only keep cells with at least `min_counts` counts or
    `min_features` genes expressed. This is to filter measurement outliers,
    i.e. “unreliable” observations.

    Only provide one of the optional parameters ``min_counts``, ``min_features``,
    ``max_counts``, ``max_features`` per call.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape ``n_obs`` × ``n_vars``.
        Rows correspond to cells and columns to genes.
    min_counts
        Minimum number of counts required for a cell to pass filtering.
    min_features
        Minimum number of genes expressed required for a cell to pass filtering.
    max_counts
        Maximum number of counts required for a cell to pass filtering.
    max_features
        Maximum number of genes expressed required for a cell to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on ``inplace``, returns the following arrays or directly subsets
    and annotates the data matrix:

    cells_subset : numpy.ndarray
        Boolean index mask that does filtering. ``True`` means that the
        cell is kept. ``False`` means the cell is removed.
    number_per_cell : numpy.ndarray
        Depending on what was tresholded (``counts`` or ``genes``), the array stores
        ``n_counts`` or ``n_cells`` per gene.
    """
    if copy:
        adata_copy = sc.pp.filter_cells(adata, min_counts, min_genes=min_peaks, max_counts=max_counts,
                                        max_genes=max_peaks, inplace=inplace, copy=copy)
        adata_copy.obs['nb_features'] = adata_copy.obs['n_genes']
        del adata_copy.obs['n_genes']
        return (adata_copy)
    else:
        sc.pp.filter_cells(adata, min_counts, min_genes=min_peaks, max_counts=max_counts,
                           max_genes=max_peaks, inplace=inplace, copy=copy)
        adata.obs['nb_features'] = adata.obs['n_genes']
        del adata.obs['n_genes']





def filter_features(adata,
                    min_counts = None,
                    min_cells = None,
                    max_counts = None,
                    max_cells = None,
                    inplace=True,
                    copy=False):
    """Filter features based on number of cells or counts.

    Keep features that have at least ``min_counts`` counts or are expressed in at
    least ``min_cells`` cells or have at most ``max_counts`` counts or are expressed
    in at most ``max_cells`` cells.

    Only provide one of the optional parameters ``min_counts``, ``min_cells``,
    ``max_counts``, ``max_cells`` per call.

    Parameters
    ----------
    data
        An annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    min_counts
        Minimum number of counts required for a gene to pass filtering.
    min_cells
        Minimum number of cells expressed required for a gene to pass filtering.
    max_counts
        Maximum number of counts required for a gene to pass filtering.
    max_cells
        Maximum number of cells expressed required for a gene to pass filtering.
    inplace
        Perform computation inplace or return result.

    Returns
    -------
    Depending on `inplace`, returns the following arrays or directly subsets
    and annotates the data matrix

    gene_subset : numpy.ndarray
        Boolean index mask that does filtering. `True` means that the
        gene is kept. `False` means the gene is removed.
    number_per_gene : numpy.ndarray
        Depending on what was tresholded (`counts` or `cells`), the array stores
        `n_counts` or `n_cells` per gene.
    """
    if copy:
        return (sc.pp.filter_genes(adata, min_counts, min_cells, max_counts,max_cells, inplace, copy))
    else:
        sc.pp.filter_genes(adata, min_counts, min_cells, max_counts,max_cells, inplace, copy)
        

def load_data(path, transpose=False,transforms=[]):
    print("Loading  data ...")
    t0 = time.time()
    if os.path.isdir(path):
        adata = read_mtx(path)
    elif os.path.isfile(path):
        adata = read_csv(path)
    else:
        raise ValueError("File {} not exists".format(path))

    if transpose:
        adata = adata.transpose()

    if type(adata.X) == np.ndarray:
        adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.var_names_make_unique()


    print('Original data contains {} cells x {} peaks'.format(*adata.shape))
    print("Finished loading takes {:.2f} min".format((time.time() - t0) / 60))
    return adata


def read_mtx(path):
    for filename in glob(path + '/*'):
        basename = os.path.basename(filename)
        if (('count' in basename) or ('matrix' in basename)) and ('mtx' in basename):
            count = mmread(filename).T.tocsr().astype('float32')
        elif 'barcode' in basename:
            barcode = pd.read_csv(filename, sep='\t', header=None)[0].values
        elif 'gene' in basename or 'peak' in basename:
            feature = pd.read_csv(filename, sep='\t', header=None).iloc[:, -1].values
    adata = ad.AnnData(count, obs=pd.DataFrame(index=barcode),
                       var=pd.DataFrame(index=peaks))
    return adata


def read_csv(path):
    if ('.txt' in path) or ('tsv' in path):
        sep = '\t'
    elif '.csv' in path:
        sep = ','
    else:
        raise ValueError("File {} not in format txt or csv".format(path))
    data = pd.read_csv(path, sep=sep, index_col=0).T.astype('float32')
    genes = data.columns.values
    barcode = data.index.values
    adata = ad.AnnData(scipy.sparse.csr_matrix(data.values), obs=pd.DataFrame(index=barcode),
                       var=pd.DataFrame(index=genes))
    return adata






