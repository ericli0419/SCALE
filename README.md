# Single-Cell ATAC-seq analysis via Latent feature Extraction

## Description

This is a SCANPY compatible version of SCALE based on [SCALE](https://github.com/jsxlei/SCALE.git).

## Installation  

#### install from GitHub

	git clone https://github.com/ericli0419/SCALE.git
	cd SCALE
	python setup.py install

This SCALE version may conflict with the original version of [SCALE](https://github.com/jsxlei/SCALE.git), it is recommended to uninstall [SCALE](https://github.com/jsxlei/SCALE.git) before install this SCALE version.

## Quick Start

#### Input
* either a **count matrix file**:  
	* row is peak and column is barcode, in **txt** / **tsv** (sep=**"\t"**) or **csv** (sep=**","**) format
* or a **folder** contains **three files**:   
	* **count file**: count in **mtx** format, filename contains key word **"count"** / **"matrix"**    
	* **peak file**: 1-column of peaks **chr_start_end**, filename contains key word **"peak"**  
	* **barcode file**: 1-column of barcodes, filename contains key word **"barcode"**

#### Run
with known cluster number k:  

    SCALE.py -d [input] -k [k]

with estimated cluster number k by SCALE if k is unknown: 

    SCALE.py -d [input]

#### Output
Output will be saved in the output folder including:
* **model.pt**:  saved model to reproduce results cooperated with option --pretrain
* **adata.h5ad**: anndata file contains filterted count matrix, cluster assignmnets, latent features and numerical /binary imputated data.
* **emb_tsne/UMAP.txt**:  2d t-SNE/UMAP embeddings of each cell
* **emb_tsne/UMAP.pdf**:  visualization of 2d t-SNE/UMAP embeddings of each cell

#### Imputation  
Get binary imputed data in binary layers of adata.h5ad with option **--binary** (recommended for saving storage)

    SCALE.py -d [input] --binary  

or get numerical imputed data in imputed layers of adata.h5ad with option **--impute**

    SCALE.py -d [input] --impute

#### Useful options  
* save results in a specific folder: [-o] or [--outdir] 
* embed feature by tSNE or UMAP: [--emb]  TSNE/UMAP
* filter rare peaks if the peaks quality if not good or too many, default is 0.01: [-x]
* filter low quality cells by valid peaks number, default 100: [--min_peaks]  
* modify the initial learning rate, default is 0.002: [--lr]  
* change iterations by watching the convergence of loss, default is 30000: [-i] or [--max_iter]  
* change random seed for parameter initialization, default is 18: [--seed]
* binarize the imputation values: [--binary]
* run with scRNA-seq dataset: [--log_transform]
	
#### Note    
If come across the nan loss, 
* try another random seed
* filter peaks with harsher threshold, e.g. -x 0.04 or 0.06
* filter low quality cells, e.g. --min_peaks 400 or 600
* change the initial learning rate, e.g. --lr 0.0002 
	

#### Help
Look for more usage of SCALE

	SCALE.py --help 

Use functions in SCALE packages.

	import scale
	from scale import *
	from scale.plot import *
	from scale.utils import *

#### Data availability  

* [Forebrain](http://zhanglab.net/SCALE_SOURCE_DATA/Forebrain.tar)
* [Splenocyte](http://zhanglab.net/SCALE_SOURCE_DATA/Splenocyte.tar)
* [mouse_atlas](http://zhanglab.net/SCALE_SOURCE_DATA/mouse_atlas.tar)
* [InSilico](http://zhanglab.net/SCALE_SOURCE_DATA/InSilico.tar)
* [Leukemia](http://zhanglab.net/SCALE_SOURCE_DATA/Leukemia.tar)
* [GM12878vsHEK](http://zhanglab.net/SCALE_SOURCE_DATA/GM12878vsHEK.tar)
* [GM12878vsHL](http://zhanglab.net/SCALE_SOURCE_DATA/GM12878vsHL.tar)
* [Breast_Tumor](http://zhanglab.net/SCALE_SOURCE_DATA/Breast_Tumor.tar)


## Reference
[Lei Xiong, Kui Xu, Kang Tian, Yanqiu Shao, Lei Tang, Ge Gao, Michael Zhang, Tao Jiang & Qiangfeng Cliff Zhang. SCALE method for single-cell ATAC-seq analysis via latent feature extraction. Nature Communications, (2019).](https://www.nature.com/articles/s41467-019-12630-7)
