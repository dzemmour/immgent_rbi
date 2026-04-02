#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
#parser.add_argument('--path_to_anndata', help='Path to combined AnnData object')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--path_to_ImmgenT', help='Path to ImmgenT MuData object')
parser.add_argument('--path_to_matrix', help='Path to query matrix')
parser.add_argument('--query_sample_path', help='Path to query sample table')
parser.add_argument('--subgroup', default=None, help='Subgroup for integration')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
#path_to_anndata = args.path_to_anndata
path_to_ImmgenT = args.path_to_ImmgenT
path_to_matrix = args.path_to_matrix
query_sample_path = args.query_sample_path
subgroup = args.subgroup
prefix = args.prefix
batchkey = args.batchkey
#confoundings = args.confoundings
if args.categorical_covariate_keys:
        # Split the string into a list by commas
        categorical_covariate_keys = args.categorical_covariate_keys.split(',')
        #print(categorical_covariate_keys)
else:
    categorical_covariate_keys = None
corrected_counts = args.corrected_counts
denoised_data = args.denoised_data
#mde_ref_file = args.mde_ref_file
#latent_key = args.latent_key

# working_dir = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/'
# path_to_mudata = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu'
# prefix = 'totalvi_igt1_56_allgenes_20240526_igtsampleregressedout'
# batchkey = 'IGT'
# categorical_covariate_keys = ['sample_id']

# working_dir='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96'
# path_to_mudata='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241006.h5mu'
# prefix='totalvi_20241006' #prefix = "totalvi_20241008_rmIGTsample"
# batchkey='IGT'
# categorical_covariate_keys='IGTHT'
# corrected_counts=False
# denoised_data=False


print(f"Working Directory: {working_dir}")
#print(f"Path to ImmgenT AnnData: {path_to_ImmgenT}")
#print(f"Path to query AnnData: {path_to_query}")
#print(f"Path to combined AnnData: {path_to_anndata}")
print(f"Path to ImmgenT AnnData: {path_to_ImmgenT}")
print(f"Path to matrix: {path_to_matrix}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
#print(f"mde_ref_file: {mde_ref_file}")
#print(f"Latent Key: {latent_key}")

print("Importing libraries")
import warnings; warnings.simplefilter('ignore')
import scvi
import os
import sys
import scanpy as sc
#import muon as mu
import mudata as mu
import anndata as AnnData
import numpy as np
#import mplscience
import scipy
import scipy.io
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import fast_matrix_market as fmm
from annoy import AnnoyIndex
from collections import Counter
import re
from pathlib import Path
#import pymde #to run MDE
#from annoy import AnnoyIndex
#from collections import Counter

print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
scvi.settings.seed = 0  # optional: ensures reproducibility
#pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("Reading mudata")
os.chdir(working_dir)
## ImmgenT
ImmgenT_mdata_all = AnnData.read(path_to_ImmgenT)
## Keep only IGT1_96
# Define the values to remove
exclude_IGT_values = [f'IGT{i}' for i in range(97, 105)]  # 105 is exclusive
All_IGT_values = [f'IGT{i}' for i in range(1, 97)]
# Filter out those cells
ImmgenT_mdata_all = ImmgenT_mdata_all[~ImmgenT_mdata_all.obs['IGT'].isin(exclude_IGT_values)].copy()
ImmgenT_mdata = ImmgenT_mdata_all[ImmgenT_mdata_all.obs['level1'] == subgroup, :].copy()
ImmgenT_cells = ImmgenT_mdata.obs_names.astype(str)
ImmgenT_mdata.X = ImmgenT_mdata.X.copy()
ImmgenT_mdata.layers["counts"] = ImmgenT_mdata.X.copy()
print(ImmgenT_mdata)
## Read in query Anndata object
#query_mdata = AnnData.read(path_to_query)

print("Create anndata object from matrix")
##matrix = scipy.io.mmread(path_to_matrix+"/matrix.mtx.gz").tocsr()
##matrix = matrix.transpose().tocsr()
##genes = pd.read_csv(path_to_matrix+"/features.tsv.gz", sep = "\t", header=None, compression='infer')
##cells = pd.read_csv(path_to_matrix+"/barcodes.tsv.gz", sep = "\t", header=None, compression='infer')

p = Path(path_to_matrix)

def find_file(stem, exts=("", ".gz")):
    for ext in exts:
        f = p / f"{stem}{ext}"
        if f.exists():
            return f
    raise FileNotFoundError(f"Could not find {stem}[.gz] in {p}")

# locate files
matrix_file   = find_file("matrix.mtx")
#barcodes_file = find_file("barcodes.tsv")
barcodes_file = (
    find_file("barcodes.tsv") if (p / "barcodes.tsv").exists() or (p / "barcodes.tsv.gz").exists()
    else find_file("cells.tsv")
)

features_file = (
    find_file("features.tsv") if (p / "features.tsv").exists() or (p / "features.tsv.gz").exists()
    else find_file("genes.tsv")
)

# read data
matrix = scipy.io.mmread(matrix_file).tocsr().T
cells  = pd.read_csv(barcodes_file, sep="\t", header=None)
genes  = pd.read_csv(features_file, sep="\t", header=None)
genes.columns = ['gene_name'] #, 'feature_type']
# Make genes column unique
def make_unique(names):
    counts = {}
    result = []
    for name in names:
        if name not in counts:
            counts[name] = 1
            result.append(name)
        else:
            new_name = f"{name}_{counts[name]}"
            while new_name in counts:
                counts[name] += 1
                new_name = f"{name}_{counts[name]}"
            result.append(new_name)
            counts[new_name] = 1
    return result

# Apply to gene_name column
genes['gene_name'] = make_unique(genes['gene_name'])

cells.columns = ['barcode']

query_mdata = AnnData.AnnData(X=matrix,
        obs=pd.DataFrame(index=cells["barcode"].values),
        var=pd.DataFrame(index=genes["gene_name"].values))
#del mdata
query_mdata.X = query_mdata.X.copy()
query_mdata.layers["counts"] = query_mdata.X.copy()
print("Pre gene filtering: ")
print(query_mdata)
sc.pp.calculate_qc_metrics(query_mdata, inplace=True)
cells_nGenes = query_mdata.obs.index[(query_mdata.obs['n_genes_by_counts'] < 300)]
sc.pp.filter_cells(query_mdata, min_genes=300)
print("Post gene filtering: ")
print(query_mdata)
## Add sample and hashtag information to query mdata before spike-in merge
query_sample_table = pd.read_csv(query_sample_path, index_col=0)
#query_sample_table.index = query_sample_table["Cell_ID"]
joint_bcs = query_sample_table.index.intersection(query_mdata.obs_names)
query_sample_table = query_sample_table.loc[joint_bcs, :].copy()
query_mdata = query_mdata[joint_bcs, :].copy()
query_mdata.obs["IGT"] = query_sample_table["batch"]
#query_mdata.obs["HT"] = query_sample_table["sample"]
#query_mdata.obs["IGTHT"] = query_mdata.obs["IGT"].astype(str) + "_" + query_mdata.obs["HT"].astype(str)
query_mdata.obs.to_csv(prefix+"/query_all_metadata.csv", index = True)
query_cells = query_mdata.obs_names.astype(str)
print("Save Anndata")
query_mdata.write_h5ad(prefix+"/adata_RNA.h5ad")

print("Concat anndata - keep only common genes")
mdata = AnnData.concat(
    [ImmgenT_mdata, query_mdata],
    axis=0,
    join="inner",
    label="origin",
    keys=["ImmgenT", "query"]
)


##mdata.obs['IGTHT'] = pd.concat([ImmgenT_IGTHT["IGTHT"],query_IGTHT["IGTHT"]]).astype(str)
##mdata.obs['IGT'] = pd.concat([ImmgenT_IGTHT["IGT"],query_IGTHT["IGT"]]).astype(str)
mdata.layers["counts"] = mdata.X.copy()
mdata.obs_names = mdata.obs_names.astype(str)
mdata.var_names = mdata.var_names.astype(str)
mdata.X = mdata.X.copy()
##mdata.var = mdata.var.reset_index(drop=True)
##mdata.mod['RNA'].var = mdata.mod['RNA'].var.reset_index(drop=True)
##mdata.mod['protein'].var = mdata.mod['protein'].var.reset_index(drop=True)
print(mdata)
##mdata.update()

print("Batch key list")
print(mdata.obs[batchkey].unique().tolist())
print("Batch key list: any NA?")
print((mdata.obs[batchkey].isna()).any())


print("categorical_covariate_keys")
if categorical_covariate_keys is not None:
    for c in categorical_covariate_keys:
        print(c)
        print(mdata.obs[c].unique().tolist())
        print("any NA?")
        print((mdata.obs[c].isna()).any())
        ##print(mdata.mod["RNA"].obs[c].head(10))
        ##mdata.mod["RNA"].obs[c] = mdata.mod["RNA"].obs[c].str.replace('.', '', regex=False)
        ##print(mdata.mod["RNA"].obs[c].head(10))

#ImmgenT_IGTHT = pd.read_csv("ImmgenT_IGTHT.csv", index_col = 0)
#ImmgenT_cells = ImmgenT_IGTHT.index
## Normalize first to use scoring
for b in mdata.obs['IGT'].unique():
    sc.pp.normalize_total(mdata[mdata.obs['IGT'] == b], target_sum=1e4)
    sc.pp.log1p(mdata[mdata.obs['IGT'] == b])

## Add T Cell score and filter based on this
## Read in gene signature list - select T cell genes
signature_list = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/LineageSpecGns072018_top27.csv", header = 0)
Tcell_genes = signature_list.loc[81:107, "Marker"].astype(str).tolist()
#Tcell_genes_clean = [g for g in Tcell_genes if g in mdata.var_names]
print(Tcell_genes)

# 3. Compute gene scores (after normalization but before batch correction)
sc.tl.score_genes(mdata, gene_list=Tcell_genes, score_name='Tcell_gene_module_score')

print("Pre Tcell filter: ")
print(mdata)
Tcell_mask = (
    (mdata.obs_names.isin(query_cells) & (mdata.obs['Tcell_gene_module_score'] > 0))
    | mdata.obs_names.isin(ImmgenT_cells)
)

non_Tcell_mask = (
    (mdata.obs_names.isin(query_cells) & (mdata.obs['Tcell_gene_module_score'] <= 0))
)

cells_Tcell_score = mdata.obs.index[non_Tcell_mask]
#mdata = mdata[Tcell_mask, :].copy()

print("save QC filtered cells (if any)")
dfs = []

if cells_Tcell_score is not None:
    #dfs.append(pd.DataFrame({"cell_id": cells_Tcell_score, "QC_filter": "non T cell"}))

if cells_nGenes is not None:
    dfs.append(pd.DataFrame({"cell_id": cells_nGenes, "QC_filter": "low gene count"}))

if dfs:
    filtered_cells = pd.concat(dfs, axis=0, ignore_index=True)
    filtered_cells.to_csv(prefix+"/QC_filtered_cells.csv", index=False)

print("Remaining Tcells: ")
print(mdata)

print("Subgroup filter: ")


print("Train new SCVI model")
#annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/annotation_table_20250505_IGT1_104.csv", delimiter = ",", index_col=0)
annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/new_annotations/GSE297097_annotation_table_20260206_IGT1_104_cleaned.csv", delimiter = ",", index_col=0)
annotation_level2 = annotation_level2.loc[annotation_level2.index.intersection(mdata.obs.index)]
mdata.obs['level1'] = "Unknown"
mdata.obs['level2'] = "Unknown"
##mdata.obs['level2.group'] = "Unknown"
mdata.obs.loc[annotation_level2.index, 'level1'] = annotation_level2['level1'].values
mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values
##mdata.obs.loc[annotation_level2.index, 'level2.group'] = annotation_level2['level2.group'].values

print("LABELS:")
print(mdata.obs["level2"].value_counts())

print("\nBATCHES:")
print(mdata.obs[batchkey].value_counts())

#print("\nCOVARIATES:")
#for cov in categorical_covariate_keys:
#    print(cov)
#    print(mdata.obs[cov].value_counts())

#df_counts = mdata.obs[categorical_covariate_keys].value_counts()
#single_values = df_counts[df_counts == 1].index.tolist()
#mdata.obs["IGTHT"] = mdata.obs["IGTHT"].astype("category").cat.add_categories("merged_batch_ref")
#mdata.obs.loc[mdata.obs["IGTHT"].isin(single_values), "IGTHT"] = "merged_batch_ref"

#print("\nCOVARIATES Singletons:")
#print(single_values)

# Comment if already trained
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model = scvi.model.SCVI(mdata, n_latent=30, n_layers=2)
scvi_model.train()
scvi_model.save(prefix+"/scvi_model/", save_anndata=True)
#scvi_model = scvi.model.SCVI.load(prefix+"/scvi_model/", adata=mdata)



SCANVI_LATENT_KEY = "X_scANVI"
print("Save latent_representation.csv")
#latent_representation = scvi_model.get_latent_representation()
SCVI_LATENT_KEY = "X_scVI"
##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
#mdata.obsm[SCVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
#latent_df_scvi = pd.DataFrame(latent_representation, index = mdata.obs.index)

#latent_df_scvi.to_csv(prefix+"/latent_scvi.csv", index=True)

print("Save umap.csv")
#sc.pp.neighbors(mdata, use_rep=SCVI_LATENT_KEY)
#sc.tl.umap(mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
#umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.obs.index)
#umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save Anndata")
#mdata.write_h5ad(prefix+"/adata_RNA.h5ad")

##mdata = AnnData.read(prefix+"_orig/adata_RNA.h5ad")
##latent_df = pd.read_csv(prefix+"_orig/latent.csv", index_col=0)

print("Predict level2 annotations with SCANVI")
## Define ref and query masks (future use)
#ref_mask = mdata.obs["origin"] == "ImmgenT"
#query_mask = mdata.obs["IGT"] == "10X_sc"
query_mask = (~mdata.obs["IGT"].str.contains("IGT", na=False))

## level2
## Training scanvi model on scvi model
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys, labels_key="level2")

#print("LABELS:")
#print(mdata.obs["level2"].value_counts())

#print("\nBATCHES:")
#print(mdata.obs[batchkey].value_counts())

#print("\nCOVARIATES:")
#for cov in categorical_covariate_keys:
#    print(cov)
#    print(mdata.obs[cov].value_counts())

#df_counts = mdata.obs[categorical_covariate_keys].value_counts()
#single_values = df_counts[df_counts == 1].index.tolist()
#mdata.obs["IGTHT"] = mdata.obs["IGTHT"].astype("category").cat.add_categories("merged_batch_ref")
#mdata.obs.loc[mdata.obs["IGTHT"].isin(single_values), "IGTHT"] = "merged_batch_ref"


n_labeled = (mdata.obs["level2"] != "not classified").sum()
print("Labeled level2 cells:", n_labeled)
if n_labeled < 2:
    print("Skipping level2 SCANVI: <2 labeled cells after filtering")
else:
    level2_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown", labels_key="level2")
    level2_model.train(50, train_size=0.9, batch_size=128, validation_size=0.1, datasplitter_kwargs={"drop_last": True})
    #level2_model = scvi.model.SCANVI.from_scvi_model(scvi_model,"Unknown", labels_key="level2")
    #level2_model.train(50)
    level2_model.save(prefix+"/scanvi_level2_model/", save_anndata=True)

#level2_model = scvi.model.SCANVI.from_scvi_model(scvi_model,"Unknown", labels_key="level2")
#level2_model.train(50)
#level2_model.save(prefix+"/scanvi_level2_model/", save_anndata=True)
#level2_model = scvi.model.SCANVI.load(prefix+"/scanvi_level2_model/", adata=mdata)


## Predictions and scores - create output file
## level2
##SCVI_LATENT_KEY = "X_SCVI"
LEVEL2_SCANVI_LATENT_KEY = "level2_X_scANVI"
LEVEL2_SCANVI_PREDICTIONS_KEY = "level2_C_scANVI"

##mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(mdata)
mdata.obsm[LEVEL2_SCANVI_LATENT_KEY] = level2_model.get_latent_representation(mdata)

mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model.predict(mdata)
output_file = pd.DataFrame(mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY], index = mdata.obs.index)
latent_df = pd.DataFrame(mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = mdata.obs.index)
latent_df.to_csv(prefix+"/latent_level2.csv", index=True)

## Get posterior probabilities for all labels
level2_probs = level2_model.predict(mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
level2_confidence = level2_probs.max(axis=1)
## Add to AnnData
output_file["level2_scanvi_confidence"] = level2_confidence
##output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
confidence_threshold = 0.85
output_file["level2_final"] = output_file[LEVEL2_SCANVI_PREDICTIONS_KEY]
#output_file.loc[output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "unclear"
output_file.loc[output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "not classified"


## Save user annotations to return to the user
##output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels']] #,'level2.group_transfer_labels']]
##output_file.index = mdata.obs.index.copy()
##output_file.to_csv(prefix+"/output_annotations.csv", index=True)
output_file.to_csv(prefix+"/predictions_output_file.csv", index=True)
user_output_file = output_file.loc[query_mask, :]
##user_output_file.index = query_mdata.obs.index
user_output_file.to_csv(prefix+"/user_predictions_output_file.csv", index=True)


