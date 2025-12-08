#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_anndata', help='Path to combined AnnData object')
#parser.add_argument('--output_file', help='Path to output file query')
parser.add_argument('--metadata_query', help='Path to metadata query')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_anndata = args.path_to_anndata
#output_file_path = args.output_file
metadata_query_path = args.metadata_query
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
print(f"Path to combined AnnData: {path_to_anndata}")
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
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import fast_matrix_market as fmm
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
#mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
#print(mdata)
mdata = AnnData.read(path_to_anndata) #load this one (results of this block)
mdata = mdata[~mdata.obs_names.str.contains("IGT")].copy()
metadata_query = pd.read_csv(metadata_query_path, index_col=0)
metadata_query = metadata_query[~metadata_query.index.str.contains("IGT")].copy()
common_id = mdata.obs_names.intersection(metadata_query.index)
mdata = mdata[common_id, :].copy()
metadata_query = metadata_query.loc[common_id, :].copy()
#mdata.obs.loc[output_file.index,  :] = output_file.loc[output_file.index, :].values
#mdata.obs = metadata_query
mdata.obs[["level2_C_scANVI"]] = metadata_query[["level2_C_scANVI"]]
mdata.obs[["level2_final"]] = metadata_query[["level2_final"]]
mdata.obs[["level2_scanvi_confidence"]] = metadata_query[["level2_scanvi_confidence"]]
mdata.obs.loc[mdata.obs["level2_final"] == "unclear", "level2_final"] = "not classified"
mdata.raw = None
print(mdata)

print("Creating RNA layer counts")
mdata.layers["counts"] = mdata.X.copy()

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

print("Train new SCVI model")
## add final classifications from medata query
print("metadata_query shape: ",  metadata_query.shape)
print("mdata: ")
print(mdata)
## use original output file

## add not classified percentage

## start for loop

## not just level1 and level2


## check percentages
#...
## < 1% not classified stop
## if > 1% and if % not classified is not prev % /2 (not 50% improvement), then stop
# if < 1% or no 50% improvement, stop
# if > 1% and 50% improvement, continue
initial_percentage = (mdata.obs['level2_final'].eq('not classified').sum() / len(mdata.obs)) * 100
original_percentage = 100
initial_difference = (original_percentage/initial_percentage)

print("original percentage and difference: ")
print(initial_percentage)
print(initial_difference)

if (initial_percentage < 1 or initial_difference < 1.5):
    output_file = metadata_query

percentage = initial_percentage
difference = initial_difference

## Separate gdT and abT mdata
#abT_mdata = mdata[mdata.obs['level1_final'].astype(str) != "gdT" , :].copy()
#gdT_mdata = mdata[mdata.obs['level1_final'] == "gdT" ,:].copy()

i = 1
while ((percentage > 1) and (difference > 1.5)):

    print("loop start : ", i)
    print(percentage)
    print(difference)

    # Comment if already trained
    scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
    scvi_model = scvi.model.SCVI(mdata, n_latent=30, n_layers=2)
    scvi_model.train()
    #scvi_model.save(prefix+"/scvi_model_not_classified/", save_anndata=True)
    #scvi_model = scvi.model.SCVI.load(prefix+"/scvi_model_not_classified/", adata=mdata)


    SCANVI_LATENT_KEY = "X_scANVI"
    print("Save latent_representation.csv")
    latent_representation = scvi_model.get_latent_representation()
    SCVI_LATENT_KEY = "X_scVI"
    ##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
    mdata.obsm[SCVI_LATENT_KEY] = latent_representation
    ##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
    latent_df_scvi = pd.DataFrame(latent_representation, index = mdata.obs.index)
    #latent_df_scvi.to_csv(prefix+"/latent_scvi_not_classified.csv", index=True)

    print("Save umap.csv")
    sc.pp.neighbors(mdata, use_rep=SCVI_LATENT_KEY)
    sc.tl.umap(mdata, min_dist=0.4)
    ##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
    umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.obs.index)
    #umap_df.to_csv(prefix+"/umap_python_not_classified.csv", index=True)

    print("Save Anndata")
    #mdata.write_h5ad(prefix+"/adata_RNA_not_classified.h5ad")

    ## level2
    ## Training scanvi model on scvi model
    mdata.obs["level2_final"] = mdata.obs["level2_final"].astype(str)
    labels_to_keep = mdata.obs["level2_final"].unique().astype(str)
    labels_to_keep = [x for x in labels_to_keep if x != "not classified"]
    mdata.obs["level2_final"] = mdata.obs["level2_final"].where(
            mdata.obs["level2_final"].isin(labels_to_keep),
            other="not classified"
            )
    mdata.obs["level2_final"] = mdata.obs["level2_final"].astype("category")

    scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys, labels_key="level2_final")
    level2_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "not classified", labels_key="level2_final")
    level2_model.train(30)
    #level2_model.save(prefix+"/scanvi_level2_model_not_classified/", save_anndata=True)
    #level2_model = scvi.model.SCANVI.load(prefix+"/scanvi_level2_model_not_classified/", adata=mdata)

    ## Predictions and scores - create output file
    ## level2
    ##SCVI_LATENT_KEY = "X_SCVI"
    LEVEL2_SCANVI_LATENT_KEY = "level2_X_scANVI"
    LEVEL2_SCANVI_PREDICTIONS_KEY = "level2_C_scANVI"

    ##mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(mdata)
    mdata.obsm[LEVEL2_SCANVI_LATENT_KEY] = level2_model.get_latent_representation(mdata)
    mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model.predict(mdata)
    #output_file = pd.concat([output_file, pd.DataFrame(mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY])], axis = 1)
    output_file = pd.DataFrame(mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY], index = mdata.obs_names)
    latent_df = pd.DataFrame(mdata.obsm[LEVEL2_SCANVI_LATENT_KEY], index = mdata.obs.index)
    #latent_df.to_csv(prefix+"/latent_level2_not_classified.csv", index=True)

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
    output_file.loc[output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "not classified"
    
    overlap = output_file.index.intersection(mdata.obs.index)
    ##level2
    convert_level2 = overlap[
            (mdata.obs.loc[overlap,"level2_final"]=="not classified") &
            (output_file.loc[overlap,"level2_final"]!="not classified") &
            (output_file.loc[overlap, "level2_scanvi_confidence"] > 0.85)
            ]

    ## Only overwrite accepted updates
    mdata.obs.loc[convert_level2, ["level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]] = output_file.loc[convert_level2, ["level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]].values
    
    output_file = mdata.obs[["level2_C_scANVI", "level2_scanvi_confidence", "level2_final"]]

    percentage_new = (output_file["level2_final"].eq("not classified").sum() / len(output_file.index)) * 100
    difference = (percentage/percentage_new)
    percentage = percentage_new

    print("loop end: ", i)
    print(percentage)
    print(difference)

    i = i + 1

    #del scvi_model
    #del level2_model


## check percentages
#...
## < 1% not classified stop
## if > 1% and if % not classified is not prev % /2 (not 50% improvement), then stop
# if < 1% or no 50% improvement, stop
# if > 1% and 50% improvement, continue

#scvi_model = scvi.model.SCVI.load(prefix+"/scvi_model_not_classified/", adata=mdata)
#level2_model = scvi.model.SCANVI.load(prefix+"/scanvi_level2_model_not_classified/", adata=mdata)
level2_model.save(prefix+"/scanvi_level2_model_not_classified/", save_anndata=True)
scvi_model.save(prefix+"/scvi_model_not_classified/", save_anndata=True)


## Save user annotations to return to the user
##output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels']] #,'level2.group_transfer_labels']]
##output_file.index = mdata.obs.index.copy()
##output_file.to_csv(prefix+"/output_annotations.csv", index=True)
output_file.to_csv(prefix+"/predictions_output_file_not_classified.csv", index=True)
user_output_file = output_file.copy()
##user_output_file.index = query_mdata.obs.index
user_output_file.to_csv(prefix+"/user_predictions_output_file_not_classified.csv", index=True)

