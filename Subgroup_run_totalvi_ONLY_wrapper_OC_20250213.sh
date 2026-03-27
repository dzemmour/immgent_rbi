#!/bin/bash
#SBATCH -p medium
#SBATCH -t 1-11:59:59
#SBATCH --mem 128G
#SBATCH -c 8
#SBATCH -o wrapper_totalvi.log
#SBATCH -e wrapper_totalvi.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch run_totalvi_wrapper_OC_20250213.sh .
# from working directory

module load gcc/14.2.0 python/3.13.1 hdf5/1.14.5 boost/1.87.0 openmpi/4.1.8 fftw/3.3.10 java/jdk-23.0.1 conda/miniforge3/24.11.3-0
#module load python
#source $(conda info --base)/etc/profile.d/conda.sh

export NUMBA_CACHE_DIR="/tmp"
#source activate /project/zemmour/david/envs/scvi120_20241008
#source ~/scvi-tools_20241105/bin/activate
#source /n/groups/cbdm_lab/odc180/Python/envs/scvi_tools_py3.9_250318/bin/activate
source /n/groups/cbdm_lab/odc180/Python/envs/RHEL_env/scvi_tools_pymde_annoy_py3.9_20250522/bin/activate

SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Scripts_TRBI_subgroup/
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_Laurent/Subgroup/
#path_to_anndata=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD4/$dataset/MERGED_RNA_${dataset}_query_ImmgenT.h5ad
path_to_matrix=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_Laurent/matrix/
path_to_samples=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Integration_Webpage/Trial_Laurent/cell_batch.csv
path_to_ImmgenT=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/new_annotations/igt1_104_withtotalvi20260212_downsampled.h5ad
prefix=Trial_Final
#prefix_SCANVI=Trial_SCANVI_Combined_Fix
batchkey=IGT
#categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
subgroup=CD8
cd $working_dir

mkdir $prefix

python3 $SCRIPT_DIR/run_SCVI_SCANVI_subgroup.py --working_dir=$working_dir --path_to_matrix=$path_to_matrix --path_to_ImmgenT=$path_to_ImmgenT --query_sample_path=$path_to_samples --prefix=$prefix --batchkey=$batchkey --corrected_counts=$corrected_counts --denoised_data=$denoised_data --subgroup=$subgroup >totalvi.log 2>totalvi.err

mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/$subgroup/mde_incremental_${subgroup}.csv
totalvi_integrated_file=$prefix/latent_level2.csv
python3 $SCRIPT_DIR/run_mde_incremental_OC_subgroup.py $working_dir $prefix $mde_ref_file $totalvi_integrated_file >mde.log 2>mde.err

## Plot results
output_file=$prefix/predictions_output_file.csv
user_output_file=$prefix/user_predictions_output_file.csv
mde_incremental=$prefix/mde_incremental.csv
query_IGTHT=$prefix/query_all_metadata.csv
path_to_anndata_not_classified=$prefix/adata_RNA.h5ad
SCRIPT_DIR_SUBGROUP=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/clean_not_classified_run/subgroup/
python3 $SCRIPT_DIR/run_SCVI_SCANVI_not_classified_subgroup.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata_not_classified --metadata_query=$output_file --prefix=$prefix --batchkey=$batchkey --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi_not_classified.log 2>totalvi_not_classified.err

conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

annotation_column=IGT
predictions_output_file_not_classified=$prefix/predictions_output_file_not_classified.csv
#subgroup=CD4
mde_plot=$SCRIPT_DIR/mde_plot_new.csv
#EXP_name=$dataset
#so_path=$working_dir/${dataset}_seurat_object.Rds
Rscript $SCRIPT_DIR/Create_final_seurat_object_multipage_subgroup.R $path_to_matrix $query_IGTHT $predictions_output_file_not_classified $mde_incremental $mde_plot $annotation_column $subgroup $prefix >make_final_seurat_objectFinal_separate_SCVIs.log 2>make_final_seurat_objectFinal_separate_SCVIs.err

