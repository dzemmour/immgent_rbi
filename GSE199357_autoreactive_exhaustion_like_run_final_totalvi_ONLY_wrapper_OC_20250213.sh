#!/bin/bash
#SBATCH -p medium
#SBATCH -t 1-11:59:59
#SBATCH --mem 64G
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


SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/clean_not_classified_run/
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/Webpage_Trial/
path_to_ImmgenT=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/ImmgenT_downsampled_20250505_adata.h5ad
path_to_query=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/query_adata.h5ad
path_to_anndata=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/MERGED_RNA_GSE199357_autoreactive_exhaustion_like_query_ImmgenT.h5ad
path_to_spikein=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/RNA_Integration_15k_Downsampled.h5ad
prefix=Not_classified_Full_Trial
#prefix_SCANVI=Trial_SCANVI_Combined_Fix
batchkey=IGT
categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
cd $working_dir

python3 $SCRIPT_DIR/run_SCVI_SCANVI.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata --path_to_ImmgenT=$path_to_ImmgenT --path_to_spikein=$path_to_spikein --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi.log 2>totalvi.err

path_to_anndata_not_classified=$prefix/adata_RNA.h5ad
predictions_output_file=$prefix/predictions_output_file.csv
python3 $SCRIPT_DIR/run_SCVI_SCANVI_not_classified.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata --metadata_query=$predictions_output_file --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi_not_classified.log 2>totalvi_not_classified.err

predictions_output_file_not_classified=$prefix/predictions_output_file_not_classified.csv
latent_level1_abT=$prefix/latent_level1_abT.csv
latent_df_scvi_abT=$prefix/latent_level2_abT.csv
latent_df_scvi_gdT=$prefix/latent_level2_gdT.csv
mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_downsampled_level1_and_level2_mde.csv
python3 $SCRIPT_DIR/run_mde_incremental_OC.py $working_dir $prefix $mde_ref_file $latent_df_scvi_abT $latent_df_scvi_gdT $latent_level1_abT $predictions_output_file_not_classified >mde.log 2>mde.err

conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

annotation_column=IGTHT
query_IGTHT=$working_dir/../query_all_metadata.csv
user_output_file=$prefix/user_output_file.csv
output_file=$prefix/output_file.csv
mde_plot=$SCRIPT_DIR/mde_plot.csv
Rscript $SCRIPT_DIR/Plotting_for_webpage_relax_threshold_samples_annotation.R $output_file $output_file $prefix $query_IGTHT $annotation_column $mde_plot >plot_samples_Final_separate_SCVIs.log 2>plot_samples_Final_separate_SCVIs.err

so_path=$working_dir/../GSE199357_autoreactive_exhaustion_like_seurat_object.Rds
output_file=$prefix/output_file.csv
Rscript $SCRIPT_DIR/Create_final_seurat_object_multipage.R $so_path $query_IGTHT $output_file $mde_plot $annotation_column $prefix $EXP_name  >make_final_seurat_objectFinal_separate_SCVIs.log 2>make_final_seurat_objectFinal_separate_SCVIs.err
