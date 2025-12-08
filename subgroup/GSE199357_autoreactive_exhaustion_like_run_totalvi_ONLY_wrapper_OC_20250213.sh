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


SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/clean_not_classified_run/subgroup/
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/
path_to_anndata=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/MERGED_RNA_GSE199357_autoreactive_exhaustion_like_query_ImmgenT.h5ad
prefix=Trial_Final
prefix_SCANVI=Trial_SCANVI_Combined_Fix
batchkey=IGT
categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
cd $working_dir
python3 $SCRIPT_DIR/run_SCVI_SCANVI.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi.log 2>totalvi.err


mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/CD8AB/mde_incremental_CD8AB.csv
totalvi_integrated_file=$prefix/latent_level2.csv
python3 $SCRIPT_DIR/run_mde_incremental_OC.py $working_dir $prefix $mde_ref_file $totalvi_integrated_file >mde.log 2>mde.err

source deactivate

## Plot results
conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env
output_file=$prefix/predictions_output_file.csv
user_output_file=$prefix/user_predictions_output_file.csv
mde_incremental=$prefix/mde_incremental.csv
query_IGTHT=$working_dir/query_all_metadata.csv
Rscript $SCRIPT_DIR/Plot_results.R $output_file $user_output_file $mde_incremental $prefix $query_IGTHT > plot.log 2 > plot.err

#so_path=$working_dir/GSE131847_circulating_tissue_resident_seurat_object.Rds
#Rscript $SCRIPT_DIR/Create_final_seurat_object_subgroup.R $so_path $query_IGTHT $output_file $mde_incremental $prefix >make_final_seurat_objectFinal_separate_SCVIs.log 2>make_final_seurat_objectFinal_separate_SCVIs.err

source /n/groups/cbdm_lab/odc180/Python/envs/RHEL_env/scvi_tools_pymde_annoy_py3.9_20250522/bin/activate

path_to_anndata_not_classified=$prefix/adata_RNA.h5ad
SCRIPT_DIR_SUBGROUP=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/clean_not_classified_run/subgroup/
python3 $SCRIPT_DIR/run_SCVI_SCANVI_not_classified_subgroup.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata --metadata_query=$output_file --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi_not_classified.log 2>totalvi_not_classified.err

conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

annotation_column=IGTHT
predictions_output_file_not_classified=$prefix/predictions_output_file_not_classified.csv
Rscript $SCRIPT_DIR/Plot_results_not_classified.R $output_file $user_output_file $mde_incremental $prefix $query_IGTHT $predictions_output_file_not_classified $annotation_column  > plot_not_classified.log 2 > plot_not_classified.err

annotation_column=IGTHT
subgroup=CD8
mde_plot=$SCRIPT_DIR/../mde_plot.csv
EXP_name=GSE131847_circulating_tissue_resident
so_path=$working_dir/GSE199357_autoreactive_exhaustion_like_seurat_object.Rds
Rscript $SCRIPT/Create_final_seurat_object_multipage_subgroup.R $so_path $query_IGTHT $predictions_output_file_not_classified $mde_incremental $mde_plot $annotation_column $subgroup $prefix $EXP_name >make_final_seurat_objectFinal_separate_SCVIs.log 2>make_final_seurat_objectFinal_separate_SCVIs.err

