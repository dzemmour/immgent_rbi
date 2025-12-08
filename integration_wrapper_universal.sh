#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-1:00:00
#SBATCH --mem 64G
#SBATCH -c 4
#SBATCH -o wrapper_format.log
#SBATCH -e wrapper_format.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch integration_wrapper_universal.sh /n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/ . 
# from Broad Run directory

### Load Relevant Modules
module load gcc/14.2.0 R/4.4.2 hdf5/1.14.5 boost/1.87.0 openmpi/4.1.8 fftw/3.3.10 java/jdk-23.0.1 conda/miniforge3/24.11.3-0 

#export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
#export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
conda activate /n/groups/cbdm_lab/odc180/Python/conda/R_4.4.2_clean_env

### Change to target directory
cd $2

so_path=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/
so_name=GSE199357_autoreactive_exhaustion_like_seurat_object.Rds
output_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE199357_autoreactive_exhaustion_like/
h5_name=GSE199357_autoreactive_exhaustion_like_query_ImmgenT.h5seurat
subgroup=CD8AB


Rscript $1/Create_mixed_anndata_object.R $so_path $so_name $output_dir $h5_name $subgroup >create.log 2>create.err

# set virtual environment
#source ~/scvi-tools_20241105/bin/activate

#python3 $1/Format_merged_mudata.py /n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/New_Model/gdT/Enxhi/ MERGED_RNA_gdT_query_ImmgenT.h5ad MERGED_ADT_gdT_query_ImmgenT.h5ad >format.log 2>format.err



#python3 $1/train_query_model_v3.py /n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/unconventional/IGT_Trial /n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/Treg/ /n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/Treg/ unconventional 2>python.err > python.log

conda deactivate

sbatch $so_path/GSE199357_autoreactive_exhaustion_like_run_totalvi_ONLY_wrapper_OC_20250213.sh
