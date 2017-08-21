################################################################################################
#on the cluster, go to a directory where monocle / xacHelper located:  
################################################################################################
#install the package: 
##make sure you install all the packages listed in the script before you run all the downstream 
#analysis: 
#Rscript install_packages.R ##can this be done well automatically?  

################################################################################################
#download the data necessary for the analysis: 
wget  http://www.gs.washington.edu/~xqiu/proj/BEAM_analysis_data.tar.gz

#untar the file and then a directory named data which included all data necessary for reproducing 
# the BEAM analysis will be provided: 
tar -zxvf BEAM_analysis_data.tar.gz 
rm BEAM_analysis_data.tar.gz

#make the directories to store the figures generatated in the script 
mkdir -p main_figures supplementary_figures supplementary_data tmp RData

################################################################################################
#run all following code in parallel: 

################################################################
#prepare the data for the lung dataset analysis
qsub -N prepare_lung_data -l mfree=100G Rscript prepare_lung_data.R 

##notee that the following dependent on the prepare_lung_data.R but each of them are independent 
#and can be parallelized

#perform the BEAM analysis for the lung dataset
qsub -N analysis_lung_data -hold_jid prepare_lung_data -l mfree=100G Rscript analysis_lung_data.R 
#run the analysis using the recovery counts
qsub -N mc_analysis_lung_data -hold_jid prepare_lung_data -l mfree=100G Rscript analysis_lung_data_mc.R 

#perform the analysis for the spike-in free algorithm
qsub -N spikein_free_algorithm_sampling -hold_jid prepare_lung_data -l mfree=100G Rscript spikein_free_algorithm_sampling.R 
#perform benchmark analysis
qsub -N deg_benchmark_analysis -hold_jid prepare_lung_data -l mfree=100G Rscript deg_benchmark_analysis.R
################################################################

################################################################
#perform the analysis for the HSMM dataset
qsub -N analysis_HSMM_data -l mfree=100G Rscript analysis_HSMM_data.R 
################################################################

################################################################
#perform the analysis for the UMI dataset 
qsub -N analysis_UMI_data -l mfree=100G Rscript analysis_UMI_data.R 
################################################################

################################################################
#perform the analysis for the Shalek dataset
qsub -N analysis_shalek_data -l mfree=100G Rscript analysis_shalek_data.R 
################################################################

################################################################
#perform the analysis for making supplementary figures 
qsub -N analysis_distribution_fitting -hold_jid prepare_lung_data -l mfree=100G Rscript analysis_distribution_fitting.R 
################################################################


################################################################################################
#when the objects are generated from the above run, the following scripts 
#can be used to generate figures in the manuscript
################################################################################################

################################################################################################
#the following script can be run in parallel or separately: 
#generate the figures: 
qsub -N gen_lung_figures  -hold_jid analysis_lung_data,spikein_free_algorithm_sampling,deg_benchmark_analysis -l mfree=100G Rscript gen_lung_figures.R
qsub -N mc_gen_lung_figures  -hold_jid mc_analysis_lung_data,spikein_free_algorithm_sampling,deg_benchmark_analysis -l mfree=100G Rscript gen_lung_figures_mc.R
qsub -N gen_shalek_figures -hold_jid analysis_shalek_data -l mfree=100G Rscript gen_shalek_figures.R
qsub -N gen_supplementary_figure -hold_jid analysis_lung_data,spikein_free_algorithm_sampling,deg_benchmark_analysis,analysis_distribution_fitting,analysis_HSMM_data -l mfree=100G Rscript gen_supplementary_figure.R 

###### The following are scripts developed for the first resubmission #####
#The following code only includes scripts for new figures added into the manuscript
#The reviewer figures analysis scripts are included but not incorporated 
#########################################################

## Analysis cell downsampling in parallel
qsub -N A_analysis_cell_downsampling -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_A.R 
qsub -N B_analysis_cell_downsampling -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_B.R 
qsub -N C_analysis_cell_downsampling -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_C.R 

qsub -N 1_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_1.R 
qsub -N 2_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_2.R 
qsub -N 3_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_3.R 
qsub -N 4_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_4.R 
qsub -N 5_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_5.R 
qsub -N 6_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_6.R 
qsub -N 7_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_7.R 
qsub -N 8_analysis_BEAM_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_BEAM_cell_downsampling_subset_8.R 

qsub -N 1_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 2_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 3_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 4_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 5_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 6_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 7_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 
qsub -N 8_analysis_cell_downsampling_subset -hold_jid analysis_shalek_data -l mfree=100G Rscript analysis_cell_downsampling_subset_1.R 

## DEG bencmark for HSMM data using bulk RNA-seq as gold standard
qsub -N bulk_deg_benchmark_analysis_HSMM -hold_jid analysis_HSMM_data -l mfree=100G Rscript deg_benchmark_analysis_HSMM_bulk.R 

## function for perform ROC analysis, required for deg_benchmark_analysis_HSMM_bulk  
qsub -N roc_curves -hold_jid deg_benchmark_analysis_HSMM_bulk -l mfree=100G Rscript roc_curves.R 

## Perform branch time point analysis (for shalek data, in particular)
qsub -N branchTimePoint -hold_jid analysis_shalek_data -l mfree=100G Rscript branchTimePoint.R 

## Perform spike-in free recovery using UMI data 
qsub -N umi_normalization -hold_jid umi_normalization -l mfree=100G Rscript umi_normalization.R 

## Compare with other existing software (MAST)
qsub -N cmpr_three_packages -hold_jid deg_benchmark_analysis -l mfree=100G Rscript cmpr_three_packages.R 

## Perform read count downsampling 
qsub -N prepare_lung_downsampling_data -hold_jid prepare_lung_data -l mfree=100G Rscript prepare_lung_downsampling_data.R 

## Making downsampled figure from lung data 
#(WARNING: this script takes about 2 days to finish on a cluster with 64 cores)
qsub -N gen_lung_downsampling_figures -hold_jid prepare_lung_downsampling_data -l mfree=100G Rscript gen_lung_downsampling_figures.R 

## Script to generate other figures not in scripts above 
qsub -N gen_lung_downsampling_figures -hold_jid gen_lung_figures_mc,gen_shalek_figures,gen_supplementary_figure -l mfree=100G Rscript numbers_in_papers.R 

###### The following are scripts developed for the second resubmission #####
##
#########################################################
# New figures added in manuscript and rebuttal figures for second resubmission 
# The following section contains scripts that generates
# figure panels from the paper.
#########################################################
qsub -N third_rebuttal -hold_jid analysis_lung_data -l mfree=100G Rscript third_rebuttal.R 
qsub -N simulator -l mfree=100G Rscript mode_simulation.R 
qsub -N allele_sandberg mfree=100G Rscript allele_sandberg.R 
qsub -N isoform_analysis -l mfree=100G Rscript isoform_analysis.R 
################################################################################################
#done 