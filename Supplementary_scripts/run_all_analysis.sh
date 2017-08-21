## note that this script is for demonstration only, the actual parallelism depends on different systems 
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
the BEAM analysis will be provided: 
tar -zxvf BEAM_analysis_data.tar.gz 
rm BEAM_analysis_data.tar.gz

#make the directories to store the figures generatated in the script 
mkdir -p main_figures supplementary_figures supplementary_data tmp RData

################################################################################################
#each of the following analysis can be run in parallel or separately: 

################################################################
#prepare the data for the lung dataset analysis
qsub -N prepare_lung_data -l mfree=100G Rscript prepare_lung_data.R 

##notee that the following dependent on the prepare_lung_data.R but each of them are independent 
#and can be parallelized

#perform the BEAM analysis for the lung dataset
qsub -N analysis_lung_data -hold_jid prepare_lung_data -l mfree=100G Rscript analysis_lung_data.R 
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
qsub -N gen_lung_figures_mc  -hold_jid analysis_lung_data,spikein_free_algorithm_sampling,deg_benchmark_analysis -l mfree=100G Rscript gen_lung_figures_mc.R
qsub -N gen_shalek_figures -hold_jid analysis_shalek_data -l mfree=100G Rscript gen_shalek_figures.R
qsub -N gen_supplementary_figure -hold_jid analysis_lung_data,spikein_free_algorithm_sampling,deg_benchmark_analysis,analysis_distribution_fitting,analysis_HSMM_data -l mfree=100G Rscript gen_supplementary_figure.R 

#########################################################
# New figures added in manuscript and rebuttal figures for first resubmission 
# All the scripts below doesn't depend on previous figure 
# excepting third_rebuttal
#########################################################

qsub -N deg_benchmark_analysis_HSMM_bulk -hold_jid analysis_HSMM_data  -l mfree=100G deg_benchmark_analysis_HSMM_bulk.R       
qsub -N branchTimePoint -hold_jid gen_shalek_figures,gen_lung_figures_mc -l mfree=100G branchTimePoint.R             
qsub -N umi_normalization -hold_jid analysis_UMI_data -l mfree=100G  umi_normalization.R
qsub -N cmpr_three_packages -hold_jid deg_benchmark_analysis -l mfree=100G  cmpr_three_packages.R         

#########################################################
# New figures added in manuscript and rebuttal figures for second resubmission 
# All the scripts below doesn't depend on previous figure 
# excepting third_rebuttal
#########################################################
qsub -N simulator -l mfree=100G simulator.R       
qsub -N isoform_analysis -l mfree=100G isoform_analysis.R             
qsub -N allele_sandberg -l mfree=100G  allele_sandberg.R
qsub -N cmpr_three_packages -l mfree=100G cmpr_three_packages.R      
qsub -N third_rebuttal -hold_jid analysis_lung_data -l mfree=100G third_rebuttal.R      

#########################################################
# Run the downsampling process in parallel 
#########################################################
qsub -N analysis_cell_downsampling_subset_1 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_1.R       
qsub -N analysis_cell_downsampling_subset_2 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_2.R       
qsub -N analysis_cell_downsampling_subset_3 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_3.R       
qsub -N analysis_cell_downsampling_subset_4 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_4.R       
qsub -N analysis_cell_downsampling_subset_5 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_5.R       
qsub -N analysis_cell_downsampling_subset_6 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_6.R       
qsub -N analysis_cell_downsampling_subset_7 -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_subset_7.R       

qsub -N analysis_BEAM_cell_downsampling_subset_1 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_1.R       
qsub -N analysis_BEAM_cell_downsampling_subset_2 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_2.R       
qsub -N analysis_BEAM_cell_downsampling_subset_3 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_3.R       
qsub -N analysis_BEAM_cell_downsampling_subset_4 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_4.R       
qsub -N analysis_BEAM_cell_downsampling_subset_5 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_5.R       
qsub -N analysis_BEAM_cell_downsampling_subset_6 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_6.R       
qsub -N analysis_BEAM_cell_downsampling_subset_7 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_7.R       
qsub -N analysis_BEAM_cell_downsampling_subset_8 -hold_jid analysis_shalek_data -l mfree=100G analysis_BEAM_cell_downsampling_subset_8.R       

qsub -N analysis_cell_downsampling_A -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_A.R       
qsub -N analysis_cell_downsampling_B -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_B.R       
qsub -N analysis_cell_downsampling_C -hold_jid analysis_shalek_data -l mfree=100G analysis_cell_downsampling_C.R       

qsub -N analysis_cell_downsampling -hold_jid analysis_cell_downsampling_subset_1,analysis_cell_downsampling_subset_2,analysis_cell_downsampling_subset_3,analysis_cell_downsampling_subset_4,analysis_cell_downsampling_subset_5,analysis_cell_downsampling_subset_6,analysis_cell_downsampling_subset_7 -l mfree=100G analysis_cell_downsampling.R       
qsub -N analysis_BEAM_cell_downsampling -hold_jid analysis_cell_downsampling,analysis_BEAM_cell_downsampling_subset_1,analysis_BEAM_cell_downsampling_subset_2,analysis_BEAM_cell_downsampling_subset_3,analysis_BEAM_cell_downsampling_subset_4,analysis_BEAM_cell_downsampling_subset_5,analysis_BEAM_cell_downsampling_subset_6,analysis_BEAM_cell_downsampling_subset_7,analysis_BEAM_cell_downsampling_subset_8 -l mfree=100G analysis_BEAM_cell_downsampling.R       

################################################################################################
#done 