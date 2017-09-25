# Analysis for "Reversed graph embedding resolves complex single-cell trajectories"

## Introduction
This distribution includes all analysis and figure generation performed for our paper, [Reversed graph embedding resolves complex single-cell trajectories](https://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4402.html). In order to run all the code included, you need to install six packages we prepared for this project (See more details in file install_packages.R).

# Jupyter notebook 
A remarkable result of Monocle 2 is its capability to automatically resolve complicate developmental trajectory. In addition to all code we have wrote, we also provided a jupyter notebook to reproduce the developmental trajectory for the Paul (which includes five branch points and six lineages) as well as the Olsson datasets (includes two branch points). See the folder `Jupyter_notebook`

# Citation
If you find Monocle 1 or Monocle 2 helps you to analyze the single cell RNA-seq dataset, please cite the following papers: 
1. Monocle 1: Trapnell, Cole, et al. "The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells." Nature biotechnology 32.4 (2014): 381-386.
APA	
2. Monocle 2: Qiu, Xiaojie, et al. "Single-cell mRNA quantification and differential analysis with Census." Nature methods 14.3 (2017): 309.
3. Monocle 2: Qiu, Xiaojie, et al. "Reversed graph embedding resolves complex single-cell trajectories." Nature methods (2017).

# Data 
Necessary datasets (for example, gene expression matrix, cell or gene annotation file, etc.) can be downloaded from here: 
http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz

# Notes
- Note that users should visit Monocle 2 website (http://cole-trapnell-lab.github.io/monocle-release/) for newest version of Monocle.  
- Since this is an extensive analysis across several different single-cell RNA-seq datasets, we have broken our analysis down into multiple scripts so that users can selectively run parts of the analysis or paralellize the entire analysis as needed.
- To reproduce exactly the result as in the paper, the same packages included here may be necessary (This is because Monocle 2 is still under active development). 
- We suggest the user to look first at the Jupyter notebook which highlights Monocle 2's capabilities to accurately and robustly resolve the developmental trajectory before digging into other complicated analysis performed in this study 

## Requirements
- Memory: We have tested all scripts with 32G of memory (RAM), but the memory usage in reality is much lower for most of the analysis
- Operating system: We performed all our analysis on a ​iMac (3.5 GHz Intel Core i7, 32 GB 1600 MHz DDR3), but Linux and other Unix-based operating systems should also work
- Tools: 
	1. You must have R version 3.3.2 (the one we used) or higher and permissions to install additional packages
	2. One library used during analysis requires rJava, which relies on Java being available (we have tested with Java SE Development Kit 7u17)
- Internet connection: Used to download packages for installation and a tarball with data required for the analysis (too large to provide with this distribution)

Before running any scripts, you must:
- Download the required data and unzip:
 1. wget http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz
 2. tar -zxvf RGE_analysis_data.tar.gz
 3. rm RGE_analysis_data.tar.gz

- Make the relevant destination directories
 1. mkdir Figure/main_figure Figure/supplementary_figures RData tmp

- Install all R packages (this may require user interaction to confirm permission to update existing packages):
 1. Rscript install_packages.R

Once this is complete, individual analyses may be run. The following is a list of all scripts and their dependencies (all dependencies must be run before running a given script):

HSMM Data
- analysis_HSMM_data.r (no dependencies)
- gen_HSMM_figures.r (analysis_HSMM_data.r)
- gen_HSMM_SI_figures.r (analysis_HSMM_data)

Lung Data Analysis
- analysis_lung_benchmark.r (depends on results from prepare_lung_data.r)

Neuron Simulation Analysis
- analysis_neuron_simulation.r (depends on results from prepare_lung_data.r)
- analysis_neuron_simulation_other_methods.r (depends on results from prepare_lung_data.r)
- gen_simulation_dpt_DDRTree_finalized_figures.r (depends on results from prepare_lung_data.r)
- analysis_neuron_simulation_other_methods.r (depends on results from prepare_lung_data.r)
- analysis_complex_tree_structure.r (depends on results from prepare_lung_data.r)

MARS-seq dataset Analysis
- MARSseq_analysis_tutorial.r (no depdencies)
- analysis_mars_seq_data.r (depends on results from MARSseq_analysis_tutorial.r)
- MARSseq_downsampling_empirical_ordering.r (depends on results from MARSseq_analysis_tutorial.r)
- analysis_mars_seq_all_dataset.r (depends on results from MARSseq_analysis_tutorial.r)

Olsson dataset Analysis
- analysis_olsson_WT_data.r (depends on results from prepare_lung_data.r)
- analysis_olsson_KO_data.r (depends on results from analysis_olsson_WT_data.r)

Benchmarking Analysis (algorithm robustness and accuracy; Parameters robustness)
- analysis_DDRTree_parameters.r
- analysis_DDRTree_all_parameters.r (depends on results from analysis_DDRTree_parameters.r)

Shalek Data
- analysis_shalek_data.R (no dependencies)
- gen_shalek_figures.R (depdends on results from analysis_shalek_data.R)

Other:
- calculate_monocle12_slicer_wishbone_dpt_downsampling.r (no dependencies)
- dpFeature_for_all_datasets.r (depends on all the above code)  
- gen_simulation_dpt_DDRTree_finalized_figure.r

## DHS analysis: 
Quake_intersect_dhs.sh and Shalek_intersect_dhs.sh are the shell scripts used to perform the DHS analysis disussed in the supplementary file. For simplicity, we only provided the output gmt files in the data folder for the current distribution.  

New scripts used to generated new figures since the first revision (all depends on script above excepting deg_benchmark_analysis_HSMM_bulk): 
# (folder: revision_1)
- revision_1_HSMM_myo_other_software.R		
- revision_1_comparison_across_software.R		
- revision_1_empirical_downsample_downsampling.R	
- revision_1_simulation_data.R
- revision_1_accuracy_na_simulation_res.R		
- revision_1_dpFeature.R				
- revision_1_new_blood_dataset.R			
- revision_1_test_wishbone_res.R
- revision_1_census.R				
- revision_1_dpt_feature_selection.R		
- revision_1_running_time.R			
- revision_1_traditional_venn_diagram.R

Note that this script will first install all necessary R packages. Once the analysis has run the following folders and files will be populated:
- data: Data to recreate the analysis (too large to provide with this distribution) are downloaded to this folder
- mat_data: Matlab file used in this study, mostly related to the simulation analysis  
- Figure
	- main_figure
	- supplementary_figure
- wishbone: folder storing the wishbone package, code on running wishbone in jupyter notebook data saved for running wishbone 
- csv_data: folder storing csv files used in running the code 
- Nature_hta_paper: folder storing information on the Olsson dataset
- RData: RData files saved for each script so data can be loaded by other scripts in the analysis.
- generate_clustering: Folder storing code/data from Ido Amit on the Paul dataset  
- tmp: Temporary helper pdf files used in figure annontation during the manuscript preparation 
