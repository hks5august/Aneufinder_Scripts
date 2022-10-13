#this script is to generate blacklist.bed file for genome (hg19), this bed file would be used for aneufinder

Rscript Final_blacklist_script.R

#this script will generate CNV from single cell whole genome (sorted bam and index file) data. This will generate outputs as MODELs and Plots

Rscript Final_aneufinder_with_backlist_GC_corr_script.R

#This script will generate plots and tables based on the CNV models generated from aneufinder

Rscript Aneufinder_Final_Plots_from_models.R
