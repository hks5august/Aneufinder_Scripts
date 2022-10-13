#this script is for the detection of CNV from single cell whole genome data (sorted bam files). with GC correction 
#It will create models for identified CNV which could be used further plot and table generation

#load library
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)

#set seed
set.seed(7)

#set path for working directory
setwd("/data/kaurh8/Tofilon_SC_data/All_files/merged_lanes_data/WGS_2/deupp1")

#load blacklist bed file (for hg19 ) as outfile from the directory. Note: hg19 blacklist bed file must be there in defined folder
outfile <- paste(getwd(), "hg19_blacklist_50k", sep="/")
outfile
blacklist_file_path <- paste0(outfile, ".bed.gz" )
blacklist_file_path

# Run aneufinder
Aneufinder(inputfolder = 'dedupp_files', # provide folder name containing sorted .bam and index files (bai)
           outputfolder = "Complete_output_50k_GC",  #define output folder name to save results
           assembly = 'hg19',  #provide genome version
           binsizes=50000, #define bin size for chromomsomes
           method = c("edivisive"),  #method for CNV detection
           correction.method = 'GC', GC.BSgenome=BSgenome.Hsapiens.UCSC.hg19, #GC correction
           chromosomes = c(1:22,'X','Y'), #select chromosome for which you want detect CNV
           blacklist = blacklist_file_path, #blacklist bed file
           pairedEndReads = TRUE, #define paired or unpaired samples (bam file)
           states = c("zero-inflation", paste0(0:10, "-somy")),  #define ploidy states
           confint = NULL, 
           refine.breakpoints = FALSE, hotspot.bandwidth = NULL, 
           hotspot.pval = 0.05, cluster.plots = TRUE)  #define cluster TRUE if you want to see clusters of samples        
           
