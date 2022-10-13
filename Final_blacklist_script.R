#this script is to generate blacklist bed file for specific genome (here, hg19)

#load library
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)


#set seed
set.seed(7)

#set path for working directory
setwd("/data/kaurh8/Tofilon_SC_data/All_files/merged_lanes_data/WGS_2/deupp1")

#load hg19 genome from aneufinder
bedfile <- system.file("extdata", "hg19_diploid.bam.bed.gz", package="AneuFinderData")

## Make 50kb fixed-width bins
bins <- binReads(bedfile, assembly='hg19', binsize=50000,chromosomes=c(1:22,'X', 'Y'))[[1]]
print(bins)

## Make a plot for visual inspection and get the blacklist
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.95)
p <- plot(bins) + coord_cartesian(ylim=c(0,50))
p <- p + geom_hline(aes(yintercept=lcutoff), color='red')
p <- p + geom_hline(aes(yintercept=ucutoff), color='red')
#print(p)

## Select regions that are above or below the cutoff as blacklist
blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist
blacklist <- reduce(blacklist)
blacklist

## Write blacklist to file
outfile <- paste(getwd(), "hg19_blacklist_50k", sep="/")
#outfile <- paste(getwd(), "hg19_blacklist_50k_Y_exc", sep="/")
outfile 

#save file using export function
exportGRanges(blacklist, filename=outfile , header=FALSE,chromosome.format='NCBI')

#confirm file saved
blacklist_file_path <- paste0(outfile, ".bed.gz" )
blacklist_file_path
