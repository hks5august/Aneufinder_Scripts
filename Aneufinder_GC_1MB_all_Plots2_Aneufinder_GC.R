#load libraries
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)

#set seed
set.seed(7)

#set path
path <- "/Users/kaurh8/Documents/Tofilon_SC_data/Aneufinder_GC_1MB/Complete_output_1MB_GC/MODELS/method-edivisive/"
#path <- "/Users/kaurh8/Documents/Tofilon_SC_data/sorted_bam/Complete_output_5k/MODELS/method-edivisive/"

setwd(path)
getwd()

#load CNV model files
files <- list.files(path=path, pattern=".RData")
files 

#veh files
files.Veh <- grep("Veh", files, ignore.case=TRUE, value=TRUE/FALSE)
files.Veh

#CFI files
files.CFI <- grep("CFI", files, ignore.case=TRUE, value=TRUE/FALSE)
files.CFI


#Cluster by quality
cl <- clusterByQuality(files, measures=c('spikiness','num.segments','entropy','bhattacharyya','sos'))

plot(cl$Mclust, what='classification')


#Save Quality control plot
jpeg(file="ClusterQuality_plot_plot.jpeg", units="in", width=10, height=10, res=350)
plot(cl$Mclust, what='classification')
dev.off()



#Genomewide Heatmap for selected cluster (if needed, but here, ww will create for all files)
cl$classification
#selected.files <- unlist(cl$classification[1:2])
selected.files <- unlist(cl$classification[1:3])
selected.files 
files


Heat_Aneu <- heatmapAneuploidies(files, as.data.frame = FALSE, cluster = F, ylabels = Sample_Labels)
Heat_Aneu

Heat_Aneu <- heatmapAneuploidies(files, as.data.frame = FALSE, cluster = F)
Heat_Aneu

files1 <- rbind(files.Veh , files.CFI )


#aneuploidy heatmap for all samples
Heat_Aneu1 <- heatmapAneuploidies(files1, as.data.frame = FALSE, cluster = F, ylabels = Sample_Labels)
Heat_Aneu1

#aneuploidy heatmap for veh samples

Sample_Labels_Veh <- sapply(strsplit(files.Veh,"_FKDN"), getElement, 1)
Heat_veh <- heatmapAneuploidies(files.Veh, as.data.frame = FALSE, cluster = F, ylabels = Sample_Labels_Veh)
Heat_veh 
jpeg(file="Aneuploidy_Veh_Heatmap.jpeg", units="in", width=25, height=10, res=350)
Heat_veh 
dev.off()


#aneuploidy heatmap for CFI samples
Sample_Labels_CFI <- sapply(strsplit(files.CFI,"_FKDN"), getElement, 1)
Heat_CFI<- heatmapAneuploidies(files.CFI, as.data.frame = FALSE, cluster = F, ylabels = Sample_Labels_CFI)
Heat_CFI 
jpeg(file="Aneuploidy_CFI_Heatmap.jpeg", units="in", width=25, height=10, res=350)
Heat_CFI 
dev.off()

#Genomewide heatmap for veh samples
Heat_veh_GN <- heatmapGenomewide(files.Veh,  cluster = F, ylabels = Sample_Labels_Veh)
Heat_veh_GN
jpeg(file="Genomewide_Veh_Heatmap.jpeg", units="in", width=25, height=10, res=350)
Heat_veh_GN
dev.off()

#Genomewide heatmap for CFI samples
Heat_CFI_GN<- heatmapGenomewide(files.CFI, cluster = F, ylabels = Sample_Labels_CFI)
Heat_CFI_GN 
jpeg(file="Genomewide_CFI_Heatmap.jpeg", units="in", width=25, height=10, res=350)
Heat_CFI_GN
dev.off()


#aneuploidy table for all samples
Heat_Aneu_df <- heatmapAneuploidies(files, as.data.frame = TRUE, ylabels = Sample_Labels)
Heat_Aneu_df
write.table(Heat_Aneu_df,file="Samplewise_Aneu_states.txt", sep='\t',  quote = F,row.names = F)

Heat_Aneu_df1 <- read.table("Samplewise_Aneu_states2.txt", sep='\t', header=T, row.names = 1)
Heat_Aneu_df1
H = sum( table(Heat_Aneu_df1[1:10] ) * 0:(length(table(Heat_Aneu_df[1:10] ))-1) ) / 28


#create Simple samples labels
Sample_Labels <- sapply(strsplit(files,"_FKDN"), getElement, 1)
Sample_Labels


#aneuploidy heatmap for all samples with clusters
Heat_Aneu_clus <- heatmapAneuploidies(files, as.data.frame = FALSE, cluster = TRUE, ylabels = Sample_Labels)
Heat_Aneu_clus

jpeg(file="Aneuploidy_Heatmap_with_clusters.jpeg", units="in", width=15, height=10, res=350)
Heat_Aneu_clus
dev.off()



#Draw genomewide heatmap without clusters 
Heat_G <- heatmapGenomewide(files, cluster = F, ylabels = Sample_Labels)
Heat_G

jpeg(file="Genome_Wide_Heatmap.jpeg", units="in", width=35, height=20, res=350)
Heat_G
dev.off()


jpeg(file="Genome_Wide_Heatmap1.jpeg", units="in", width=30, height=15, res=350)
Heat_G
dev.off()

jpeg(file="Genome_Wide_Heatmap2.jpeg", units="in", width=30, height=12, res=350)
Heat_G
dev.off()



pdf(file="Genome_Wide_Heatmap.pdf", width = 30, height = 15)
Heat_G
dev.off()

#Draw genomewide heatmap with clusters 

Heat_G_with_clus <- heatmapGenomewide(files, cluster = T, ylabels = Sample_Labels)
Heat_G_with_clus 

jpeg(file="GenomeWide_Heatmap_with_Clusters.jpeg", units="in", width=35, height=15, res=350)
Heat_G_with_clus 
dev.off()

pdf(file="GenomeWide_Heatmap_with_Clusters.pdf",  width=35, height=15)
Heat_G_with_clus 
dev.off()

#PCA plots
results <- sapply(files, function(x) mget(load(x)), simplify = TRUE) 
results 
plot_pca(results )

classes1 <- c(rep('CFI', length(files.CFI)), rep('Veh', length(files.Veh)))
classes1 


#Save PCA plot
jpeg(file="PCA_plot.jpeg", units="in", width=10, height=10, res=350)
plot_pca(results, colorBy=classes1, PC1=1, PC2=2)
dev.off()

results <- sapply(files, function(x) mget(load(x)), simplify = TRUE) 
results 
plot_pca(results )

classes1 <- c(rep('Veh', length(files.Veh)), rep('CFI', length(files.CFI)))
classes1 


#Save PCA plot
jpeg(file="PCA_plot_with_Sample_labels.jpeg", units="in", width=10, height=10, res=350)
plot_pca(results, colorBy=Sample_Labels, PC1=1, PC2=2)
#plot_pca(results, colorBy=classes1, PC1=1, PC2=2) 
dev.off()


#Save PCA plot
jpeg(file="PCA_plot.jpeg", units="in", width=10, height=10, res=350)
plot_pca(results, colorBy=classes1, PC1=1, PC2=2)
dev.off()


## Get karyotype measures
k.Veh <- karyotypeMeasures(files.Veh)
k.CFI <- karyotypeMeasures(files.CFI)
## Print the scores in one data.frame
df <- rbind(Veh = k.Veh$genomewide, CFI = k.CFI$genomewide)
head(df)

write.table(df,file="aneu_heter_score_classwise.txt", sep='\t',  quote = F,row.names = T)
df<- read.table("aneu_heter_score_classwise2.txt", header=T, sep='\t')
df
jpeg(file="Aneu_hetero_group_plot_scatterplot_new.jpeg", units="in", width=10, height=10, res=350)
ggplot(df, aes(x=Aneuploidy, y=Heterogeneity, color=Class)) + geom_point()
dev.off()


#compute heterogeneity and aneuploidy per group per chromosome
df1 <- rbind(Veh = k.Veh$per.chromosome, CFI = k.CFI$per.chromosome)
print(df1)
plot(df1)
write.table(df1,file="aneu_heter_score_classwise_per_chr.txt", sep='\t',  quote = F,row.names = T)


H_aneu_plot <- plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_aneu_plot 
H_aneu_plot$data

#save chromosome_wise heterogeneity and aneuploidy of groups
chromosome_wise_df <- as.data.frame(H_aneu_plot$data)
write.table(chromosome_wise_df,file="chromosome_wise_df_class_wise.txt", sep='\t',  quote = F,row.names = T)

#plot aneuploidy vs heterogeneity chromosome wise for veh group
H_aneu_plot_veh_group <- plotHeterogeneity(hmms.list = list(Veh=files.Veh))
H_aneu_plot_veh_group$data

jpeg(file="Veh_group_Aneuplody_heterogeniety_plot1.jpeg", units="in", width=10, height=10, res=350)
#plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_aneu_plot_veh_group
dev.off()

#plot aneuploidy vs heterogeneity chromosome wise for CFI group
H_aneu_plot_CFI_group <- plotHeterogeneity(hmms.list = list(CFI=files.CFI))
H_aneu_plot_CFI_group
H_aneu_plot_CFI_group$data

jpeg(file="CFI_group_Aneuplody_heterogeniety_plot1.jpeg", units="in", width=10, height=10, res=350)
#plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_aneu_plot_CFI_group
dev.off()


H_df_Veh_CFI <-plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_df_Veh_CFI$data

jpeg(file="Aneuplody_heterogeniety_plot.jpeg", units="in", width=10, height=10, res=350)
#plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_aneu_plot 
dev.off()

jpeg(file="Group_wise_Aneuplody_heterogeniety_plot1.jpeg", units="in", width=15, height=10, res=350)
#plotHeterogeneity(hmms.list = list(Veh=files.Veh, CFI=files.CFI))
H_aneu_plot 
dev.off()


pdf(file="Veh_Aneuplody_heterogeniety_plot.pdf",  width=10, height=10)
H_aneu_plot
dev.off()


########## Compute and Plot aneuploidy and heterogeneity per samples per chromosome ############

files.CFI <- grep("CFI", files, ignore.case=TRUE, value=TRUE/FALSE)
files.CFI
CFI1_1_f <- grep("CFI1_1_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_3_f <- grep("CFI1_3_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_4_f <- grep("CFI1_4_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_5_f <- grep("CFI1_5_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_6_f <- grep("CFI1_6_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_7_f <- grep("CFI1_7_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_8_f <- grep("CFI1_8_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_11_f <- grep("CFI1_11_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_13_f <- grep("CFI1_13_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_14_f <- grep("CFI1_14_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_15_f <- grep("CFI1_15_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_18_f <- grep("CFI1_18_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_21_f <- grep("CFI1_21_", files, ignore.case=TRUE, value=TRUE/FALSE)
CFI1_22_f <- grep("CFI1_22_", files, ignore.case=TRUE, value=TRUE/FALSE)

## Get karyotype measures for each sample of CFI group #####
k.CFI1_1 <- karyotypeMeasures(CFI1_1_f)
k.CFI1_3 <- karyotypeMeasures(CFI1_3_f)
k.CFI1_4 <- karyotypeMeasures(CFI1_4_f)
k.CFI1_5 <- karyotypeMeasures(CFI1_5_f)
k.CFI1_6 <- karyotypeMeasures(CFI1_6_f)
k.CFI1_7 <- karyotypeMeasures(CFI1_7_f)
k.CFI1_8 <- karyotypeMeasures(CFI1_8_f)
k.CFI1_11 <- karyotypeMeasures(CFI1_11_f)
k.CFI1_13 <- karyotypeMeasures(CFI1_13_f)
k.CFI1_14 <- karyotypeMeasures(CFI1_14_f)
k.CFI1_15 <- karyotypeMeasures(CFI1_15_f)
k.CFI1_18 <- karyotypeMeasures(CFI1_18_f)
k.CFI1_21 <- karyotypeMeasures(CFI1_21_f)
k.CFI1_22 <- karyotypeMeasures(CFI1_22_f)
head(k.CFI1_1)



## Print the scores in one data.frame
df_CFI <- rbind(CFI_1_sample = k.CFI1_1$genomewide, CFI_3_sample = k.CFI1_3$genomewide, 
                CFI_4_sample = k.CFI1_4$genomewide, CFI_5_sample = k.CFI1_5$genomewide,
                CFI_6_sample = k.CFI1_6$genomewide, CFI_7_sample = k.CFI1_7$genomewide,
                CFI_8_sample = k.CFI1_8$genomewide, CFI_11_sample = k.CFI1_11$genomewide,
                CFI_13_sample = k.CFI1_13$genomewide, CFI_14_sample = k.CFI1_14$genomewide,
                CFI_15_sample = k.CFI1_15$genomewide, CFI_18_sample = k.CFI1_18$genomewide,
                CFI_21_sample = k.CFI1_21$genomewide, CFI_22_sample = k.CFI1_22$genomewide)
print(df_CFI)
write.table(df_CFI,file="df_CFI.txt", sep='\t',  quote = F,row.names = T)

#CFI1_1
#Het_plot_3samples <- plotHeterogeneity(hmms.list = list(CFI1_1=CFI1_1, CFI1_3=CFI1_3, CFI1_5=CFI1_5 ))
Het_aneu_plot_CFIsamples <- plotHeterogeneity(hmms.list = list(CFI1_1=CFI1_1_f,  CFI1_3=CFI1_3_f, 
                                                               CFI1_4=CFI1_4_f, CFI1_5=CFI1_5_f, CFI1_6=CFI1_6_f, 
                                                               CFI1_7=CFI1_7_f,  CFI1_8=CFI1_8_f,CFI1_11=CFI1_11_f,
                                                               CFI1_13=CFI1_13_f , CFI1_14=CFI1_14_f, CFI1_15=CFI1_15_f,
                                                               CFI1_18=CFI1_18_f,  CFI1_21=CFI1_21_f, CFI1_22=CFI1_22_f))

Het_aneu_plot_CFIsamples
chromosome_wise_aneu_CFI <- as.data.frame(Het_aneu_plot_CFIsamples$data)
write.table(chromosome_wise_aneu_CFI,file="chromosome_wise_aneu_CFI.txt", sep='\t',  quote = F,row.names = T)


jpeg(file="Aneup_hetero_plot_CFI_samples.jpeg", units="in", width=20, height=10, res=350)
Het_aneu_plot_CFIsamples
dev.off()

##### Compute & plot aneuplody for each samples of  Veh group

files.Veh <- grep("Veh", files, ignore.case=TRUE, value=TRUE/FALSE)
files.Veh

Veh_1_f <- grep("Veh1_1_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_2_f <- grep("Veh1_2_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_3_f <- grep("Veh1_3_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_4_f <- grep("Veh1_4_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_5_f <- grep("Veh1_5_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_7_f <- grep("Veh1_7_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_8_f <- grep("Veh1_8_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_9_f <- grep("Veh1_9_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_10_f <- grep("Veh1_10_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_12_f <- grep("Veh1_12_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_14_f <- grep("Veh1_14_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_16_f <- grep("Veh1_16_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_18_f <- grep("Veh1_18_", files, ignore.case=TRUE, value=TRUE/FALSE)
Veh_22_f <- grep("Veh1_22_", files, ignore.case=TRUE, value=TRUE/FALSE)
## Get karyotype measures
k.Veh_1 <- karyotypeMeasures(Veh_1_f)
k.Veh_2 <- karyotypeMeasures(Veh_2_f)
k.Veh_3 <- karyotypeMeasures(Veh_3_f)
k.Veh_4 <- karyotypeMeasures(Veh_4_f)
k.Veh_5 <- karyotypeMeasures(Veh_5_f)
k.Veh_7 <- karyotypeMeasures(Veh_7_f)
k.Veh_8 <- karyotypeMeasures(Veh_8_f)
k.Veh_9 <- karyotypeMeasures(Veh_9_f)
k.Veh_10 <- karyotypeMeasures(Veh_10_f)
k.Veh_12 <- karyotypeMeasures(Veh_12_f)
k.Veh_14 <- karyotypeMeasures(Veh_14_f)
k.Veh_16 <- karyotypeMeasures(Veh_16_f)
k.Veh_18 <- karyotypeMeasures(Veh_18_f)
k.Veh_22 <- karyotypeMeasures(Veh_22_f)

head(k.Veh1_1)



## Print the scores in one data.frame
df_veh <- rbind(Veh_1_sample = k.Veh_1$genomewide, Veh_2_sample = k.Veh_2$genomewide, 
                Veh_3_sample = k.Veh_3$genomewide, Veh_4_sample = k.Veh_4$genomewide, 
                Veh_5_sample = k.Veh_5$genomewide, Veh_7_sample = k.Veh_7$genomewide,
                Veh_8_sample = k.Veh_8$genomewide, Veh_9_sample = k.Veh_9$genomewide,
                Veh_10_sample = k.Veh_10$genomewide, Veh_12_sample = k.Veh_12$genomewide,
                Veh_14_sample = k.Veh_14$genomewide, Veh_16_sample = k.Veh_16$genomewide,
                Veh_18_sample = k.Veh_18$genomewide, Veh_22_sample = k.Veh_22$genomewide)
print(df_veh)
write.table(df_veh,file="df_veh.txt", sep='\t',  quote = F,row.names = T)



#Veh1_1
#Het_plot_3samples <- plotHeterogeneity(hmms.list = list(Veh1_1=Veh1_1, Veh1_3=Veh1_3, Veh1_5=Veh1_5 ))
Het_aneu_plot_Veh_samples <- plotHeterogeneity(hmms.list = list(Veh_1=Veh_1_f, Veh_2=Veh_2_f, Veh_3=Veh_3_f, 
                                                                Veh_4=Veh_4_f, Veh_5=Veh_5_f,  Veh_7=Veh_7_f, 
                                                                Veh_8=Veh_8_f,Veh_9=Veh_9_f, Veh_10=Veh_10_f , 
                                                                Veh_12=Veh_12_f, Veh_14=Veh_14_f, Veh_16=Veh_16_f, 
                                                                Veh_18=Veh_18_f, Veh_22=Veh_22_f))

Het_aneu_plot_Veh_samples
head(Het_aneu_plot_Veh_samples$data)
chromosome_wise_aneu_veh <- as.data.frame(Het_aneu_plot_Veh_samples$data)
write.table(chromosome_wise_aneu_veh,file="chromosome_wise_aneu_veh.txt", sep='\t',  quote = F,row.names = T)

jpeg(file="Aneup_hetero_plot_Veh_samples.jpeg", units="in", width=20, height=10, res=350)
Het_aneu_plot_Veh_samples
dev.off()


new_df2<- read.table("aneu_heter_score_classwise_per_chr3.txt", header=T, sep='\t')
head(new_df2)

#p1 <-ggplot(new_df2, aes(x=Aneuploidy, y=Heterogeneity, color=Sample, shape=Class)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
p1 <-ggplot(new_df2, aes(x=Aneuploidy, y=Heterogeneity, color=Class)) + geom_point() + theme(axis.text.x = element_text(angle = 90))

pp1 <- p1 + geom_text(label=new_df2$Chr,  check_overlap = T) # add labels
jpeg(file="Aneup_hetero_group_wise_sactterplot.jpeg", units="in", width=15, height=10, res=350)
pp1
dev.off()


### sample_wise aneuploidy barplot
new_df3 <- read.table("Sample_wise_aneuploidy.txt", header=T, sep='\t')
new_df3


#p <-ggplot(new_df2, aes(x=Sample, y=Aneuploidy, color=Class)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
p <-ggplot(new_df3, aes(x=sample, y=Aneuploidy, color=Class)) + geom_point() + theme(axis.text.x = element_text(angle = 90))

p
p2 <- p  + ylab("Aneuploidy Score") + xlab("Samples")

jpeg(file="Barplot_samples_wise_aneuploidy_plot.jpeg", units="in", width=10, height=10, res=350)
p2
dev.off()

pp <-ggplot(new_df3, aes(x=sample, y=Aneuploidy, size=Aneuploidy,  color=Class)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
pp2 <- pp  + ylab("Aneuploidy Score") + xlab("Samples")

jpeg(file="samples_wise_aneuploidy_plot2.jpeg", units="in", width=10, height=10, res=350)
pp2
dev.off()


df_var<- read.table("aneuploidy_vs_variance", header=T, sep='\t') ## this is computed in excel based on aneuploidy value per sample per chr
head(df_var)

jpeg(file="Aneu_variance_sample_wise_scatterplot_new.jpeg", units="in", width=10, height=10, res=350)
ggplot(df_var, aes(x=Aneuploidy, y=Variance, color=Class)) + geom_point() + geom_text(label=df_var$Sample,  check_overlap = F,  size=3) # add labels
dev.off()



############# Draw plots for Profile, karygram and histogram per samples ############
plot(files[1], type='profile')

#profile
plot(files[1], type='profile')

#Karyogram
plot(files[1], type='karyogram')





##### Karyograms #######

#profile
jpeg(file="Karyogram_CF1_1.jpeg", units="in", width=20, height=10, res=350)
plot(files[1], type='karyogram')
dev.off()

plot(files[2], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_11.jpeg", units="in", width=20, height=10, res=350)
plot(files[2], type='karyogram')
dev.off()

plot(files[3], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_13.jpeg", units="in", width=20, height=10, res=350)
plot(files[3], type='karyogram')
dev.off()

plot(files[4], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[4], type='karyogram')
dev.off()


plot(files[5], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_15.jpeg", units="in", width=20, height=10, res=350)
plot(files[5], type='karyogram')
dev.off()


files[6]
plot(files[6], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[6], type='karyogram')
dev.off()


files[7]
plot(files[7], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_21.jpeg", units="in", width=20, height=10, res=350)
plot(files[7], type='karyogram')
dev.off()


files[8]
plot(files[8], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_22.jpeg", units="in", width=20, height=10, res=350)
plot(files[8], type='karyogram')
dev.off()

files[9]
plot(files[9], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[9], type='karyogram')
dev.off()


files[10]
plot(files[10], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[10], type='karyogram')
dev.off()


files[11]
plot(files[11], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[11], type='karyogram')
dev.off()


files[12]
plot(files[12], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_6.jpeg", units="in", width=20, height=10, res=350)
plot(files[12], type='karyogram')
dev.off()

files[13]
plot(files[13], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[13], type='karyogram')
dev.off()

files[14]
plot(files[14], type='karyogram')
#profile
jpeg(file="Karyogram_CF1_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[14], type='karyogram')
dev.off()

files[15]
plot(files[15], type='karyogram')
#profile
jpeg(file="Karyogram_Veh1.jpeg", units="in", width=20, height=10, res=350)
plot(files[15], type='karyogram')
dev.off()

files[16]
plot(files[16], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_10.jpeg", units="in", width=20, height=10, res=350)
plot(files[16], type='karyogram')
dev.off()


files[17]
plot(files[17], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_12.jpeg", units="in", width=20, height=10, res=350)
plot(files[17], type='karyogram')
dev.off()

files[18]
plot(files[18], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[18], type='karyogram')
dev.off()

files[19]
plot(files[19], type='karyogram')
#profile
jpeg(file="Karyogram_Veh1_16.jpeg", units="in", width=20, height=10, res=350)
plot(files[19], type='karyogram')
dev.off()

files[20]
plot(files[20], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[20], type='karyogram')
dev.off()

files[21]
plot(files[21], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_2.jpeg", units="in", width=20, height=10, res=350)
plot(files[21], type='karyogram')
dev.off()

files[22]
plot(files[22], type='karyogram')
#profile
jpeg(file="Karyogram_Veh22.jpeg", units="in", width=20, height=10, res=350)
plot(files[22], type='karyogram')
dev.off()

files[23]
plot(files[23], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[23], type='karyogram')
dev.off()

files[24]
plot(files[24], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[24], type='karyogram')
dev.off()


files[25]
plot(files[25], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[25], type='karyogram')
dev.off()

files[26]
plot(files[26], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[26], type='karyogram')
dev.off()


files[27]
plot(files[27], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[27], type='karyogram')
dev.off()


files[28]
plot(files[28], type='karyogram')
#profile
jpeg(file="Karyogram_Veh_9.jpeg", units="in", width=20, height=10, res=350)
plot(files[28], type='karyogram')
dev.off()







########## profiles plots

#profile
jpeg(file="profile_CF1_1.jpeg", units="in", width=20, height=10, res=350)
plot(files[1], type='profile')
dev.off()

plot(files[2], type='profile')
#profile
jpeg(file="profile_CF1_11.jpeg", units="in", width=20, height=10, res=350)
plot(files[2], type='profile')
dev.off()

plot(files[3], type='profile')
#profile
jpeg(file="profile_CF1_13.jpeg", units="in", width=20, height=10, res=350)
plot(files[3], type='profile')
dev.off()

plot(files[4], type='profile')
#profile
jpeg(file="profile_CF1_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[4], type='profile')
dev.off()


plot(files[5], type='profile')
#profile
jpeg(file="profile_CF1_15.jpeg", units="in", width=20, height=10, res=350)
plot(files[5], type='profile')
dev.off()


files[6]
plot(files[6], type='profile')
#profile
jpeg(file="profile_CF1_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[6], type='profile')
dev.off()


files[7]
plot(files[7], type='profile')
#profile
jpeg(file="profile_CF1_21.jpeg", units="in", width=20, height=10, res=350)
plot(files[7], type='profile')
dev.off()


files[8]
plot(files[8], type='profile')
#profile
jpeg(file="profile_CF1_22.jpeg", units="in", width=20, height=10, res=350)
plot(files[8], type='profile')
dev.off()

files[9]
plot(files[9], type='profile')
#profile
jpeg(file="profile_CF1_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[9], type='profile')
dev.off()


files[10]
plot(files[10], type='profile')
#profile
jpeg(file="profile_CF1_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[10], type='profile')
dev.off()


files[11]
plot(files[11], type='profile')
#profile
jpeg(file="profile_CF1_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[11], type='profile')
dev.off()


files[12]
plot(files[12], type='profile')
#profile
jpeg(file="profile_CF1_6.jpeg", units="in", width=20, height=10, res=350)
plot(files[12], type='profile')
dev.off()

files[13]
plot(files[13], type='profile')
#profile
jpeg(file="profile_CF1_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[13], type='profile')
dev.off()

files[14]
plot(files[14], type='profile')
#profile
jpeg(file="profile_CF1_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[14], type='profile')
dev.off()

files[15]
plot(files[15], type='profile')
#profile
jpeg(file="profile_Veh1.jpeg", units="in", width=20, height=10, res=350)
plot(files[15], type='profile')
dev.off()

files[16]
plot(files[16], type='profile')
#profile
jpeg(file="profile_Veh_10.jpeg", units="in", width=20, height=10, res=350)
plot(files[16], type='profile')
dev.off()


files[17]
plot(files[17], type='profile')
#profile
jpeg(file="profile_Veh_12.jpeg", units="in", width=20, height=10, res=350)
plot(files[17], type='profile')
dev.off()

files[18]
plot(files[18], type='profile')
#profile
jpeg(file="profile_Veh_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[18], type='profile')
dev.off()

files[19]
plot(files[19], type='profile')
#profile
jpeg(file="profile_Veh1_16.jpeg", units="in", width=20, height=10, res=350)
plot(files[19], type='profile')
dev.off()

files[20]
plot(files[20], type='profile')
#profile
jpeg(file="profile_Veh_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[20], type='profile')
dev.off()

files[21]
plot(files[21], type='profile')
#profile
jpeg(file="profile_Veh_2.jpeg", units="in", width=20, height=10, res=350)
plot(files[21], type='profile')
dev.off()

files[22]
plot(files[22], type='profile')
#profile
jpeg(file="profile_Veh22.jpeg", units="in", width=20, height=10, res=350)
plot(files[22], type='profile')
dev.off()

files[23]
plot(files[23], type='profile')
#profile
jpeg(file="profile_Veh_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[23], type='profile')
dev.off()

files[24]
plot(files[24], type='profile')
#profile
jpeg(file="profile_Veh_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[24], type='profile')
dev.off()


files[25]
plot(files[25], type='profile')
#profile
jpeg(file="profile_Veh_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[25], type='profile')
dev.off()

files[26]
plot(files[26], type='profile')
#profile
jpeg(file="profile_Veh_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[26], type='profile')
dev.off()


files[27]
plot(files[27], type='profile')
#profile
jpeg(file="profile_Veh_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[27], type='profile')
dev.off()


files[28]
plot(files[28], type='profile')
#profile
jpeg(file="profile_Veh_9.jpeg", units="in", width=20, height=10, res=350)
plot(files[28], type='profile')
dev.off()




######## Histograms

#Karygram
jpeg(file="histogram_CF1_1.jpeg", units="in", width=20, height=10, res=350)
plot(files[1], type='histogram')
dev.off()

plot(files[2], type='histogram')
#Karygram
jpeg(file="histogram_CF1_11.jpeg", units="in", width=20, height=10, res=350)
plot(files[2], type='histogram')
dev.off()

plot(files[3], type='histogram')
#Karygram
jpeg(file="histogram_CF1_13.jpeg", units="in", width=20, height=10, res=350)
plot(files[3], type='histogram')
dev.off()

plot(files[4], type='histogram')
#Karygram
jpeg(file="histogram_CF1_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[4], type='histogram')
dev.off()


plot(files[5], type='histogram')
#Karygram
jpeg(file="histogram_CF1_15.jpeg", units="in", width=20, height=10, res=350)
plot(files[5], type='histogram')
dev.off()


files[6]
plot(files[6], type='histogram')
#Karygram
jpeg(file="histogram_CF1_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[6], type='histogram')
dev.off()


files[7]
plot(files[7], type='histogram')
#Karygram
jpeg(file="histogram_CF1_21.jpeg", units="in", width=20, height=10, res=350)
plot(files[7], type='histogram')
dev.off()


files[8]
plot(files[8], type='histogram')
#Karygram
jpeg(file="histogram_CF1_22.jpeg", units="in", width=20, height=10, res=350)
plot(files[8], type='histogram')
dev.off()

files[9]
plot(files[9], type='histogram')
#Karygram
jpeg(file="histogram_CF1_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[9], type='histogram')
dev.off()


files[10]
plot(files[10], type='histogram')
#Karygram
jpeg(file="histogram_CF1_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[10], type='histogram')
dev.off()


files[11]
plot(files[11], type='histogram')
#Karygram
jpeg(file="histogram_CF1_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[11], type='histogram')
dev.off()


files[12]
plot(files[12], type='histogram')
#Karygram
jpeg(file="histogram_CF1_6.jpeg", units="in", width=20, height=10, res=350)
plot(files[12], type='histogram')
dev.off()

files[13]
plot(files[13], type='histogram')
#Karygram
jpeg(file="histogram_CF1_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[13], type='histogram')
dev.off()

files[14]
plot(files[14], type='histogram')
#Karygram
jpeg(file="histogram_CF1_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[14], type='histogram')
dev.off()

files[15]
plot(files[15], type='histogram')
#Karygram
jpeg(file="histogram_Veh1.jpeg", units="in", width=20, height=10, res=350)
plot(files[15], type='histogram')
dev.off()

files[16]
plot(files[16], type='histogram')
#Karygram
jpeg(file="histogram_Veh_10.jpeg", units="in", width=20, height=10, res=350)
plot(files[16], type='histogram')
dev.off()


files[17]
plot(files[17], type='histogram')
#Karygram
jpeg(file="histogram_Veh_12.jpeg", units="in", width=20, height=10, res=350)
plot(files[17], type='histogram')
dev.off()

files[18]
plot(files[18], type='histogram')
#Karygram
jpeg(file="histogram_Veh_14.jpeg", units="in", width=20, height=10, res=350)
plot(files[18], type='histogram')
dev.off()

files[19]
plot(files[19], type='histogram')
#Karygram
jpeg(file="histogram_Veh1_16.jpeg", units="in", width=20, height=10, res=350)
plot(files[19], type='histogram')
dev.off()

files[20]
plot(files[20], type='histogram')
#Karygram
jpeg(file="histogram_Veh_18.jpeg", units="in", width=20, height=10, res=350)
plot(files[20], type='histogram')
dev.off()

files[21]
plot(files[21], type='histogram')
#Karygram
jpeg(file="histogram_Veh_2.jpeg", units="in", width=20, height=10, res=350)
plot(files[21], type='histogram')
dev.off()

files[22]
plot(files[22], type='histogram')
#Karygram
jpeg(file="histogram_Veh22.jpeg", units="in", width=20, height=10, res=350)
plot(files[22], type='histogram')
dev.off()

files[23]
plot(files[23], type='histogram')
#Karygram
jpeg(file="histogram_Veh_3.jpeg", units="in", width=20, height=10, res=350)
plot(files[23], type='histogram')
dev.off()

files[24]
plot(files[24], type='histogram')
#Karygram
jpeg(file="histogram_Veh_4.jpeg", units="in", width=20, height=10, res=350)
plot(files[24], type='histogram')
dev.off()


files[25]
plot(files[25], type='histogram')
#Karygram
jpeg(file="histogram_Veh_5.jpeg", units="in", width=20, height=10, res=350)
plot(files[25], type='histogram')
dev.off()

files[26]
plot(files[26], type='histogram')
#Karygram
jpeg(file="histogram_Veh_7.jpeg", units="in", width=20, height=10, res=350)
plot(files[26], type='histogram')
dev.off()


files[27]
plot(files[27], type='histogram')
#Karygram
jpeg(file="histogram_Veh_8.jpeg", units="in", width=20, height=10, res=350)
plot(files[27], type='histogram')
dev.off()


files[28]
plot(files[28], type='histogram')
#Karygram
jpeg(file="histogram_Veh_9.jpeg", units="in", width=20, height=10, res=350)
plot(files[28], type='histogram')
dev.off()



