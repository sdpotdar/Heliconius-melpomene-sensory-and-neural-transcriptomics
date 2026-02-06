### Install and load packages required for analysis

# Install BiocManager (tool to install and update Bioconductor packages)
#if (!requireNamespace("BiocManager", quietly = TRUE)) # This only needs to be done once
#install.packages("BiocManager")


# Install DESeq2 and apeglm packages (for differential expression analysis)
BiocManager::install("DESeq2") # This only needs to be done once
BiocManager::install("apeglm") # This only needs to be done once


# Load DESeq2 and apeglm packages
library(DelayedArray)
library(DESeq2)
library(apeglm)

# View DESeq2 user guide
browseVignettes("DESeq2")


################################################################################################################################################


### Convert STAR GeneCounts files to match format of HTseq files (i.e., extract stranded=reverse counts); this only needs to be done once

# Make variable containing path to directory with counts files
directory <- "~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads"

# Make a list of all the counts files in that directory
sampleFiles <- grep("ReadsPerGene.out.tab",list.files(directory),value=TRUE)

# Make new counts files to be used for differential expression analysis
setwd(directory)
for(file in sampleFiles){
  a <- read.table(paste(directory, "/", file, sep=""), header = F) # Read in each file in directory
  a <- a[,c(1,4)] # Extract column 1 (gene name) and column 4 (reverse counts)
  a <- a[c(-1:-4),] # Remove the first 4 rows (stats)
  write.table(a, paste(directory, "/", file, "_REVERSE-ONLY.txt", sep = ""),  sep="\t", col.names = F, row.names = F, quote=F) # Write new file
}

################################################################################################################################################


### DIFFERENTIAL EXPRESSION ANALYSIS

# Make list of new counts files (those we created above) 
sampleFiles2 <- grep("REVERSE-ONLY.txt",list.files(directory),value=TRUE)

# Collect sample names
sampleName <- sub("_Reads.*","", sampleFiles2)

# Collect sample condition for each sample (i.e., Rosina or Malleti and Control or Trained)
sampleCondition <- sub("_M_","_Malleti_", sampleFiles2)
sampleCondition <- sub("_R_","_Rosina_", sampleCondition)
sampleCondition <- sub("_C_","_Control_", sampleCondition)
sampleCondition <- sub("_T_","_Trained_", sampleCondition)
sampleCondition <- sub("_ReadsPerGene.*","", sampleCondition)
sampleCondition <- sub("SP-*0.","", sampleCondition)
sampleCondition <- sub("SP-*1.","", sampleCondition)
sampleCondition <- sub("SP-*2.","", sampleCondition)
sampleCondition <- sub("SP-*3.","", sampleCondition)
sampleCondition <- sub("SP-*4.","", sampleCondition)
sampleCondition <- sub("SP-*5.","", sampleCondition)
sampleCondition <- sub("SP-*6.","", sampleCondition)
sampleCondition <- sub("SP-*7.","", sampleCondition)
sampleCondition <- sub("SP-*8.","", sampleCondition)
sampleCondition <- sub("_Malleti","Malleti", sampleCondition)
sampleCondition <- sub("_Rosina","Rosina", sampleCondition)

sampleCondition <- sub("_AN","", sampleCondition)


#Collect treatment condition for each treatment i.e Control and Trained
treatment<-sub(".*_Control","Control", sampleCondition)
treatment<-sub(".*_Trained","Trained", treatment)

#Collect group condition for each sample i.e whether its rosina or malleti
group<-sub("Malleti_Control", "Malleti", sampleCondition)
group<-sub("Malleti_Trained", "Malleti", group)
group<-sub("Rosina_Control", "Rosina", group)
group<-sub("Rosina_Trained", "Rosina", group)
group<-as.factor(group)

#Check whether everything is correct by sampling with head and tail
head(treatment)
head(sampleCondition)
tail(treatment)
tail(sampleCondition)
head(group)
tail(group)

##Create vector with courtship information in the same order as the samples appear 
courtship<-c("NA", "NA", "Yes", "NA", "No", "NA", "NA", "Yes", "NA", "NA",
             "No", "NA", "NA", "No", "NA", "No", "Yes", "NA", "NA", "Yes",
             "No", "NA", "NA", "NA", "Yes", "Yes", "NA", "NA", "Yes", "No",
             "NA", "NA", "Yes", "NA", "No", "Yes", "NA")
             
# Create table containing information on all samples
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles2,
                          condition = sampleCondition,
                          group = group,
                          treatment = treatment,
                          courtship = courtship)


# Build DESeqDataSet using counts files we created
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

# Set rosina control as reference level
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "Rosina_Control")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Run differential expression analysis
dds <- DESeq(ddsHTSeq)

# Obtain results for Differentially Expressed Genes (DEGs)
resultsNames(dds) # View coefficients
res.lfcShrink <- lfcShrink(dds, coef="condition_Malleti_Control_vs_Rosina_Control", type="apeglm") # Shink LFC for genes with small counts
summary(res.lfcShrink, alpha = 0.05) # Summary stats of results with FDR < 0.05
sum(res.lfcShrink$padj < 0.05, na.rm=TRUE) # Total number of DEGs (FDR < 0.05)
res.lfcShrink.ordered <- res.lfcShrink[order(res.lfcShrink$padj),] # All results ordered by padj (i.e., FDR)
res.lfcShrink_05 <- subset(res.lfcShrink.ordered, padj <= 0.05) # Make subset of all significant DEGs (FDR < 0.05)
res.lfcShrink_05_lfc <- res.lfcShrink_05[order(-abs(res.lfcShrink_05$log2FoldChange)),] # Order significant DEGs by abs(LFC)
res.lfcShrink_05_up <- subset(res.lfcShrink_05, log2FoldChange > 0) # Make subset of all upregulated significant DEGs (FDR < 0.05)
res.lfcShrink_05_down <- subset(res.lfcShrink_05, log2FoldChange < 0) # Make subset of all downregulated significant DEGs (FDR < 0.05)
write.csv(res.lfcShrink_05, file = "AN_malleticontrol_vs_rosinacontrol_DEG.csv", row.names = TRUE) #CSV file of differentially expressed genes



#######################################################################################################################################




### Visualization

# Make MA plot
plotMA(res.lfcShrink, alpha = 0.05, ylim=c(-16,16))


# View genes with lowest FDR
head(res.lfcShrink_05)
tail(res.lfcShrink_05)

# Plot counts of gene with lowest FDR
plotCounts(dds, "HMEL016620g1", intgroup = "treatment", normalized = T)
# Perform variance stabilizing transformation for PCA 
vsd <- vst(dds, blind=FALSE)
deg.mat <- assay(vsd)
# Calculate Z-scores
deg.mat <- t(scale(t(deg.mat)))

# Plot PCA
plotPCA(vsd, intgroup=c("group", "treatment"))

#################################################### Heatmaps #########################################
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendsort)
library(tibble)
library(cowplot)

# Rosina trained vs Rosina Control
#DEGs

sig_genes_AN <- read.delim("rosinatrained_vs_rosinacontrol_degs.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat <- as.data.frame(deg.mat[grepl(paste(sig_genes_AN, collapse='|'),
                                       row.names(deg.mat)), ]) %>%
  dplyr::select(matches(c("*._R_C_.*", "*._R_T_.*")))

# Specify column annotation data for heatmap
annot2 <- colData(dds)
annot2 <- as.data.frame(annot2) %>% dplyr::select(treatment, courtship)
colnames(annot2)[1] <- "Treatment"
colnames(annot2)[2] <- "Courtship"

# Create list with colors for each annotation
annotation_colors2 <- list(Treatment = c(Trained="#994F00", Control="grey"),
                           Courtship = c(Yes="darkgreen", No="#D41159", "NA"="#FEFE62")) # col colours

# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat)/paletteLength, max(deg.mat), length.out=floor(paletteLength/2)))
  
# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

rosinatrained_control_heatmap  <- pheatmap(deg.mat,
                                   main = "Rosina Antenna Trained vs Control DEGS",
                                   color = my_palette,
                                   border_color = NA,
                                   annotation_colors = annotation_colors2,
                                   cluster_rows=TRUE, 
                                   cluster_cols=TRUE, 
                                   annotation_col=annot2,
                                   annotation_names_col = FALSE, # deletes row/ col names
                                   annotation_names_row = FALSE,
                                   clustering_callback = callback,
                                   show_colnames = FALSE, # deletes individual row/ col cluster names
                                   show_rownames = FALSE,
                                   breaks = myBreaks, 
                                   fontsize = 8)# affects scale of exported heatmap png,

png(file= "Rosina_AN_trainedvscontrol_DEGs.png", height = 1600, width = 1920, res=300)
rosinatrained_control_heatmap
dev.off()


# Malleti Control vs Rosina Control
#DEGs

sig_genes2_AN <- read.delim("malleticontrol_vs_rosinacontrol_degs.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat <- as.data.frame(deg.mat[grepl(paste(sig_genes2_AN, collapse='|'),
                                       row.names(deg.mat)), ]) %>%
  dplyr::select(matches(c("*._M_C_AN", "*._R_C_AN")))

# Specify column annotation data for heatmap
annot <- colData(dds)
annot <- as.data.frame(annot) %>% dplyr::select(group)
colnames(annot)[1] <- "Group"

# Create list with colors for each annotation
annotation_colors3 <- list(Group = c(Malleti="orange", Rosina="black")) # col colours

# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat)/paletteLength, max(deg.mat), length.out=floor(paletteLength/2)))

# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

control_malleti_rosina_heatmap  <- pheatmap(deg.mat,
                                           main = "Control Antenna Malleti vs Rosina DEGS",
                                           color = my_palette,
                                           border_color = NA,
                                           annotation_colors = annotation_colors3,
                                           cluster_rows=TRUE, 
                                           cluster_cols=TRUE, 
                                           annotation_col=annot,
                                           annotation_names_col = FALSE, # deletes row/ col names
                                           annotation_names_row = FALSE,
                                           clustering_callback = callback,
                                           show_colnames = FALSE, # deletes individual row/ col cluster names
                                           show_rownames = FALSE,
                                           breaks = myBreaks, 
                                           fontsize = 8)# affects scale of exported heatmap png,

png(file= "Control_AN_malletivsrosina_DEGs.png", height = 1600, width = 1920, res=300)
control_malleti_rosina_heatmap
dev.off()

################## Setting malleti control as reference #############################

# Build DESeqDataSet using counts files we created
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

# Set malleti control as reference level
ddsHTSeq2$condition <- relevel(ddsHTSeq$condition, ref = "Malleti_Control")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq2)) >= 10
ddsHTSeq2 <- ddsHTSeq2[keep,]

# Run differential expression analysis
dds2 <- DESeq(ddsHTSeq2)

# Obtain results
resultsNames(dds2) # View coefficients
res.lfcShrink <- lfcShrink(dds2, coef="condition_Malleti_Trained_vs_Malleti_Control", type="apeglm") # Shink LFC for genes with small counts
summary(res.lfcShrink, alpha = 0.05) # Summary stats of results with FDR < 0.05
sum(res.lfcShrink$padj < 0.05, na.rm=TRUE) # Total number of DEGs (FDR < 0.05)
res.lfcShrink.ordered <- res.lfcShrink[order(res.lfcShrink$padj),] # All results ordered by padj (i.e., FDR)
res.lfcShrink_05 <- subset(res.lfcShrink.ordered, padj <= 0.05) # Make subset of all significant DEGs (FDR < 0.05)
res.lfcShrink_05_lfc <- res.lfcShrink_05[order(-abs(res.lfcShrink_05$log2FoldChange)),] # Order significant DEGs by abs(LFC)
res.lfcShrink_05_up <- subset(res.lfcShrink_05, log2FoldChange > 0) # Make subset of all upregulated significant DEGs (FDR < 0.05)
res.lfcShrink_05_down <- subset(res.lfcShrink_05, log2FoldChange < 0) # Make subset of all downregulated significant DEGs (FDR < 0.05)

write.csv(res.lfcShrink_05, file = "AN_rosinatrained_vs_rosinacontrol_DEG.csv", row.names = TRUE)


#######################################################################################################################################




### Visualization

# Make MA plot
plotMA(res.lfcShrink, alpha = 0.05, ylim=c(-16,16))


# View genes with lowest FDR
head(res.lfcShrink_05)
tail(res.lfcShrink_05)

# Plot counts of gene with lowest FDR
plotCounts(dds, "HMEL016620g1", intgroup = "group", normalized = T)

# Perform variance stabilizing transformation for PCA 
vsd2 <- vst(dds2, blind=FALSE)
deg.mat2 <- assay(vsd2)
# Calculate Z-scores
deg.mat2 <- t(scale(t(deg.mat2)))

# Plot PCA
plotPCA(vsd, intgroup=c("group", "treatment"))

#################################################### Heatmaps #########################################
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendsort)
library(tibble)
library(cowplot)

# Rosina trained vs Rosina Control
#DEGs

sig_genes3_AN <- read.delim("malletitrained_vs_malleticontrol_degs.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat2 <- as.data.frame(deg.mat2[grepl(paste(sig_genes3_AN, collapse='|'),
                                       row.names(deg.mat2)), ]) %>%
  dplyr::select(matches(c("*._M_C_.*", "*._M_T_.*")))

# Specify column annotation data for heatmap
annot4 <- colData(dds2)
annot4 <- as.data.frame(annot4) %>% dplyr::select(treatment, courtship)
colnames(annot4)[1] <- "Treatment"
colnames(annot4)[2] <- "Courtship"
# Create list with colors for each annotation
annotation_colors4 <- list(Treatment = c(Trained="#994F00", Control="grey"),
                           Courtship=c(Yes="darkgreen", No="#D41159", "NA"="#FEFE62")) # col colours

# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat2)/paletteLength, max(deg.mat2), length.out=floor(paletteLength/2)))

# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

malletitrained_control_heatmap  <- pheatmap(deg.mat2,
                                           main = "Malleti Antenna Trained vs Control DEGS",
                                           color = my_palette,
                                           border_color = NA,
                                           annotation_colors = annotation_colors4,
                                           cluster_rows=TRUE, 
                                           cluster_cols=TRUE, 
                                           annotation_col=annot4,
                                           annotation_names_col = FALSE, # deletes row/ col names
                                           annotation_names_row = FALSE,
                                           clustering_callback = callback,
                                           show_colnames = FALSE, # deletes individual row/ col cluster names
                                           show_rownames = FALSE,
                                           breaks = myBreaks, 
                                           fontsize = 8)# affects scale of exported heatmap png,

png(file= "Malleti_AN_trainedvscontrol_DEGs.png", height = 1600, width = 1920, res=300)
malletitrained_control_heatmap
dev.off()


################## Setting rosina trained as reference #############################

# Build DESeqDataSet using counts files we created
ddsHTSeq3 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~ condition)

# Set malleti control as reference level
ddsHTSeq3$condition <- relevel(ddsHTSeq$condition, ref = "Rosina_Trained")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq3)) >= 10
ddsHTSeq3 <- ddsHTSeq3[keep,]

# Run differential expression analysis
dds3 <- DESeq(ddsHTSeq3)

# Obtain results
resultsNames(dds3) # View coefficients
res.lfcShrink <- lfcShrink(dds3, coef="condition_Malleti_Trained_vs_Rosina_Trained", type="apeglm") # Shink LFC for genes with small counts
summary(res.lfcShrink, alpha = 0.05) # Summary stats of results with FDR < 0.05
sum(res.lfcShrink$padj < 0.05, na.rm=TRUE) # Total number of DEGs (FDR < 0.05)
res.lfcShrink.ordered <- res.lfcShrink[order(res.lfcShrink$padj),] # All results ordered by padj (i.e., FDR)
res.lfcShrink_05_03 <- subset(res.lfcShrink.ordered, padj <= 0.05) # Make subset of all significant DEGs (FDR < 0.05)
res.lfcShrink_05_lfc <- res.lfcShrink_05[order(-abs(res.lfcShrink_05$log2FoldChange)),] # Order significant DEGs by abs(LFC)
res.lfcShrink_05_up <- subset(res.lfcShrink_05, log2FoldChange > 0) # Make subset of all upregulated significant DEGs (FDR < 0.05)
res.lfcShrink_05_down <- subset(res.lfcShrink_05, log2FoldChange < 0) # Make subset of all downregulated significant DEGs (FDR < 0.05)

write.csv(res.lfcShrink_05_03, file = "AN_malletitrained_vs_rosinatrained_DEG.csv", row.names = TRUE)


#######################################################################################################################################




### Visualization

# Make MA plot
plotMA(res.lfcShrink, alpha = 0.05, ylim=c(-16,16))


# View genes with lowest FDR
head(res.lfcShrink_05)
tail(res.lfcShrink_05)

# Plot counts of gene with lowest FDR
plotCounts(dds, "HMEL016620g1", intgroup = "group", normalized = T)

# Perform variance stabilizing transformation for PCA 
vsd3 <- vst(dds3, blind=FALSE)
deg.mat3 <- assay(vsd3)
# Calculate Z-scores
deg.mat3 <- t(scale(t(deg.mat3)))

# Plot PCA
plotPCA(vsd, intgroup=c("group", "treatment"))

#################################################### Heatmaps #########################################
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendsort)
library(tibble)
library(cowplot)

# MAlleti trained vs Rosina Trained
#DEGs

sig_genes4_AN <- read.delim("malletitrained_vs_rosinatrained_degs.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat3 <- as.data.frame(deg.mat3[grepl(paste(sig_genes4_AN, collapse='|'),
                                         row.names(deg.mat3)), ]) %>%
  dplyr::select(matches(c("*._M_T_AN", "*._R_T_AN")))

# Specify column annotation data for heatmap
annot5 <- colData(dds3)
annot5 <- as.data.frame(annot5) %>% dplyr::select(group, courtship)
colnames(annot5)[1] <- "Group"
colnames(annot5)[2] <- "Courtship"
# Create list with colors for each annotation
annotation_colors5 <- list(Group = c(Malleti="orange", Rosina="black"),
                           Courtship = c(Yes="darkgreen", No="#D41159")) # col colours

# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat3)/paletteLength, max(deg.mat3), length.out=floor(paletteLength/2)))

# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

trained_malletirosina_heatmap  <- pheatmap(deg.mat3,
                                            main = "Trained Antenna Malleti vs Rosina DEGS",
                                            color = my_palette,
                                            border_color = NA,
                                            annotation_colors = annotation_colors5,
                                            cluster_rows=TRUE, 
                                            cluster_cols=TRUE, 
                                            annotation_col=annot5,
                                            annotation_names_col = FALSE, # deletes row/ col names
                                            annotation_names_row = FALSE,
                                            clustering_callback = callback,
                                            show_colnames = FALSE, # deletes individual row/ col cluster names
                                            show_rownames = FALSE,
                                            breaks = myBreaks, 
                                            fontsize = 8)# affects scale of exported heatmap png,
trained_malletirosina_heatmap
png(file= "Trained_AN_malletivsrosina_DEGs.png", height = 1600, width = 1920, res=300)
trained_malletirosina_heatmap
dev.off()
