### Install and load packages required for analysis

# Install BiocManager (tool to install and update Bioconductor packages)
#if (!requireNamespace("BiocManager", quietly = TRUE)) # This only needs to be done once
#install.packages("BiocManager")


# Install DESeq2 and apeglm packages (for differential expression analysis)
#BiocManager::install("DESeq2") # This only needs to be done once
#BiocManager::install("apeglm") # This only needs to be done once


# Load DESeq2 and apeglm packages
library(DESeq2)
library(apeglm)

# View DESeq2 user guide
#browseVignettes("DESeq2")


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

#Check whether everything is correct by sampling with head and tail
head(treatment)
head(sampleCondition)
tail(treatment)
tail(sampleCondition)
head(group)


# Create table containing information on all samples
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles2,
                          condition = sampleCondition,
                          treatment = treatment,
                          group = group)

# Build DESeqDataSet using counts files we created
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        #design = ~ family + group)
                                        design = ~ treatment+ group)


############################################ Set rosina as reference level ######################################
ddsHTSeq1$group <- relevel(ddsHTSeq1$group, ref = "Rosina")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq1)) >= 10
ddsHTSeq1 <- ddsHTSeq1[keep,]

# Run differential expression analysis
dds2 <- DESeq(ddsHTSeq1)

# Obtain results for Rosina Trained vs Rosina Control
resultsNames(dds2) # View coefficients
res.lfcShrink <- lfcShrink(dds2, coef="group_Malleti_vs_Rosina", type="apeglm") # Shink 
summary(res.lfcShrink, alpha = 0.05) # Summary stats of results with FDR < 0.05

write.csv(res.lfcShrink, file = "AN_malleti_vs_rosina_AEG.csv", row.names = TRUE)

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
vsd <- vst(dds2, blind=FALSE)
deg.mat2 <- assay(vsd2)
# Calculate Z-scores
deg.mat2 <- t(scale(t(deg.mat2)))

# Plot PCA
plotPCA(vsd, intgroup=c("group", "treatment"))
