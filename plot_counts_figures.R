### Install and load packages required for analysis

# Install BiocManager (tool to install and update Bioconductor packages)
#if (!requireNamespace("BiocManager", quietly = TRUE)) # This only needs to be done once
#install.packages("BiocManager")


# Install DESeq2 and apeglm packages (for differential expression analysis)
#BiocManager::install("DESeq2") # This only needs to be done once
#BiocManager::install("apeglm") # This only needs to be done once


# Load DESeq2 and apeglm packages
library(DelayedArray)
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
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
treatment<-as.factor(treatment)

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
courtship<-as.factor(courtship)
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
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "Rosina_Trained")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Run differential expression analysis
dds <- DESeq(ddsHTSeq)

# Obtain results for Differentially Expressed Genes (DEGs)
resultsNames(dds) # View coefficients
res.lfcShrink <- lfcShrink(dds, coef="condition_Malleti_Trained_vs_Rosina_Trained", type="apeglm") # Shink LFC for genes with small counts


### Visualization that is plotting counts using ggplot2

# Make MA plot
plotMA(res.lfcShrink, alpha = 0.05, ylim=c(-16,16))


# View genes with lowest FDR
head(res.lfcShrink_05)
tail(res.lfcShrink_05)

#Plot in ggplot2 (just change gene=; ggtitle; png(name))
pc<-plotCounts(dds, gene="HMEL005376g1", c("treatment","group"), returnData = TRUE) %>%
  ggplot() + aes(treatment, count) + geom_boxplot(aes(fill=group)) + ylim(0, 4000) + 
  theme_classic()+
  scale_fill_manual("Phenotype", values = c("Malleti"="orange", "Rosina"="darkgrey"))+
  xlab("")+ylab("Normalized counts")+ggtitle("alcohol dehydrogenase (HMEL005376g1)")+
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18),
        legend.position = "none", legend.title = element_text(size = 18), legend.text = element_text(size = 15),
        plot.title = element_text(size = 18, hjust = 0.5))
pc
png("/home/sushant/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_DEG_figures/AN_alcdehyd_trained.png",
    units="in",width=5, height=5, res=1280, bg="transparent")  
pc
dev.off()



