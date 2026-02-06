setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA")
### Filter and perform variance stabilizing transformation on count data ###

library(DelayedArray)
library(DESeq2)
library(apeglm)

### Convert STAR GeneCounts files to match format of HTseq files (i.e., extract stranded=reverse counts); this only needs to be done once

# Make variable containing path to directory with counts files
directory <- "~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads"

# Make a list of all the counts files in that directory
sampleFiles <- grep("ReadsPerGene.out.tab",list.files(directory),value=TRUE)

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
#treatment<-as.factor(treatment)
head(treatment)
#Collect group condition for each sample i.e whether its rosina or malleti
group<-sub("Malleti_Control", "Malleti", sampleCondition)
group<-sub("Malleti_Trained", "Malleti", group)
group<-sub("Rosina_Control", "Rosina", group)
group<-sub("Rosina_Trained", "Rosina", group)
#group<-as.factor(group)

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
write.csv(sampleTable, "~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA/AN_sample_traits_from_SampleTable.csv")

# Build DESeqDataSet using counts files we created
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ group + treatment + group:treatment)

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]

# Set rosina control as reference level
ddsHTSeq$group <- relevel(ddsHTSeq$group, ref = "Rosina")
ddsHTSeq$treatment <- relevel(ddsHTSeq$treatment, ref = "Control")
# Run differential expression analysis
dds <- DESeq(ddsHTSeq)

# Run variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

write.csv(assay(vsd), "~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA/AN_variance_stabilizing_transformation.csv")

#=====================================================================================#
# #=====================================================================================#
# 
# ### BELOW IS CODE THAT WAS USED TO MAKE PAIRWISE BINARY INDICATOR VARIABLES FOR ANTENNAE (I.E., MAKE VARIABLES THAT CONTRAST EACH GROUP AGAINST EACH OTHER GROUP)
# ### See https://peterlangfelder.com/2018/11/25/working-with-categorical-variables/
# 
library(WGCNA)
library(tidyr)
library(dplyr)
# 
# 
# # Define a categorical variable with 2 levels
antenna_samples <- read.csv("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA/AN_sample_traits.csv") %>% 
   select(sampleName)

antenna_samples2 <- antenna_samples %>%
 {gsub("_.*", "", .$sampleName)}
# 
# # Binarize it into pairwise indicators
out = binarizeCategoricalVariable(antenna_samples2,
                                    includePairwise = TRUE,
                                    includeLevelVsAll = FALSE)
# 
# # Print the variable and the indicators
datTraits <- data.frame(antenna_samples, out)
# 
# # Write to file
write.csv(datTraits, "ANTENNA_sample_traits.csv", row.names = F) ### ADDED group AND treatment MANUALLY AFTERWARD

##################################################################################################################
##################################################################################################################

### WCGNA ANALYSIS ###

#=====================================================================================
#
#  Set working directory and variables
#
#=====================================================================================

# Set working directory
setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA")


# Set variables
count_data <- "AN_variance_stabilizing_transformation.csv" # filtered, vst count data
sample_traits <- "ANTENNA_sample_traits.csv" # description of sample traits

#mytrait0 <- "TMB.vs.TFB" # trait of interest (column in datTraits)
#module = "tan" # module of interest
#module = "midnightblue" # module of interest
#modules = c("magenta", "midnightblue")
functional_annot <- "Hmel_blast_header_edited_table.csv" # B2G functional annotation
network_type <- "signed"

#=====================================================================================
#
#  Load required packages and set options for WGCNA
#
#=====================================================================================

library(WGCNA)
library(dplyr)
options(stringsAsFactors = FALSE)

#=====================================================================================
#
#  Load variance stabilizing transformed counts and transpose for WGCNA
#
#=====================================================================================

# Load brain vsd dataset from DESeq2
Data <- read.csv(count_data, row.names = 1)

# Take a quick look at what is in the data set
dim(Data)
names(Data) # column names (i.e., samples)

# Transpose the expression data (i.e., genes as columns and samples as rows)
datExpr0 = as.data.frame(t(Data))

#=====================================================================================
#
#  Filter out bad genes/samples and obvious outliers
#
#=====================================================================================

# Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# If above not "TRUE", remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
pdf(file = "ANTENNA_sampleClustering_2.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 90, col = "red") # Plot a line to show the cut
dev.off()

# Remove outliers, if necessary (this is directly from tutorial; cutHeight must be modified for data if used!)
# Plot a line to show the cut
#abline(h = 100, col = "red")
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)
datExpr <- datExpr0 # IF NO OUTLIERS, SO JUST RENAME TO datExpr0

### datExpr is now ready for network analysis ###

#=====================================================================================
#
#  Load sample characteristics dataset and cluster samples
#
#=====================================================================================

# Load sample characteristics
datTraits = read.csv(sample_traits, row.names = NULL) %>%
  dplyr::select(MT.vs.MC,
                RC.vs.MC,
                RT.vs.MT,
                RT.vs.RC,
                group, treatment) # select only comparisons of interest
dim(datTraits)
names(datTraits)

collectGarbage()

# Visualize how traits related to the sample dendrogram
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
pdf(file = "ANTENNA_sampleClustering_with_traits.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Save the relevant expression and trait data for use in the next steps 
save(datExpr, datTraits, file = "ANTENNA-dataInput.RData")

#=====================================================================================
#
#  Select soft-thresholding power and calculate adjacencies
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = network_type)

# Plot the results:
sizeGrWindow(9, 5)
pdf(file = "ANTENNA_soft_threshold_and_mean_connectivity.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Calculate adjacencies using softPower threshold based on soft_thresholding_and_mean_connectivity.pdf plot
softPower = 9 # soft-threshold power
adjacency = adjacency(datExpr, power = softPower, type = network_type)

#=====================================================================================
#
#  Topological Overlap Matrix (TOM)
#
#=====================================================================================

# Turn adjacency into topological overlap (to minimize effects of noise and spurious associations)
TOM = TOMsimilarity(adjacency, TOMType = network_type)
dissTOM = 1-TOM

collectGarbage()

#=====================================================================================
#
#  Clustering using TOM and plotting dendrograms
#
#=====================================================================================

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)

sizeGrWindow(12,9)
pdf(file = "ANTENNA_dendogram.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


# Set the minimum module size relatively high:
minModuleSize = 30

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
collectGarbage()
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

write.csv(table(dynamicColors), "ANTENNA_number_genes_per_module.csv")

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf(file = "ANTENNA_dendogram_with_colors.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#Removing any NaNs columns
#MEs<-MEs %>% select_if(~ !any(is.na(.)))

# Calculate dissimilarity of module eigengenes
MEDiss = 1-WGCNA::cor(MEs) #adding use = 'p' as I have a column of NaNs (MEgrey) as suggested by the WGCNA package

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
pdf(file = "ANTENNA_clustering_of_module_eigengenes.pdf", width = 12, height = 9)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Set height cut
MEDissThres = 0.25 # HEIGHT CUT CORRESPONDING CORRELATION OF 0.75
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Plot dendrogram with original and merged modules
sizeGrWindow(12, 9)
pdf(file = "ANTENNA_dendrogram_original_and_merged_modules.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ANTENNA-networkConstruction-stepByStep.RData")

#=====================================================================================
#
#  Quantifying module-trait associations 
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#moduleTraitCor = cor(MEs, datTraits, use = "p") # INCLUDE ALL PAIRWISE COMPARISONS
moduleTraitCor = WGCNA::cor(MEs, datTraits, use = "p") # INCLUDING ONLY PAIRWISE COMPARISONS OF INTEREST
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue_FDR <- p.adjust(moduleTraitPvalue, method = "fdr")  %>% # Adjust p-value using fdr method
  matrix(., nrow = nrow(moduleTraitPvalue), byrow = FALSE, 
         dimnames = list(row.names(moduleTraitPvalue), colnames(moduleTraitPvalue)))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n",
                   #"p = ", signif(moduleTraitPvalue, 1), "\n", # THESE ARE UNADJUSTED P-VALUES
                   "FDR = ", signif(moduleTraitPvalue_FDR, 1), sep = "") # THESE ARE ADJUSTED P-VALUES (FDR)
dim(textMatrix) = dim(moduleTraitCor)


# Display the correlation values within a heatmap plot
#pdf(file = "BRAIN_module-trait_relationships_heatmap_withouttext.pdf", width = 18, height = 15)
png(filename = "ANTENNA-Module-Trait Relationship.png", width = 18, height = 15, res=300, units = "in")
par(mar = c(6, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               #xLabels = names(datTraits), # FOR ALL PAIRWISE COMPARISONS
               xLabels = names(datTraits), # FOR ONLY PAIRWISE COMPARISONS OF INTEREST
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab = 1.5,
               #main = paste("Module-trait relationships"),
               zlim = c(-1,1))
dev.off()


# Save all variables and data for further analysis with "WGCNA_extract_data_for_important_modules.R"
save.image("Antenna_network_data.RData")
