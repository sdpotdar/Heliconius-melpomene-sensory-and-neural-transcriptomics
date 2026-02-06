setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads")
library(dplyr)

#We will load gene categories file (this has only seqname and category)
AN_MRtrain_unique_categories<-read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_categories.csv")


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

### Create courtship vector and add it to the sampleTable below
courtship<-c("NA", "NA", "Yes", "NA", "No", "NA", "NA", "Yes", "NA", "NA",
             "No", "NA", "NA", "No", "NA", "No", "Yes", "NA", "NA", "Yes",
             "No", "NA", "NA", "NA", "Yes", "Yes", "NA", "NA", "Yes", "No",
             "NA", "NA", "Yes", "NA", "No", "Yes", "NA")

# Create table containing information on all samples
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles2,
                          condition = sampleCondition,
                          treatment = treatment,
                          group = group,
                          courtship = courtship)

# Build DESeqDataSet using counts files we created
ddsHTSeq1 <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = directory,
                                        #design = ~ family + group)
                                        design = ~ condition)

ddsHTSeq1$condition <- relevel(ddsHTSeq1$condition, ref = "Rosina_Trained")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq1)) >= 10
ddsHTSeq1 <- ddsHTSeq1[keep,]

# Run differential expression analysis
dds <- DESeq(ddsHTSeq1)

# Obtain results for Malleti Trained vs Malleti Control
resultsNames(dds) # View coefficients
res.lfcShrink1 <- lfcShrink(dds, coef="condition_Malleti_Trained_vs_Rosina_Trained", type="apeglm") # Shink LFC for genes with small counts

# Perform variance stabilizing transformation for PCA 
vsd <- vst(dds, blind=FALSE)
deg.mat <- assay(vsd)
# Calculate Z-scores
deg.mat <- t(scale(t(deg.mat)))


######### Heatmap for EY_malletitrained_vs_rosinatrained_unique ###############
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendsort)
library(tibble)
library(cowplot)

sig_genes <- read.delim("malletitrained_vs_rosinatrained_unique_degs.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat2 <- as.data.frame(deg.mat[grepl(paste(sig_genes, collapse='|'),
                                         row.names(deg.mat)), ]) %>%
  dplyr::select(matches(c("*._M_T_.*", "*._R_T_.*")))

# Add row category annotations (DOUBLE CHECK THAT THESE ARE CORRECT!!!)
#deg.mat2 <- as.data.frame(deg.mat2) %>% 
#  rownames_to_column(., var = "SeqName") %>% 
#  left_join(., EY_MRtrain_unique_gene_category[,c("SeqName", "Category")], by = "SeqName")

#deg.mat2<-deg.mat2[-c(23), ] #Do only once.. Had to do this cause the genes in my set were not matching with the ones obtained here.

#deg.mat2 <- as.data.frame(deg.mat[grepl(paste(sig_genes, collapse='|'),
#                                       row.names(deg.mat)), ]) %>%
#  dplyr::select(matches(c("*._M_T_.*", "*._R_T_.*")))

# Specify column annotation data for heatmap
annot <- colData(dds)
annot <- as.data.frame(annot) %>% dplyr::select(courtship, group)
colnames(annot)[1] <- "Courtship"
colnames(annot)[2] <- "Group"

# Specify row annotation data for heatmap
annotrow<- read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_categories.csv", row.names=1) %>% # row.names=1 is making first column into row names instead of 1, 2, 3...
  dplyr::select(Category) ## you have to make sure that rownames(row of this) == rownames(col of perm_deg_mat). This means exact same ORDER of genes (ascending usually)
annotrow$Category[annotrow$Category == "Learning"]<- "Learning"
annotrow$Category[annotrow$Category == "Eye"] <- "Eye"
annotrow$Category[annotrow$Category == "Olfaction"] <- "Olfaction"
annotrow$Category[annotrow$Category == "Neural"] <- "Neural"
annotrow$Category[annotrow$Category == "Hormones"] <- "Hormone"
annotrow$Category[annotrow$Category == "Zinc"] <- "Zinc"
annotrow$Category[annotrow$Category == "Other"] <- "Other"
annotrow$Category[annotrow$Category == "Sex"] <- "Sex"
annotrow$Category[annotrow$Category == "Wing"] <- "Wing"
annotrow$Category[annotrow$Category == "Uncharacterized"] <- "Uncharacterized"
colnames(annotrow)[1] <- "Gene categories"

# Create list with colors for each annotation
annotation_colors2 <- list(Courtship = c(Yes="darkgreen", No="#D41159", "NA"="#FEFE62"),
                           Group = c(Malleti="orange", Rosina="black")) # col colours
annotation_colors3<- list("Gene categories" = c("Learning"="green", 
                                                "Eye" ="orangered2", 
                                                "Olfaction"="mediumorchid", 
                                                "Neural"="yellow1", 
                                                "Hormone"="royalblue1", 
                                                "Zinc"="deeppink", 
                                                "Other"="black", 
                                                "Sex"="tan1",
                                                "Wing"="darkslateblue",
                                                "Uncharacterized"="grey")) # row colours
annotation_colors4<- c(annotation_colors2, annotation_colors3)


# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat2)/paletteLength, max(deg.mat2), length.out=floor(paletteLength/2)))

# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

malleti_vs_rosina_trained_unique <-pheatmap(deg.mat2,
                                            #main = "Malleti vs Rosina Trained Unique \n",
                                            color = my_palette,
                                            border_color = NA,
                                            annotation_colors = annotation_colors4,
                                            cluster_rows=TRUE, 
                                            cluster_cols=TRUE, 
                                            cellwidth = 8,
                                            cellheight = NA,
                                            cutree_rows = 2,
                                            cutree_cols = 2,
                                            treeheight_row = 50,
                                            treeheight_col = 30,
                                            annotation_col=annot,
                                            annotation_row=annotrow,
                                            annotation_names_col = FALSE, # deletes row/ col names
                                            annotation_names_row = FALSE,
                                            annotation_legend = F,
                                            clustering_callback = callback,
                                            show_colnames = FALSE, # deletes individual row/ col cluster names
                                            show_rownames = FALSE,
                                            breaks = myBreaks,
                                            fontsize = 15)# affects scale of exported heatmap png,
malleti_vs_rosina_trained_unique
png(file= "malleti_vs_rosina_trained_unique_heatmap.png", height = 1700, width = 1920, res=300)
malleti_vs_rosina_trained_unique
dev.off()

########################################################################################################
################# Malleti control vs rosina control unique

ddsHTSeq1$condition <- relevel(ddsHTSeq1$condition, ref = "Rosina_Control")

# Pre-filter to keep only rows (i.e., genes) that have at least 10 reads total
keep <- rowSums(counts(ddsHTSeq1)) >= 10
ddsHTSeq1 <- ddsHTSeq1[keep,]

# Run differential expression analysis
dds <- DESeq(ddsHTSeq1)

# Obtain results for Malleti Trained vs Malleti Control
resultsNames(dds) # View coefficients
res.lfcShrink1 <- lfcShrink(dds, coef="condition_Malleti_Control_vs_Rosina_Control", type="apeglm") # Shink LFC for genes with small counts

# Perform variance stabilizing transformation for PCA 
vsd <- vst(dds, blind=FALSE)
deg.mat <- assay(vsd)
# Calculate Z-scores
deg.mat <- t(scale(t(deg.mat)))


######### Heatmap for EY_malleticontrol_vs_rosinacontrol_unique ###############
library(pheatmap)
library(dendsort)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dendsort)
library(tibble)
library(cowplot)

sig_genes <- read.delim("malleticontrol_vs_rosinacontrol_degs_unique.txt", sep = "\t", header = FALSE) %>%
  .[,1]  ### reading in txt files, [,1] shows column 1
### please remember to change the dplyr::select part to respective file names #
deg.mat2 <- as.data.frame(deg.mat[grepl(paste(sig_genes, collapse='|'),
                                        row.names(deg.mat)), ]) %>%
  dplyr::select(matches(c("*._M_C_.*", "*._R_C_.*")))

# Add row category annotations (DOUBLE CHECK THAT THESE ARE CORRECT!!!)
#deg.mat2 <- as.data.frame(deg.mat2) %>% 
#  rownames_to_column(., var = "SeqName") %>% 
#  left_join(., EY_MRtrain_unique_gene_category[,c("SeqName", "Category")], by = "SeqName")

#deg.mat2<-deg.mat2[-c(23), ] #Do only once.. Had to do this cause the genes in my set were not matching with the ones obtained here.

#deg.mat2 <- as.data.frame(deg.mat[grepl(paste(sig_genes, collapse='|'),
#                                       row.names(deg.mat)), ]) %>%
#  dplyr::select(matches(c("*._M_T_.*", "*._R_T_.*")))

# Specify column annotation data for heatmap
annot <- colData(dds)
annot <- as.data.frame(annot) %>% dplyr::select(group)
#colnames(annot)[1] <- "Courtship"
colnames(annot)[1] <- "Group"

# Specify row annotation data for heatmap
# annotrow<- read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_categories.csv", row.names=1) %>% # row.names=1 is making first column into row names instead of 1, 2, 3...
#   dplyr::select(Category) ## you have to make sure that rownames(row of this) == rownames(col of perm_deg_mat). This means exact same ORDER of genes (ascending usually)
# annotrow$Category[annotrow$Category == "Learning"]<- "Learning"
# annotrow$Category[annotrow$Category == "Eye"] <- "Eye"
# annotrow$Category[annotrow$Category == "Olfaction"] <- "Olfaction"
# annotrow$Category[annotrow$Category == "Neural"] <- "Neural"
# annotrow$Category[annotrow$Category == "Hormones"] <- "Hormone"
# annotrow$Category[annotrow$Category == "Zinc"] <- "Zinc"
# annotrow$Category[annotrow$Category == "Other"] <- "Other"
# annotrow$Category[annotrow$Category == "Sex"] <- "Sex"
# annotrow$Category[annotrow$Category == "Wing"] <- "Wing"
# annotrow$Category[annotrow$Category == "Uncharacterized"] <- "Uncharacterized"
# colnames(annotrow)[1] <- "Gene categories"

# Create list with colors for each annotation
annotation_colors2 <- list(Courtship = c(Yes="darkgreen", No="#D41159", "NA"="#FEFE62"),
                           Group = c(Malleti="orange", Rosina="black")) # col colours
# annotation_colors3<- list("Gene categories" = c("Learning"="green", 
#                                                 "Eye" ="orangered2", 
#                                                 "Olfaction"="mediumorchid", 
#                                                 "Neural"="yellow1", 
#                                                 "Hormone"="royalblue1", 
#                                                 "Zinc"="deeppink", 
#                                                 "Other"="black", 
#                                                 "Sex"="tan1",
#                                                 "Wing"="darkslateblue",
#                                                 "Uncharacterized"="grey")) # row colours
# annotation_colors4<- c(annotation_colors2, annotation_colors3)


# Create color palette and set breaks so that midpoint (0) is white
paletteLength <- 1000
my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(deg.mat2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(deg.mat2)/paletteLength, max(deg.mat2), length.out=floor(paletteLength/2)))

# Plot heatmap

callback = function(hc, ...){dendsort(hc)}

malleti_vs_rosina_control_unique <-pheatmap(deg.mat2,
                                            #main = "Malleti vs Rosina Trained Unique \n",
                                            color = my_palette,
                                            border_color = NA,
                                            annotation_colors = annotation_colors2,
                                            cluster_rows=TRUE, 
                                            cluster_cols=TRUE, 
                                            cellwidth = 8,
                                            cellheight = NA,
                                            cutree_rows = 2,
                                            cutree_cols = 2,
                                            treeheight_row = 50,
                                            treeheight_col = 30,
                                            annotation_col=annot,
                                            #annotation_row=annotrow,
                                            annotation_names_col = FALSE, # deletes row/ col names
                                            annotation_names_row = FALSE,
                                            annotation_legend = F,
                                            clustering_callback = callback,
                                            show_colnames = FALSE, # deletes individual row/ col cluster names
                                            show_rownames = FALSE,
                                            breaks = myBreaks,
                                            fontsize = 15)# affects scale of exported heatmap png,
malleti_vs_rosina_trained_unique
png(file= "malleti_vs_rosina_control_unique_heatmap.png", height = 1700, width = 1920, res=300)
malleti_vs_rosina_control_unique
dev.off()
