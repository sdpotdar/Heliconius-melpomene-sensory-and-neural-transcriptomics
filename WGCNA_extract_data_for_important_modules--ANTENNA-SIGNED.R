setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA")
wd<-setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/AN_WGCNA")

library(WGCNA)
library(dplyr)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


### Load network analysis (from "BRAIN_binarizeCategoricalVariable_for_pairwise_comparisons_WGCNA_2.R")
load(file = "Antenna_network_data.RData")
annot <- read.csv("Hmel_description.csv", sep = "\t", header = TRUE) 
#=====================================================================================
#
#  Gene relationship to trait and important modules: Gene Significance and Module Membership
#
#=====================================================================================
# Define variable containing the specified trait column of datTrait
mytrait0 <- "MT.vs.MC" #trait of interest; column in datTraits
module = "turquoise" #module of interest

# Create function to extract module-trait info and create Cytoscape output <<< NEED TO CHECK THAT THIS IS CORRECT!!!
#traitmodule <- function(mytrait0, module){
  #=====================================================================================
  ### Gene relationship to trait and important modules: Gene Significance and Module Membership ###
  # GS = the absolute value of the correlation between the gene and the trait
  # MM = the correlation of the module eigengene and the gene expression profile
  # Define variable containing the specified trait column of datTrait
  mytrait = as.data.frame(datTraits[,mytrait0])
  names(mytrait) = mytrait0
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(cor(datExpr, mytrait, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", mytrait0, sep="")
  names(GSPvalue) = paste("p.GS.", mytrait0, sep="")
  
  #=====================================================================================
  ### Intramodular analysis: identifying genes with high GS and MM ###
  # Creates scatterplot to determine if GS and MM are highly correlated
  # Genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  sizeGrWindow(7, 7)
  pdf(file = paste0("ANTENNA_signed", mytrait0, "_", module, "_module-MM_vs_GS_scatterplot.pdf"), width = 7, height = 7)
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module membership in", module, "module"),
                     ylab = paste("Gene significance for", mytrait0),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
  #=====================================================================================
  ### Summary output of network analysis results ###
  #names(datExpr) # print all gene names in datExpr (all that were included in the network analysis)
  print(paste0("Genes in ", module, " module:"))
  print(names(datExpr)[moduleColors==module]) # Print names of genes in module of interest
  

  ##Check the functional annot file 
  dim(annot)
  names(annot)
  probes = names(datExpr)
  probes2annot = match(probes, annot$SeqName)
  
  # The following is the number or probes without annotation and should return 0
  probe_no_anno <- sum(is.na(probes2annot))
  ifelse(probe_no_anno == 0,
         print("All genes were annotated successfully!"),
         print("WARNING: Some genes were not annotated!"))
  
  # Create a data frame holding gene ID, gene description, best blast hit, module color, GS for the trait of interest,
  # and module membership and p-values in all modules. Module are ordered by GS for the trait of interest, with the
  # most significant ones to the left.
  # Create the starting geneInfo data frame
  geneInfo0 = data.frame(SeqName = probes,
                         Description = annot$Description[probes2annot],
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  
  # Order modules by their significance for trait of interest
  modOrder = order(-abs(cor(MEs, mytrait, use = "p")))
  
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]])
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[, paste0("GS.", mytrait0)]))
  geneInfo = geneInfo0[geneOrder, ]
  
  # Write gene information to file
  write.csv(geneInfo, file = paste0("ANTENNA_geneInfo_", mytrait0, ".csv"))
  #=====================================================================================
  
  ### Exporting files for Cytoscape ###
  # Recalculate topological overlap if needed
  #TOM = TOMsimilarityFromExpr(datExpr, power = softPower)
  
  # Select module(s) of interest
  module = module
  
  # Select module probes
  probes = names(datExpr)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modGenes = annot$SeqName[match(modProbes, annot$SeqName)]
  
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("ANTENNA_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("ANTENNA_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])

  ### Identify hub genes for each module (DE)
  colorh = moduleColors
  hub = chooseTopHubInEachModule(datExpr, colorh, type = network_type, power = softPower)
  # Get annotations for hub genes for each module
  #read.csv(functional_annot) %>% filter(SeqName == hubs[["midnightblue"]])
  hub2 <- as.data.frame(hub) %>%
    tibble::rownames_to_column("Module") %>%
    rename(SeqName = hub) %>%
    inner_join(., annot, by = "SeqName")
write.csv(hub2, file = "ANTENNA_geneInfo_hubgenes.csv")
  


  