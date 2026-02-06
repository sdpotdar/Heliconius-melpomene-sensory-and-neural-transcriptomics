#Annotate DEGs with Hmel_2.5_blastp_uniport.csv annotation file.
library(dplyr)
require(data.table)
#Read in the blastp file, which will be common for all the annotations we do from now on and need not be loaded everytime we do annotation for different gene lists

functional_annot <- "~/Documents/uark/phd/projects/melpomene_rnaseq/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.proteins.fa.blastp.uniprot_sprot_modified.1e-10.tsv" # Blast uniport annotation file from lepbase
annot<-as.data.frame(fread(functional_annot)) #.tsv file is too huge to be loaded into R. No of rows = 1883420 hence using fread
#Adding column names to the data frame annot
names(annot)[1]<-"SeqName"
names(annot)[2]<-"Organism"
names(annot)[3]<-"Blast_percent"
names(annot)[11]<-"Evalue"
names(annot)[15]<-"Description"
head(annot) #Check to see if the new column names appear
need<-annot%>%dplyr::select(SeqName, Description, Blast_percent, Evalue, Organism) 

#Read in the best blast hit file that was generated using Blast2go and DIAMOND
best_blast_hit<- "~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/blastx_results_header_edited_xml/Hmel_blast_header_edited_table.tsv" #Blastx file from Blast2Go and DIAMOND
blastx<-as.data.frame(fread(best_blast_hit)) #.tsv file is too huge to be loaded into R. No of rows = 20097 hence using fread

head(blastx)
names(blastx)[6]<-"Best_blast_hit"
names(blastx)[10]<-"GO_IDs"
names(blastx)[11]<-"GO_Names"
need<-annot%>%dplyr::select(SeqName, Description, Blast_percent, Evalue, Organism) 
need2<-blastx%>%dplyr::select(SeqName, Best_blast_hit, GO_IDs, GO_Names)

#Read in the DEG/AEG file that needs to be annotated
#read in DEG file
#setwd
setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads")
AN_AEGs<-read.csv("AN_AEGs.csv")#Make sure that it has the SeqName as a header for the first column
head(AN_AEGs)#Check the data
AN_AEGs<-data.frame(AN_AEGs) #Change it into a dataframe
names(AN_AEGs)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
n2<-left_join(AN_AEGs, need, by = "SeqName") # adds functional annot (description and blast hits)
write.table(n2,file="AN_AEGs_annotated.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the annotations into a csv file.

#Read in the AN_malletitrained_vs_rosinatrained_DEGs.csv file
AN_MRtrain<-read.csv("AN_malletitrained_vs_rosinatrained_DEG.csv")
head(AN_MRtrain)
AN_MRtrain<-data.frame(AN_MRtrain) #Change it into a dataframe
names(AN_AEGs)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
n3<-left_join(AN_MRtrain, need, by = "SeqName") # adds functional annot (description and blast hits)
write.table(n3,file="AN_malletitrained_vs_rosinatrained_DEG_annotated.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the annotations into a csv file.


#Read in the AN_malleticontrol_vs_rosinacontrol_DEGs.csv file
AN_MRcontrol<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG.csv")
head(AN_MRcontrol)
AN_MRcontrol<-data.frame(AN_MRcontrol) #Change it into a dataframe
names(AN_MRcontrol)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
n4<-left_join(AN_MRcontrol, need, by = "SeqName") # adds functional annot (description and blast hits)
write.table(n4,file="AN_malleticontrol_vs_rosinacontrol_DEG_annotated.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the annotations into a csv file.

#Merge DEGs SeqName unique to trained treatment with Log2foldchange, pvalue, padj and the 'need' columns
AN_MRtrain_unique<-read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique.csv")
names(AN_MRtrain_unique)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
#Read in the annotated_removed_duplicate csv file for malleti trained vs rosina trained
AN_MRtrain_annoduprem<-read.csv("AN_malletitrained_vs_rosinatrained_DEG_annotated_removed_duplicate.csv")
AN_MRtrain_unique<-data.frame(AN_MRtrain_unique)
need2<-AN_MRtrain_annoduprem%>%dplyr::select(SeqName, baseMean, log2FoldChange, lfcSE, pvalue, padj, Description, Blast_percent, Evalue, Organism)
n5<-left_join(AN_MRtrain_unique, need2, by = "SeqName") # adds functional annot (description and blast hits)
write.table(n5,file="AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the annotations into a csv file.


#Merge DEGs SeqName unique to control treatment with Log2foldchange, pvalue, padj and the 'need' columns
AN_MRcontrol_unique<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG_unique.csv")
names(AN_MRcontrol_unique)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
#Read in the annotated_removed_duplicate csv file for malleti control vs rosina control
AN_MRcontrol_annoduprem<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG_annotated_removed_duplicate.csv")
need3<-AN_MRcontrol_annoduprem%>%dplyr::select(SeqName, baseMean, log2FoldChange, lfcSE, pvalue, padj, description, Blast.description, Blast_percent, Evalue, Organism)
AN_MRcontrol_unique<-data.frame(AN_MRcontrol_unique)
n7<-left_join(AN_MRcontrol_unique, need3, by = "SeqName") # adds functional annot (description and blast hits)
write.table(n7,file="AN_malleticontrol_vs_rosinacontrol_DEG_unique_annotated_removed_duplicate.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the annotations into a csv file.


#Add gene categories from EY and BR files
#Read in the Brain category file as well, so that we can add the common gene categories between BR and EY
BR_gene_category<-read.csv("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/BR_aligned_reads/BR_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_categories.csv") # modified interproscan file from lepbase
need3<-BR_gene_category%>%dplyr::select(SeqName, Category)
#Read in the Eye category file as well, so that we can add the common gene categories between BR and EY
EY_gene_category<-read.csv("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/EY_aligned_reads/EY_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_gene_categories.csv") # modified interproscan file from lepbase
need4<-EY_gene_category%>%dplyr::select(SeqName, Category)

###Add Best blast hit, GO IDs and GO names to the important already annotated by Hmel annot files
#Read in the AN_malletitrained_vs_rosinatrained_DEGs.csv unique file with added annotation from the genomefile
AN_MRtrain_unique_annotated_rem_dup<-read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate.csv")
head(AN_MRtrain_unique_annotated_rem_dup)
AN_MRtrain_unique_annotated_rem_dup<-data.frame(AN_MRtrain_unique_annotated_rem_dup) #Change it into a dataframe
n3<-left_join(AN_MRtrain_unique_annotated_rem_dup, need2, by = "SeqName") # adds functional annot (description and blast hits)
n4<-left_join(n3, need3, by = "SeqName") # adds BR gene categories 
n5<-left_join(n4, need4, by = "SeqName") # adds EY gene categories
write.table(n5,file="AN_malletitrained_vs_rosinatrained_DEG_unique_annotated_removed_duplicate_bestblasthit_GOterms.xlsx",sep="$",row.names=FALSE, quote=FALSE) #Writes the best blast hit and GO names into a new csv file.

###Add Best blast hit, GO IDs and GO names to the important already annotated by Hmel annot files
#Read in the BR_malletitrained_vs_rosinatrained_DEGs.csv file with added annotation from the genome
EY_MRtrain_annotated_rem_dup<-read.csv("EY_malletitrained_vs_rosinatrained_DEG_annotated_removed_duplicate.csv")
head(EY_MRtrain_annotated_rem_dup)
EY_MRtrain_annotated_rem_dup<-data.frame(EY_MRtrain_annotated_rem_dup) #Change it into a dataframe
n5<-left_join(EY_MRtrain_annotated_rem_dup, need2, by = "SeqName") # adds functional annot (description and blast hits)
n6<-left_join(n5, need4, by = "SeqName")
write.table(n6,file="EY_malletitrained_vs_rosinatrained_DEG_annotated_removed_duplicate_bestblasthit_GOterms.xlsx",sep="$",row.names=FALSE, quote=FALSE) #Writes the best blast hit and GO names into a new csv file.



###Add Best blast hit, GO IDs and GO names to the important already annotated by Hmel annot files
#Read in the AN_malleticontrol_vs_rosinacontrol_DEGs.csv unique file with added annotation from the genomefile
AN_MRcontrol_unique_annotated_rem_dup<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG_unique_annotated_removed_duplicate.csv")
head(AN_MRcontrol_unique_annotated_rem_dup)
AN_MRcontrol_unique_annotated_rem_dup<-data.frame(AN_MRcontrol_unique_annotated_rem_dup) #Change it into a dataframe
n3<-left_join(AN_MRcontrol_unique_annotated_rem_dup, need2, by = "SeqName") # adds functional annot (description and blast hits)
n4<-left_join(n3, need3, by = "SeqName") # adds BR gene categories 
n5<-left_join(n4, need4, by = "SeqName") # adds EY gene categories
write.table(n3,file="AN_malleticontrol_vs_rosinacontrol_DEG_unique_annotated_removed_duplicate_bestblasthit_GOterms.xlsx",sep="$",row.names=FALSE, quote=FALSE) #Writes the best blast hit and GO names into a new csv file.
