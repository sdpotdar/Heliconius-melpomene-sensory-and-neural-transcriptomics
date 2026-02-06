#Add interproscan and GO IDs to the _annotated_removed_duplicate AND _unique_annotated_removed_duplicate DEG files (These are the files in which annotation from blastp uniport is already present)
library(dplyr)
require(data.table)

#Read in the interproscan file (_modified_interproscan.tsv), which will be common for all the additions we do from now on and need not be loaded everytime we do annotation for different gene lists
interproscan <- "~/Documents/uark/phd/projects/melpomene_rnaseq/H_melpomene_genome/Heliconius_melpomene_melpomene_Hmel2.5.proteins.fa_modified_interproscan.tsv" # modified interproscan file from lepbase
inter_annot<-as.data.frame(fread(interproscan)) #.tsv file is too huge to be loaded into R. No of rows = 24545 hence using fread
#Adding column names to the data frame annot
names(inter_annot)[1]<-"SeqName"
names(inter_annot)[3]<-"length"
names(inter_annot)[4]<-"database"
names(inter_annot)[12]<-"Interproscan_ID"
names(inter_annot)[13]<-"Interproscan_name"
names(inter_annot)[14]<-"GO_ID"
names(inter_annot)[15]<-"Reactome_ID"
head(inter_annot) #Check to see if the new column names appear
need<-inter_annot%>%dplyr::select(SeqName, database, Interproscan_ID, Interproscan_name, GO_ID, Reactome_ID) #These are the columns that we need to add to the annotated DEG/AEG file.

#Read in the DEG/AEG file that needs to be annotated
#read in DEG file
#setwd
setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads")
AN_AEGs<-read.csv("AN_AEGs_annotated.csv")#Make sure that it has the SeqName as a header for the first column
head(AN_AEGs)#Check the data
AN_AEGs<-data.frame(AN_AEGs) #Change it into a dataframe
names(AN_AEGs)[1]<-"SeqName" #Adding SeqName as header as all DEG/AEG files do not have that!
n2<-left_join(AN_AEGs, need, by = "SeqName", relationship = "many-to-many") # adds interproscan_id, GO_ID, and Reactome_ID to the annotated AEG/DEG file based on SeqName)
write.table(n2,file="AN_AEGs_annotated_interproscan_goid.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the file with interproscan and GO_ID into a csv file.
