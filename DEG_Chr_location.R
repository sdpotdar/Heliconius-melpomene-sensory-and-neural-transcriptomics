#Annotate DEGs with Chr locations.
library(dplyr)
require(data.table)
#Read in the genes_location csv, which will be common for all the Chr locations we do from now on and need not be loaded everytime

gene_location<-read.csv("Hmel2.5_genes_locations.csv", header = T) #Should have 20097 genes. Chr numbered from 0-21, with 0 being no chr assigned

#Read in the DEG file that needs to have Chr locations
#read in DEG file
AN_MRtrain<-read.csv("AN_malletitrained_vs_rosinatrained_DEG.csv")
AN_MRtrain_unique<-read.csv("AN_malletitrained_vs_rosinatrained_DEG_unique.csv")
AN_MRcontrol<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG.csv")
AN_MRcontrol_unique<-read.csv("AN_malleticontrol_vs_rosinacontrol_DEG_unique.csv")
AN_AEG<-read.csv("AN_AEGs.csv")

#Join the Chr location file with DEG file using the common column SeqName

n1<-left_join(AN_MRtrain, gene_location, by = "SeqName") # adds Gene location info to DEG file 
write.table(n1,file="AN_malletitrained_vs_rosinatrained_Chrlocations.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the Chr location per gene in a csv file.

n2<-left_join(AN_MRtrain_unique, gene_location, by = "SeqName") # adds Gene location info to DEG file 
write.table(n2,file="AN_malletitrained_vs_rosinatrained_unique_Chrlocations.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the Chr location per gene in a csv file.

n3<-left_join(AN_MRcontrol, gene_location, by = "SeqName") # adds Gene location info to DEG file 
write.table(n3,file="AN_malleticontrol_vs_rosinacontrol_Chrlocations.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the Chr location per gene in a csv file.

n4<-left_join(AN_MRcontrol_unique, gene_location, by = "SeqName") # adds Gene location info to DEG file 
write.table(n4,file="AN_malleticontrol_vs_rosinacontrol_unique_Chrlocations.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the Chr location per gene in a csv file.

n5<-left_join(AN_AEG, gene_location, by = "SeqName") # adds Gene location info to DEG file 
write.table(n5,file="AN_AEGs_annotated_with_Chrlocations.csv",sep=",",row.names=FALSE, quote=FALSE) #Writes the Chr location per gene in a csv file.


################################ Plot graphs #################################
library (ggplot2)
bar<-read.csv("Chromosome_specific_gene_expression.csv")
bar$Chromosome<-as.factor(bar$Chromosome)
aeg<-read.csv("Chromosome_specific_AEG.csv")
aeg$Chromosome<-as.factor(aeg$Chromosome)


a<-ggplot(data = aeg, aes(x=Chromosome, y=per_exp))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_classic()+xlab("Chromosomes")+ylab("% Expressed genes")+
  scale_x_discrete(breaks = aeg$Chromosome)+
  theme(axis.text = element_text(size = 15, colour = "black"), axis.title = element_text(size = 25, color = "black"),
        legend.position = "right", legend.title = element_text(size = 0), legend.text = element_text(size = 25))
a
png("/home/sushant/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/Chr_specific_DEGs/AN_AEGs.png",
    units="in",width=11.8, height=9, res=1280, bg="transparent")
a
dev.off()

g<-ggplot(data = bar, aes(x=Chromosome, y=percentage, fill=Treatment))+
  geom_bar(stat = "identity", position = "dodge")+ 
  scale_fill_manual("Treatment", values = c("Control"="sienna4", "Unique control" = "sienna1", 
                                            "Trained"="blue", "Unique trained"="cyan"))+
  theme_classic()+xlab("Chromosomes")+ylab("% Differentially expressed")+
  theme(axis.text = element_text(size = 15, colour = "black"), axis.title = element_text(size = 25, color = "black"),
        legend.position = "right", legend.title = element_text(size = 0), legend.text = element_text(size = 0))
g
png("/home/sushant/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads/Chr_specific_DEGs/AN_DEGs.png",
    units="in",width=15.8, height=9, res=1280, bg="transparent")
g
dev.off()

