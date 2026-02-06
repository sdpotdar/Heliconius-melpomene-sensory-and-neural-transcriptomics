setwd("~/Documents/uark/phd/projects/melpomene_rnaseq/bioinformatics/AN_aligned_reads")
###Venn diamgrams
library(VennDiagram)

MRcontrol<- read.csv("AN_malleticontrol_vs_rosinacontrol_DEG.csv")
MRtrain<- read.csv("AN_malletitrained_vs_rosinatrained_DEG.csv")
MMct<- read.csv("AN_malletitrained_vs_malleticontrol_DEG.csv")
RRct<- read.csv("AN_rosinatrained_vs_rosinacontrol_DEG.csv")
names(MRcontrol)[1] <- "SeqName"
names(MRtrain)[1] <- "SeqName"
names(MMct)[1] <- "SeqName"
names(RRct)[1] <- "SeqName"


## for general (all) genes
set1<-MRcontrol$SeqName
set2<-MRtrain$SeqName
set3<-MMct$SeqName
set4<-RRct$SeqName

grid.newpage()
v<-venn.diagram(x=list(set1, set2, set3, set4),
                category.names = c("Malleti vs Rosina\ncontrol", "Malleti vs Rosina\ntrained", "Control vs Trained\nMalleti", "Control vs Trained\nRosina"),
                fontfamily="sans",
                cex=1.6,
                # customise cat names and fonts
                cat.fontface="bold",
                cat.default.pos="outer",
                cat.dist = c(0.08, 0.08, 0.1, 0.1),
                cat.cex=1,
                # fill of venn diagram circles
                fill=c("#BAB6B6", "#FFD370", "#6EBEEB", "#50A68E"),
                filename="AN_venn_diagram_all_combinations.png",# saves diagram - if not, write = NULL)
                height = 8,
                width = 8,
                units = "in",
                resolution = 300)
v
grid.draw(v) ## can use grid.draw if filename is NULL, with no height, width, units etc with it
