
rm(list = ls())
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/metaR")
library(clusterSim)
library(tidyr)

# ###############ref1
# ###############ref1
# ###############ref1
# ###############ref1
# ###############ref1
kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,3:19),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,dim(phylum)[2]-meta.col+2:dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,dim(class)[2]-meta.col+2:dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,dim(ord)[2]-meta.col+2:dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,dim(family)[2]-meta.col+2:dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,dim(genus)[2]-meta.col+2:dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,dim(species)[2]-meta.col+2:dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref1 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref1.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/Ref1.meta.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref1.meta),paste("ref1.",rownames(All.tax.ref1),sep="")) != c(1:dim(All.tax.ref1)[1]))

colnames(All.tax.ref1) <- colnames(All.tax.ref1) %>% gsub("d__","k__",.) 
All.tax.ref1.nor.com <- cbind(ref1.meta[,c("Sex","Age","IBS_severity","HADS_distress","Depression")],All.tax.ref1)
All.tax.ref1.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref1.nor.com$Depression )
All.tax.ref1.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref1.nor.com$Depression )

colnames(All.tax.ref1.nor.com)[1:5] <- c("Gender","Age","IBS","Distress","Depression")

colnames(All.tax.ref1.nor.com) <- gsub(";",".",colnames(All.tax.ref1.nor.com))
colnames(All.tax.ref1.nor.com) <- gsub("-","_",colnames(All.tax.ref1.nor.com))
colnames(All.tax.ref1.nor.com) <- gsub("\\[","",colnames(All.tax.ref1.nor.com))
colnames(All.tax.ref1.nor.com) <- gsub("\\]","",colnames(All.tax.ref1.nor.com))

library(metamicrobiomeR) 
taxacom.alltax.ref1<-taxa.compare(taxtab=All.tax.ref1.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Gender","Age","IBS"),longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref1,"Ref1.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref1.crude<-taxa.compare(taxtab=All.tax.ref1.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")


write.table(taxacom.alltax.ref1.crude,"Ref1.metaR.genus.silva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref1, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref1, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref1, tax.lev="l7",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)




###############ref23
###############ref23
###############ref23
###############ref23
###############ref23
kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-1.csv",header = T,check.names = F,colClasses=c("index"="character"))
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-2.csv",header = T,check.names = F,colClasses=c("index"="character"))
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-3.csv",header = T,check.names = F,colClasses=c("index"="character"))
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-4.csv",header = T,check.names = F,colClasses=c("index"="character"))
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-5.csv",header = T,check.names = F,colClasses=c("index"="character"))
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-7.csv",header = T,check.names = F,colClasses=c("index"="character"))

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,6:18),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref23 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref23.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/Ref23.meta.dada.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref23.meta),paste("Ref23.",rownames(All.tax.ref23),sep="")) != c(1:dim(All.tax.ref23)[1]))

colnames(All.tax.ref23) <- colnames(All.tax.ref23) %>% gsub("d__","k__",.) 
All.tax.ref23.nor.com <- cbind(ref23.meta[,c("Sex","Age","BMI","Depression")],All.tax.ref23)

All.tax.ref23.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref23.nor.com$Depression )
All.tax.ref23.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref23.nor.com$Depression )

colnames(All.tax.ref23.nor.com)[1:4] <- c("Gender","Age","BMI","Depression")
library(metamicrobiomeR) 
colnames(All.tax.ref23.nor.com) <- gsub(";",".",colnames(All.tax.ref23.nor.com))
colnames(All.tax.ref23.nor.com) <- gsub("-","_",colnames(All.tax.ref23.nor.com))
colnames(All.tax.ref23.nor.com) <- gsub("\\[","",colnames(All.tax.ref23.nor.com))
colnames(All.tax.ref23.nor.com) <- gsub("\\]","",colnames(All.tax.ref23.nor.com))
taxacom.alltax.ref23<-taxa.compare(taxtab=All.tax.ref23.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Gender","Age","BMI"),longitudinal="no",p.adjust.method="fdr")

write.table(taxacom.alltax.ref23,"ref23.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref23.crude<-taxa.compare(taxtab=All.tax.ref23.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref23.crude,"Ref23.metaR.genus.silva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref23, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref23, tax.lev="l2",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref23, tax.lev="l7",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)




###############ref27
###############ref27
###############ref27
###############ref27
###############ref27
kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-1.csv",header = T,check.names = F,colClasses=c("index"="character"))
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-2.csv",header = T,check.names = F,colClasses=c("index"="character"))
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-3.csv",header = T,check.names = F,colClasses=c("index"="character"))
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-4.csv",header = T,check.names = F,colClasses=c("index"="character"))
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-5.csv",header = T,check.names = F,colClasses=c("index"="character"))
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-7.csv",header = T,check.names = F,colClasses=c("index"="character"))

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index



length.1 <- dim(kingdom)[2]
colnames(kingdom)
ref27.meta <- kingdom[,c(1,5:6)]
kingdom <- kingdom[,-c(1,5:6),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref27 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

sum(pmatch(rownames(ref27.meta),rownames(All.tax.ref27)) != c(1:dim(All.tax.ref27)[1]))

ref27.meta <- ref27.meta[!(ref27.meta$Group == "BD" & ref27.meta$Time ==2),]
All.tax.ref27 <- All.tax.ref27[rownames(ref27.meta),]
colnames(All.tax.ref27) <- colnames(All.tax.ref27) %>% gsub("d__","k__",.) 

All.tax.ref27.nor.com <- cbind(ref27.meta[,-1],All.tax.ref27)


All.tax.ref27.nor.com$Depression <- "1_depression"
All.tax.ref27.nor.com$Depression[All.tax.ref27.nor.com$Group == "H"] <- "0_depression"

library(metamicrobiomeR) 
colnames(All.tax.ref27.nor.com) <- gsub(";",".",colnames(All.tax.ref27.nor.com))
colnames(All.tax.ref27.nor.com) <- gsub("-","_",colnames(All.tax.ref27.nor.com))
colnames(All.tax.ref27.nor.com) <- gsub("\\[","",colnames(All.tax.ref27.nor.com))
colnames(All.tax.ref27.nor.com) <- gsub("\\]","",colnames(All.tax.ref27.nor.com))

taxacom.tax_BD_H<-taxa.compare(taxtab=All.tax.ref27.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")

taxacom.alltax.ref27 <- taxacom.tax_BD_H
write.table(taxacom.alltax.ref27,"ref27.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)


taxcomtab.show(taxcomtab=taxacom.alltax.ref27, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref27, tax.lev="l2",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref27, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)




###############ref40
###############ref40
###############ref40
###############ref40
###############ref40


kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,5:59),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref40 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref40.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/Ref40.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref40.meta),paste("Ref40.",rownames(All.tax.ref40),sep="")) != c(1:dim(All.tax.ref40)[1]))

colnames(All.tax.ref40) <- colnames(All.tax.ref40) %>% gsub("d__","k__",.) 

All.tax.ref40.nor.com <- cbind(ref40.meta[,c("Sex","Age","BMI","Depression")],All.tax.ref40)
All.tax.ref40.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref40.nor.com$Depression )
All.tax.ref40.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref40.nor.com$Depression )

colnames(All.tax.ref40.nor.com)[1:4] <- c("Gender","Age","BMI","Depression")
colnames(All.tax.ref40.nor.com) <- gsub(";",".",colnames(All.tax.ref40.nor.com))
colnames(All.tax.ref40.nor.com) <- gsub("-","_",colnames(All.tax.ref40.nor.com))
colnames(All.tax.ref40.nor.com) <- gsub("\\[","",colnames(All.tax.ref40.nor.com))
colnames(All.tax.ref40.nor.com) <- gsub("\\]","",colnames(All.tax.ref40.nor.com))
library(metamicrobiomeR) 
taxacom.alltax.ref40<-taxa.compare(taxtab=All.tax.ref40.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI"),longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref40,"ref40.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref40.crude<-taxa.compare(taxtab=All.tax.ref40.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref40.crude,"ref40.metaR.genus.silva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref40, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref40, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref40, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)


###############ref60
###############ref60
###############ref60
###############ref60
###############ref60


kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,5:26),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref60 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref60.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/Ref60.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref60.meta),paste("Ref60.",rownames(All.tax.ref60),sep="")) != c(1:dim(All.tax.ref60)[1]))

colnames(All.tax.ref60) <- colnames(All.tax.ref60) %>% gsub("d__","k__",.) 

All.tax.ref60.nor.com <- cbind(ref60.meta[,c("Sex","Age","BMI","Depression","Antidepression")],All.tax.ref60)
All.tax.ref60.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref60.nor.com$Depression )
All.tax.ref60.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref60.nor.com$Depression )

colnames(All.tax.ref60.nor.com)[1:5] <- c("Gender","Age","BMI","Depression","Antidepression")
colnames(All.tax.ref60.nor.com) <- gsub(";",".",colnames(All.tax.ref60.nor.com))
colnames(All.tax.ref60.nor.com) <- gsub("-","_",colnames(All.tax.ref60.nor.com))
colnames(All.tax.ref60.nor.com) <- gsub("\\[","",colnames(All.tax.ref60.nor.com))
colnames(All.tax.ref60.nor.com) <- gsub("\\]","",colnames(All.tax.ref60.nor.com))
library(metamicrobiomeR) 
taxacom.alltax.ref60<-taxa.compare(taxtab=All.tax.ref60.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI","Gender","Antidepression"),longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref60,"ref60.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref60.crude<-taxa.compare(taxtab=All.tax.ref60.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref60.crude,"ref60.metaR.genus.silva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref60, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref60, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref60, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)


###############ref65
###############ref65
###############ref65
###############ref65
###############ref65


kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,5:19),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref65 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref65.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/ref65.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref65.meta),paste("Ref65.",rownames(All.tax.ref65),sep="")) != c(1:dim(All.tax.ref65)[1]))

sum(pmatch(rownames(ref65.meta),paste("Ref65.",rownames(All.tax.ref65),sep="")) != c(1:dim(All.tax.ref65)[1]))
pid <- pmatch(rownames(ref65.meta),paste("Ref65.",rownames(All.tax.ref65),sep=""))
All.tax.ref65 <- All.tax.ref65[pid,]
sum(pmatch(rownames(ref65.meta),paste("Ref65.",rownames(All.tax.ref65),sep="")) != c(1:dim(All.tax.ref65)[1]))
rownames(All.tax.ref65) <- paste("Ref65.",rownames(All.tax.ref65),sep="")

colnames(All.tax.ref65) <- colnames(All.tax.ref65) %>% gsub("d__","k__",.) 

All.tax.ref65.nor.com <- cbind(ref65.meta[,c("Sex","Age","BMI","Depression")],All.tax.ref65)
All.tax.ref65.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref65.nor.com$Depression )
All.tax.ref65.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref65.nor.com$Depression )

colnames(All.tax.ref65.nor.com)[1:4] <- c("Gender","Age","BMI","Depression")
colnames(All.tax.ref65.nor.com) <- gsub(";",".",colnames(All.tax.ref65.nor.com))
colnames(All.tax.ref65.nor.com) <- gsub("-","_",colnames(All.tax.ref65.nor.com))
colnames(All.tax.ref65.nor.com) <- gsub("\\[","",colnames(All.tax.ref65.nor.com))
colnames(All.tax.ref65.nor.com) <- gsub("\\]","",colnames(All.tax.ref65.nor.com))
library(metamicrobiomeR) 
taxacom.alltax.ref65<-taxa.compare(taxtab=All.tax.ref65.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI","Gender"),longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref65,"ref65.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref65.crude<-taxa.compare(taxtab=All.tax.ref65.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref65.crude,"ref65.metaR.genus.silva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref65, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref65, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref65, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)


###############ref77
###############ref77
###############ref77
###############ref77
###############ref77


kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,3:6),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref77 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref77.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/Ref77.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref77.meta),paste("Ref77.",rownames(All.tax.ref77),sep="")) != c(1:dim(All.tax.ref77)[1]))

colnames(All.tax.ref77) <- colnames(All.tax.ref77) %>% gsub("d__","k__",.) 

All.tax.ref77.nor.com <- cbind(ref77.meta[,c("Sex","Age","BMI","Depression")],All.tax.ref77)
All.tax.ref77.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref77.nor.com$Depression )
All.tax.ref77.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref77.nor.com$Depression )

colnames(All.tax.ref77.nor.com)[1:4] <- c("Gender","Age","BMI","Depression")
colnames(All.tax.ref77.nor.com) <- gsub(";",".",colnames(All.tax.ref77.nor.com))
colnames(All.tax.ref77.nor.com) <- gsub("-","_",colnames(All.tax.ref77.nor.com))
colnames(All.tax.ref77.nor.com) <- gsub("\\[","",colnames(All.tax.ref77.nor.com))
colnames(All.tax.ref77.nor.com) <- gsub("\\]","",colnames(All.tax.ref77.nor.com))
library(metamicrobiomeR) 
taxacom.alltax.ref77<-taxa.compare(taxtab=All.tax.ref77.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI","Gender"),longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref77,"ref77.metaR.genus.sliva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxacom.alltax.ref77.crude<-taxa.compare(taxtab=All.tax.ref77.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref77.crude,"ref77.metaR.genus.sliva138.result.crude.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.alltax.ref77, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref77, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref77, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)


###############ref79
###############ref79
###############ref79
###############ref79
###############ref79


kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-1.csv",header = T,check.names = F)
phylum <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-2.csv",header = T,check.names = F)
class <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-3.csv",header = T,check.names = F)
ord <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-4.csv",header = T,check.names = F)
family <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-5.csv",header = T,check.names = F)
genus <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-6.csv",header = T,check.names = F)
species <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-7.csv",header = T,check.names = F)

rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index
rownames(species) <- species$index

length.1 <- dim(kingdom)[2]
colnames(kingdom)
kingdom <- kingdom[,-c(1,3:4),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
length.2 <- dim(kingdom)[2]

meta.col <- length.1-length.2

phylum <- phylum[,-c(1,(dim(phylum)[2]-meta.col+2):dim(phylum)[2]),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
class <- class[,-c(1,(dim(class)[2]-meta.col+2):dim(class)[2]),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
ord <- ord[,-c(1,(dim(ord)[2]-meta.col+2):dim(ord)[2]),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
family <- family[,-c(1,(dim(family)[2]-meta.col+2):dim(family)[2]),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
genus <- genus[,-c(1,(dim(genus)[2]-meta.col+2):dim(genus)[2]),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")
species <- species[,-c(1,(dim(species)[2]-meta.col+2):dim(species)[2]),drop=F]
species.nor <- data.Normalization(species,type = "n10", normalization = "row")


All.tax.ref79 <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor,species.nor)

ref79.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/ref79.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref79.meta),paste("Ref79.",rownames(All.tax.ref79),sep="")) != c(1:dim(All.tax.ref79)[1]))

colnames(All.tax.ref79) <- colnames(All.tax.ref79) %>% gsub("d__","k__",.) 

All.tax.ref79.nor.com <- cbind(ref79.meta[,c("Depression"),drop=F],All.tax.ref79)
All.tax.ref79.nor.com$Depression <- gsub(1,"1_depression",All.tax.ref79.nor.com$Depression )
All.tax.ref79.nor.com$Depression <- gsub(0,"0_depression",All.tax.ref79.nor.com$Depression )

colnames(All.tax.ref79.nor.com) <- gsub(";",".",colnames(All.tax.ref79.nor.com))
colnames(All.tax.ref79.nor.com) <- gsub("-","_",colnames(All.tax.ref79.nor.com))
colnames(All.tax.ref79.nor.com) <- gsub("\\[","",colnames(All.tax.ref79.nor.com))
colnames(All.tax.ref79.nor.com) <- gsub("\\]","",colnames(All.tax.ref79.nor.com))


library(metamicrobiomeR) 
t1 <- grep("s__Ruminococcus_torques",colnames(All.tax.ref79.nor.com))

id <- colnames(All.tax.ref79.nor.com[,t1])[1]
id.ab <- rowSums(All.tax.ref79.nor.com[,t1])
All.tax.ref79.nor.com <- All.tax.ref79.nor.com[,-t1]
All.tax.ref79.nor.com <- cbind(All.tax.ref79.nor.com,id.ab)
colnames(All.tax.ref79.nor.com)[704] <- id
taxacom.alltax.ref79<-taxa.compare(taxtab=All.tax.ref79.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="no",p.adjust.method="fdr")
write.table(taxacom.alltax.ref79,"ref79.metaR.genus.silva138.result.txt",sep="\t",quote = F,row.names = F,col.names = T)


taxcomtab.show(taxcomtab=taxacom.alltax.ref79, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref79, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
taxcomtab.show(taxcomtab=taxacom.alltax.ref79, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)

#####combined result
#####combined result
#####combined result
#####combined result
#####combined result
###########crude model
rm(list = ls())
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/metaR")
library(clusterSim)
library(tidyr)
library(metamicrobiomeR)
source("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/tmp/script_final/taxonomy_name_fill.R")
taxacom.alltax.ref1.crude <- read.table("Ref1.metaR.genus.silva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref23.crude <- read.table("Ref23.metaR.genus.silva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref27 <- read.table("ref27.metaR.genus.silva138.result.txt",sep="\t",head=T)
taxacom.alltax.ref40.crude <- read.table("Ref40.metaR.genus.silva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref60.crude <- read.table("Ref60.metaR.genus.silva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref65.crude <- read.table("Ref65.metaR.genus.silva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref77.crude <- read.table("ref77.metaR.genus.sliva138.result.crude.txt",sep="\t",head=T)
taxacom.alltax.ref79 <- read.table("ref79.metaR.genus.silva138.result.txt",sep="\t",head=T)

rr.id <- c(colnames(taxacom.alltax.ref1.crude),colnames(taxacom.alltax.ref23.crude),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40.crude),colnames(taxacom.alltax.ref60.crude),colnames(taxacom.alltax.ref65.crude),colnames(taxacom.alltax.ref77.crude),colnames(taxacom.alltax.ref79))
rr.id.uniq <-unique(c(colnames(taxacom.alltax.ref1.crude),colnames(taxacom.alltax.ref23.crude),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40.crude),colnames(taxacom.alltax.ref60.crude),colnames(taxacom.alltax.ref65.crude),colnames(taxacom.alltax.ref77.crude),colnames(taxacom.alltax.ref79)))

file.id <- ls(pattern="taxacom.alltax.ref\\d*.crude")
file.id <- c(file.id,"taxacom.alltax.ref27","taxacom.alltax.ref79")
file.id <- file.id[c(1:2,7,3:6,8)]
dat.crude.cb <- c()
for(i in file.id){
    txt <- paste("dat <- ",i,sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat))
    )
    
    dat.crude.cb <- dplyr::bind_rows(dat.crude.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.crude.cb,2,function(x){sum(!is.na(x))})

colnames(dat.crude.cb)[colnames(dat.crude.cb)=="i"]="study"
id_title <- read.table("../ID2Title_meta_analysis_study.txt",sep="\t",head=T)

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.crude.cb$study))
dat.crude.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.crude.cb$pop)
source("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/tmp/script_final/metaR.tax.correct.R")
source("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/tmp/script_final/metaR.meta.niceplot.revised.R")

dat.crude.cb$id <- taxname.correct(dat.crude.cb$id)

meta.depress.crude.allstudy.rr <- meta.taxa(taxcomdat=dat.crude.cb,se.pattern = "Std..Error.",
                                            summary.measure="RR", pool.var="id", studylab="study",
                                            backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")

level6.05.rr<- metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb, tax.lev="l6",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")
level7.05.rr<- metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb, tax.lev="l7",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")

write.table(level6.05.rr,"metaR.genus.allstudy.occ0.5.silva138.p005.gamlss.crude.random.ef.txt",sep="\t",quote=F)
write.table(level7.05.rr,"metaR.species.allstudy.occ0.5.silva138.p005.gamlss.crude.random.ef.txt",sep="\t",quote=F)


###random effect
metadat.rr.crude.allstudy.l2<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.crude.allstudy.l6<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.crude.allstudy.l7<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

unclass.list <- grep("__$",metadat.rr.crude.allstudy.l6$taxsig$id,value = T)
metadat.rr.crude.allstudy.l6$taxsig.all <-metadat.rr.crude.allstudy.l6$taxsig.all[!metadat.rr.crude.allstudy.l6$taxsig.all$id %in% unclass.list,]
metadat.rr.crude.allstudy.l6$taxsig <-metadat.rr.crude.allstudy.l6$taxsig[!metadat.rr.crude.allstudy.l6$taxsig$id %in% unclass.list,]

unclass.list <- grep("__$",metadat.rr.crude.allstudy.l7$taxsig$id,value = T)
metadat.rr.crude.allstudy.l7$taxsig.all <-metadat.rr.crude.allstudy.l7$taxsig.all[!metadat.rr.crude.allstudy.l7$taxsig.all$id %in% unclass.list,]
metadat.rr.crude.allstudy.l7$taxsig <-metadat.rr.crude.allstudy.l7$taxsig[!metadat.rr.crude.allstudy.l7$taxsig$id %in% unclass.list,]

mark.list <- c(as.character(metadat.rr.crude.allstudy.l2$taxsig$id),
               as.character(metadat.rr.crude.allstudy.l6$taxsig$id),
               as.character(metadat.rr.crude.allstudy.l7$taxsig$id))
###remove 65
dat.crude.rm65.cb <- c()
for(i in file.id[-6]){
    txt <- paste("dat <- ",i,sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat))
    )
    
    dat.crude.rm65.cb <- dplyr::bind_rows(dat.crude.rm65.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.crude.rm65.cb,2,function(x){sum(!is.na(x))})

colnames(dat.crude.rm65.cb)[colnames(dat.crude.rm65.cb)=="i"]="study"

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.crude.rm65.cb$study))
dat.crude.rm65.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.crude.rm65.cb$pop)


dat.crude.rm65.cb$id <- taxname.correct(dat.crude.rm65.cb$id)

meta.depress.crude.rm65.rr <- meta.taxa(taxcomdat=dat.crude.rm65.cb,se.pattern = "Std..Error.",
                                        summary.measure="RR", pool.var="id", studylab="study",
                                        backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
level6.05.rr<- metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb, tax.lev="l6",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")
level7.05.rr<- metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb, tax.lev="l7",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")

write.table(level6.05.rr,"metaR.genus.rm65.occ0.5.silva138.p005.gamlss.crude.random.ef.txt",sep="\t",quote=F)
write.table(level7.05.rr,"metaR.species.rm65.occ0.5.silva138.p005.gamlss.crude.random.ef.txt",sep="\t",quote=F)

###random effect
metadat.rr.crude.rm65.l2<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.crude.rm65.l6<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.crude.rm65.l7<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

###add marker in all study
unclass.list <- grep("__$",metadat.rr.crude.rm65.l6$taxsig$id,value = T)
metadat.rr.crude.rm65.l6$taxsig.all <-metadat.rr.crude.rm65.l6$taxsig.all[!metadat.rr.crude.rm65.l6$taxsig.all$id %in% unclass.list,]
metadat.rr.crude.rm65.l6$taxsig <-metadat.rr.crude.rm65.l6$taxsig[!metadat.rr.crude.rm65.l6$taxsig$id %in% unclass.list,]

unclass.list <- grep("__$",metadat.rr.crude.rm65.l7$taxsig$id,value = T)
metadat.rr.crude.rm65.l7$taxsig.all <-metadat.rr.crude.rm65.l7$taxsig.all[!metadat.rr.crude.rm65.l7$taxsig.all$id %in% unclass.list,]
metadat.rr.crude.rm65.l7$taxsig <-metadat.rr.crude.rm65.l7$taxsig[!metadat.rr.crude.rm65.l7$taxsig$id %in% unclass.list,]

mark.rm65.list <- c(as.character(metadat.rr.crude.rm65.l2$taxsig$id),
                    as.character(metadat.rr.crude.rm65.l6$taxsig$id),
                    as.character(metadat.rr.crude.rm65.l7$taxsig$id))

all.mark.list <- unique(c(mark.rm65.list,mark.list))
write.table(all.mark.list,"crude.model.marker.list.txt",quote = F,row.names = F,col.names = F)


##plot in all study
metadat.rr.crude.allstudy.l2<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.allstudy.l6<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.allstudy.l7<-metatab.show(metatab=meta.depress.crude.allstudy.rr$random,com.pooled.tab=dat.crude.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.allstudy.l2$taxsig <- metadat.rr.crude.allstudy.l2$taxsig[metadat.rr.crude.allstudy.l2$taxsig$id %in% all.mark.list,]
metadat.rr.crude.allstudy.l2$taxsig.all <- metadat.rr.crude.allstudy.l2$taxsig.all[metadat.rr.crude.allstudy.l2$taxsig.all$id %in% all.mark.list,]
metadat.rr.crude.allstudy.l6$taxsig.all <- metadat.rr.crude.allstudy.l6$taxsig.all[metadat.rr.crude.allstudy.l6$taxsig.all$id %in% all.mark.list,] 
metadat.rr.crude.allstudy.l6$taxsig <- metadat.rr.crude.allstudy.l6$taxsig[metadat.rr.crude.allstudy.l6$taxsig$id %in% all.mark.list,] 
metadat.rr.crude.allstudy.l7$taxsig.all <- metadat.rr.crude.allstudy.l7$taxsig.all[metadat.rr.crude.allstudy.l7$taxsig.all$id %in% all.mark.list,] 
metadat.rr.crude.allstudy.l7$taxsig <- metadat.rr.crude.allstudy.l7$taxsig[metadat.rr.crude.allstudy.l7$taxsig$id %in% all.mark.list,] 

p.crude.l2 <- meta.niceplot.rev(metadat=metadat.rr.crude.allstudy.l2,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")


p.crude.l6 <- meta.niceplot.rev(metadat=metadat.rr.crude.allstudy.l6,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.crude.l7 <- meta.niceplot.rev(metadat=metadat.rr.crude.allstudy.l7,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

##plot in remove 65 study
metadat.rr.crude.rm65.l2<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.rm65.l6<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.rm65.l7<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.crude.rm65.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.crude.rm65.l2$taxsig <- metadat.rr.crude.rm65.l2$taxsig[metadat.rr.crude.rm65.l2$taxsig$id %in% all.mark.list,]
metadat.rr.crude.rm65.l2$taxsig.all <- metadat.rr.crude.rm65.l2$taxsig.all[metadat.rr.crude.rm65.l2$taxsig.all$id %in% all.mark.list,]
metadat.rr.crude.rm65.l6$taxsig.all <- metadat.rr.crude.rm65.l6$taxsig.all[metadat.rr.crude.rm65.l6$taxsig.all$id %in% all.mark.list,] 
metadat.rr.crude.rm65.l6$taxsig <- metadat.rr.crude.rm65.l6$taxsig[metadat.rr.crude.rm65.l6$taxsig$id %in% all.mark.list,] 
metadat.rr.crude.rm65.l7$taxsig.all <- metadat.rr.crude.rm65.l7$taxsig.all[metadat.rr.crude.rm65.l7$taxsig.all$id %in% all.mark.list,] 
metadat.rr.crude.rm65.l7$taxsig <- metadat.rr.crude.rm65.l7$taxsig[metadat.rr.crude.rm65.l7$taxsig$id %in% all.mark.list,] 

p.crude.rm65.l2 <- meta.niceplot.rev(metadat=metadat.rr.crude.rm65.l2,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.crude.rm65.l6 <- meta.niceplot.rev(metadat=metadat.rr.crude.rm65.l6,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.crude.rm65.l7 <- meta.niceplot.rev(metadat=metadat.rr.crude.rm65.l7,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")


###########adjust model
###########adjust model
###########adjust model
###########adjust model
#rm(list = ls())
taxacom.alltax.ref1 <- read.table("Ref1.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref23 <- read.table("Ref23.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref27 <- read.table("Ref27.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref40 <- read.table("Ref40.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref60 <- read.table("Ref60.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref65 <- read.table("Ref65.metaR.genus.silva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref77 <- read.table("ref77.metaR.genus.sliva138.result.txt",sep="\t",head = T)
taxacom.alltax.ref79 <- read.table("Ref79.metaR.genus.silva138.result.txt",sep="\t",head = T)

rr.id <- c(colnames(taxacom.alltax.ref1),colnames(taxacom.alltax.ref23),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40),colnames(taxacom.alltax.ref60),colnames(taxacom.alltax.ref65),colnames(taxacom.alltax.ref77),colnames(taxacom.alltax.ref79))
rr.id.uniq <-unique(c(colnames(taxacom.alltax.ref1),colnames(taxacom.alltax.ref23),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40),colnames(taxacom.alltax.ref60),colnames(taxacom.alltax.ref65),colnames(taxacom.alltax.ref77),colnames(taxacom.alltax.ref79)))

table(grep("Estimate",c(colnames(taxacom.alltax.ref1),colnames(taxacom.alltax.ref23),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40),colnames(taxacom.alltax.ref60),colnames(taxacom.alltax.ref65),colnames(taxacom.alltax.ref77),colnames(taxacom.alltax.ref79)),value = T))

grep("Estimate.GenderMale",rr.id)
length(c(colnames(taxacom.alltax.ref1),colnames(taxacom.alltax.ref23),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40),colnames(taxacom.alltax.ref60),colnames(taxacom.alltax.ref65),colnames(taxacom.alltax.ref77),colnames(taxacom.alltax.ref79)))

table(grep("Estimate",c(colnames(taxacom.alltax.ref1),colnames(taxacom.alltax.ref23),colnames(taxacom.alltax.ref27),colnames(taxacom.alltax.ref40),colnames(taxacom.alltax.ref60),colnames(taxacom.alltax.ref65),colnames(taxacom.alltax.ref77),colnames(taxacom.alltax.ref79)),value = T))
id <- c("Age","BMI","Gender")
###use all study
file.id <- ls(pattern="taxacom.alltax.ref\\d*$")
dat.adjust.cb <- c()
for(i in file.id){
    txt <- paste("dat <- ",i,sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,
                   grep("Age",colnames(dat)),
                   grep("BMI",colnames(dat)),
                   grep("Gender",colnames(dat)),
                   grep("Depression",colnames(dat))
    )
    
    dat.adjust.cb <- dplyr::bind_rows(dat.adjust.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.adjust.cb,2,function(x){sum(!is.na(x))})
colnames(dat.adjust.cb)[colnames(dat.adjust.cb)=="i"]="study"

id_title <- read.table("../ID2Title_meta_analysis_study.txt",sep="\t",head=T)

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.adjust.cb$study))
dat.adjust.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.adjust.cb$pop)

dat.adjust.cb$id <- taxname.correct(dat.adjust.cb$id)

meta.depress.rr <- meta.taxa(taxcomdat=dat.adjust.cb,se.pattern = "Std..Error.",
                             summary.measure="RR", pool.var="id", studylab="study",
                             backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
meta.depress.allstudy.rr <- meta.depress.rr
level6.05.rr<- metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.cb, tax.lev="l6",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")

write.table(level6.05.rr,"metaR.genus.allstudy.occ0.5.silva138.p005.gamlss.adjust.random.ef.txt",sep="\t",quote=F)

###random effect
metadat.rr.adjust.allstudy.l2<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.adjust.allstudy.l6<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.adjust.allstudy.l7<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

unclass.list <- grep("__$",metadat.rr.adjust.allstudy.l6$taxsig$id,value = T)
metadat.rr.adjust.allstudy.l6$taxsig.all <-metadat.rr.adjust.allstudy.l6$taxsig.all[!metadat.rr.adjust.allstudy.l6$taxsig.all$id %in% unclass.list,]
metadat.rr.adjust.allstudy.l6$taxsig <-metadat.rr.adjust.allstudy.l6$taxsig[!metadat.rr.adjust.allstudy.l6$taxsig$id %in% unclass.list,]

unclass.list <- grep("__$",metadat.rr.adjust.allstudy.l7$taxsig$id,value = T)
metadat.rr.adjust.allstudy.l7$taxsig.all <-metadat.rr.adjust.allstudy.l7$taxsig.all[!metadat.rr.adjust.allstudy.l7$taxsig.all$id %in% unclass.list,]
metadat.rr.adjust.allstudy.l7$taxsig <-metadat.rr.adjust.allstudy.l7$taxsig[!metadat.rr.adjust.allstudy.l7$taxsig$id %in% unclass.list,]

mark.list <- c(as.character(metadat.rr.adjust.allstudy.l2$taxsig$id),
               as.character(metadat.rr.adjust.allstudy.l6$taxsig$id),
               as.character(metadat.rr.adjust.allstudy.l7$taxsig$id))

all.list <- unique(c(mark.list,all.mark.list))


###################remove ref65
###################remove ref65
dat.adjust.rm65.cb <- c()
for(i in file.id[-6]){
    txt <- paste("dat <- ",i,sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,
                   grep("Age",colnames(dat)),
                   grep("BMI",colnames(dat)),
                   grep("Gender",colnames(dat)),
                   grep("Depression",colnames(dat))
    )
    
    dat.adjust.rm65.cb <- dplyr::bind_rows(dat.adjust.rm65.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.adjust.rm65.cb,2,function(x){sum(!is.na(x))})
colnames(dat.adjust.rm65.cb)[colnames(dat.adjust.rm65.cb)=="i"]="study"

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.adjust.rm65.cb$study))
dat.adjust.rm65.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.adjust.rm65.cb$pop)
dat.adjust.rm65.cb$id <- taxname.correct(dat.adjust.rm65.cb$id)


meta.depress.rr <- meta.taxa(taxcomdat=dat.adjust.rm65.cb,se.pattern = "Std..Error.",
                             summary.measure="RR", pool.var="id", studylab="study",
                             backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")
meta.depress.rm65.rr <- meta.depress.rr

level6.05.rr<- metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.rm65.cb, tax.lev="l6",se.pattern="Std..Error.",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="table",highest.lev = "s")

metadat.rr<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,display="data",highest.lev = "s")

write.table(level6.05.rr,"metaR.genus.rm65.occ0.5.silva138.p005.gamlss.adjust.random.ef.txt",sep="\t",quote=F)

###random effect
metadat.rr.adjust.rm65.l2<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.adjust.rm65.l6<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

metadat.rr.adjust.rm65.l7<-metatab.show(metatab=meta.depress.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=0.05,highest.lev = "s",display="data")

unclass.list <- grep("__$",metadat.rr.adjust.rm65.l6$taxsig$id,value = T)
metadat.rr.adjust.rm65.l6$taxsig.all <-metadat.rr.adjust.rm65.l6$taxsig.all[!metadat.rr.adjust.rm65.l6$taxsig.all$id %in% unclass.list,]
metadat.rr.adjust.rm65.l6$taxsig <-metadat.rr.adjust.rm65.l6$taxsig[!metadat.rr.adjust.rm65.l6$taxsig$id %in% unclass.list,]


unclass.list <- grep("__$",metadat.rr.adjust.rm65.l7$taxsig$id,value = T)
metadat.rr.adjust.rm65.l7$taxsig.all <-metadat.rr.adjust.rm65.l7$taxsig.all[!metadat.rr.adjust.rm65.l7$taxsig.all$id %in% unclass.list,]
metadat.rr.adjust.rm65.l7$taxsig <-metadat.rr.adjust.rm65.l7$taxsig[!metadat.rr.adjust.rm65.l7$taxsig$id %in% unclass.list,]

mark.rm65.list <- c(as.character(metadat.rr.adjust.rm65.l2$taxsig$id),
                    as.character(metadat.rr.adjust.rm65.l6$taxsig$id),
                    as.character(metadat.rr.adjust.rm65.l7$taxsig$id))

all.list <- unique(c(mark.rm65.list,all.list ))

write.table(all.list,"crude.and.adjust.model.marker.list.txt",quote = F)

###plot in all study

metadat.rr.adjust.allstudy.l2<-metatab.show(metatab=meta.depress.allstudy.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.adjust.allstudy.l6<-metatab.show(metatab=meta.depress.allstudy.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.adjust.allstudy.l7<-metatab.show(metatab=meta.depress.allstudy.rr$random,com.pooled.tab=dat.adjust.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.adjust.allstudy.l2$taxsig <- metadat.rr.adjust.allstudy.l2$taxsig[metadat.rr.adjust.allstudy.l2$taxsig$id %in% all.list,]
metadat.rr.adjust.allstudy.l2$taxsig.all <- metadat.rr.adjust.allstudy.l2$taxsig.all[metadat.rr.adjust.allstudy.l2$taxsig.all$id %in% all.list,]
metadat.rr.adjust.allstudy.l6$taxsig.all <- metadat.rr.adjust.allstudy.l6$taxsig.all[metadat.rr.adjust.allstudy.l6$taxsig.all$id %in% all.list,] 
metadat.rr.adjust.allstudy.l6$taxsig <- metadat.rr.adjust.allstudy.l6$taxsig[metadat.rr.adjust.allstudy.l6$taxsig$id %in% all.list,] 
metadat.rr.adjust.allstudy.l7$taxsig.all <- metadat.rr.adjust.allstudy.l7$taxsig.all[metadat.rr.adjust.allstudy.l7$taxsig.all$id %in% all.list,] 
metadat.rr.adjust.allstudy.l7$taxsig <- metadat.rr.adjust.allstudy.l7$taxsig[metadat.rr.adjust.allstudy.l7$taxsig$id %in% all.list,] 

p.adjust.l2 <- meta.niceplot.rev(metadat=metadat.rr.adjust.allstudy.l2,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.adjust.l6 <- meta.niceplot.rev(metadat=metadat.rr.adjust.allstudy.l6,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.adjust.l7 <- meta.niceplot.rev(metadat=metadat.rr.adjust.allstudy.l7,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

###plot in remove 65 study

metadat.rr.adjust.rm65.l2<-metatab.show(metatab=meta.depress.rm65.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l2",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.adjust.rm65.l6<-metatab.show(metatab=meta.depress.rm65.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l6",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")

metadat.rr.adjust.rm65.l7<-metatab.show(metatab=meta.depress.rm65.rr$random,com.pooled.tab=dat.adjust.rm65.cb,se.pattern="Std..Error.", tax.lev="l7",showvar="Depression1_depression",p.cutoff.type="p", p.cutoff=1,highest.lev = "s",display="data")


metadat.rr.adjust.rm65.l2$taxsig <- metadat.rr.adjust.rm65.l2$taxsig[metadat.rr.adjust.rm65.l2$taxsig$id %in% all.list,]
metadat.rr.adjust.rm65.l2$taxsig.all <- metadat.rr.adjust.rm65.l2$taxsig.all[metadat.rr.adjust.rm65.l2$taxsig.all$id %in% all.list,]
metadat.rr.adjust.rm65.l6$taxsig.all <- metadat.rr.adjust.rm65.l6$taxsig.all[metadat.rr.adjust.rm65.l6$taxsig.all$id %in% all.list,] 
metadat.rr.adjust.rm65.l6$taxsig <- metadat.rr.adjust.rm65.l6$taxsig[metadat.rr.adjust.rm65.l6$taxsig$id %in% all.list,] 
metadat.rr.adjust.rm65.l7$taxsig.all <- metadat.rr.adjust.rm65.l7$taxsig.all[metadat.rr.adjust.rm65.l7$taxsig.all$id %in% all.list,] 
metadat.rr.adjust.rm65.l7$taxsig <- metadat.rr.adjust.rm65.l7$taxsig[metadat.rr.adjust.rm65.l7$taxsig$id %in% all.list,] 

p.adjust.rm65.l2 <- meta.niceplot.rev(metadat=metadat.rr.adjust.rm65.l2,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.adjust.rm65.l6 <- meta.niceplot.rev(metadat=metadat.rr.adjust.rm65.l6,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

p.adjust.rm65.l7 <- meta.niceplot.rev(metadat=metadat.rr.adjust.rm65.l7,sumtype="taxa",level="sub",p="p", p.adjust="p.adjust",phyla.col="rainbow",p.sig.heat="yes", heat.forest.width.ratio =c(1.5,1), leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10,  forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

pdf("metaR.allstudy.pgs.cb.pdf",width = 30,height = 12,onefile=FALSE)
ggarrange(p.crude.l2$heatmap,p.crude.l2$forest,p.crude.rm65.l2$forest,
          p.crude.l6$heatmap,p.crude.l6$forest,p.crude.rm65.l6$forest,
          p.crude.l7$heatmap,p.crude.l7$forest,p.crude.rm65.l7$forest,
          
          p.adjust.l2$heatmap,p.adjust.l2$forest,p.adjust.rm65.l2$forest,
          p.adjust.l6$heatmap,p.adjust.l6$forest,p.adjust.rm65.l6$forest,
          p.adjust.l7$heatmap,p.adjust.l7$forest,p.adjust.rm65.l7$forest,
          ncol = 3,nrow = 6,align = "hv",heights=c(1.3,2,6,1.4,4,6),widths=c(1.8,2,2),legend="right",common.legend=T,legend.grob = get_legend(p.adjust.rm65.l7))
dev.off()
