###Review 65 dada2
rm(list=ls())
library(ggplot2)
library(car)
library(gridExtra)
library(grid)
library(ggpubr)
library(lattice)
library("kableExtra")
library(dplyr)
library(MASS)
library(RVAideMemoire)
library(clusterSim)
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/")

###marker in systematic review
decrease <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/decrease.maker.list.txt",sep="\t")
increase <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/increase.maker.list.txt",sep="\t")
##metadata

meta1 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/DataDownload/Ref65/Cohort1_mappingfile.txt",sep="\t",head=T,comment.char = "")
meta2 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/DataDownload/Ref65/map_01.11.18.txt",sep="\t",row.names = 1,head=T,comment.char = "")
meta1[1:4,1:7]
meta2[1:4,1:7]
commom.sample <- intersect(rownames(meta2),meta1$X.SampleID)
metadata.all <- cbind(meta2[pmatch(commom.sample,rownames(meta2)),],
                  meta1[pmatch(commom.sample,meta1$X.SampleID),])

pmatch(rownames(metadata.all),metadata.all$X.SampleID)==c(1:513)

metadata.duplicat <- unique(metadata.all[,-c(13,14)])
colnames(metadata.all)
metadata <- metadata.duplicat[,c(1,5:8,14,23,25,41,39,20,22)]
metadata$depression <- 0
metadata$depression[metadata$Depression_overall == "Higher_Than_7"] <- 1

meta.Colon_Left <- metadata[metadata$AnatomicLocation=="Colon_Left",]
table(meta.Colon_Left$depression)
dup.ind <- names(which(table(meta.Colon_Left$patient_number)>1))
dup.sample <-rownames(meta.Colon_Left[meta.Colon_Left$patient_number %in% dup.ind ,])
dup.sample.recod <- metadata.all[rownames(metadata.all) %in% dup.sample,]

meta.Colon_Right <- metadata[metadata$AnatomicLocation=="Colon_Right",]
table(meta.Colon_Right$depression)

meta.Colon_Transverse <- metadata[metadata$AnatomicLocation=="Colon_Transverse",]
table(meta.Colon_Transverse$depression)

meta.Ileum <- metadata[metadata$AnatomicLocation=="Ileum",]
table(meta.Ileum$depression)

meta.Rectum <- metadata[metadata$AnatomicLocation=="Rectum",]
table(meta.Rectum$depression)
######################
######################
##1.alpha diversity

#observe otu
obt <- read.table("dada2/observed_otus_vector.qza/35bff01d-d3df-46fe-8f41-1ebdf33a5a9a/data/alpha-diversity.tsv")

#chao1
chao1.d <- read.table("dada2/chao1_alpha_diverisyt.qza/84ef786f-1900-4b22-a9d4-bab6c2af2a66/data/alpha-diversity.tsv")

#shannon
shannon <- read.table("dada2/shannon_vector.qza/118a0ccf-4c62-4a44-83dd-aa7eeb3540b7/data/alpha-diversity.tsv")

#simpson
simpson <- read.table("dada2/simpson_alpha_diverisyt.qza/2562d1ab-969f-4e7b-a3ab-58f154ad9647/data/alpha-diversity.tsv")

##2.taxonomy
kingdom <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-1.csv")
phylum <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-2.csv")
class <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-3.csv")
ord <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-4.csv")
family <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-5.csv")
genus <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-6.csv")
species <- read.csv("dada2/taxarare-bar-plots.qzv/ca080795-8899-46f5-9233-d54c1c55a7ac/data/level-7.csv")
rownames(kingdom) <- kingdom$index
rownames(phylum) <- phylum$index
rownames(class) <- class$index
rownames(ord) <- ord$index
rownames(family) <- family$index
rownames(genus) <- genus$index

###change tax name
#phylum
py.id <- grep("\\.__",colnames(phylum))
colnames(phylum)[py.id] <- gsub("\\.__",".D_1__",colnames(phylum)[py.id])

#class
cl.id1 <- grep("\\.__\\.__",colnames(class))
colnames(class)[cl.id1] <- gsub("\\.__\\.__",".D_1__.D_2__",colnames(class)[cl.id1])
cl.id2 <- grep("\\.__",colnames(class))
colnames(class)[cl.id2] <- gsub("\\.__",".D_2__",colnames(class)[cl.id2])

#order
od.id1 <- grep("\\.__\\.__\\.__",colnames(ord))
colnames(ord)[od.id1] <- gsub("\\.__\\.__\\.__",".D_1__.D_2__.D_3__",colnames(ord)[od.id1])

od.id2 <- grep("\\.__\\.__",colnames(ord))
colnames(ord)[od.id2] <- gsub("\\.__\\.__",".D_2__.D_3__",colnames(ord)[od.id2])

od.id3 <- grep("\\.__",colnames(ord))
colnames(ord)[od.id3] <- gsub("\\.__",".D_3__",colnames(ord)[od.id3])

#family
family.id0 <- grep("\\.__\\.__\\.__\\.__",colnames(family))
colnames(family)[family.id0] <- gsub("\\.__\\.__\\.__\\.__",".D_1__.D_2__.D_3__.D_4__",colnames(family)[family.id0])

family.id1 <- grep("\\.__\\.__\\.__",colnames(family))
colnames(family)[family.id1] <- gsub("\\.__\\.__\\.__",".D_2__.D_3__.D_4__",colnames(family)[family.id1])

family.id2 <- grep("\\.__\\.__",colnames(family))
colnames(family)[family.id2] <- gsub("\\.__\\.__",".D_3__.D_4__",colnames(family)[family.id2])

family.id3 <- grep("\\.__",colnames(family))
colnames(family)[family.id3] <- gsub("\\.__",".D_4__",colnames(family)[family.id3])

#genus
genus.id0 <- grep("\\.__\\.__\\.__\\.__\\.__",colnames(genus))
colnames(genus)[genus.id0] <- gsub("\\.__\\.__\\.__\\.__\\.__",".D_1__.D_2__.D_3__.D_4__.D_5__",colnames(genus)[genus.id0])

genus.id0 <- grep("\\.__\\.__\\.__\\.__",colnames(genus))
colnames(genus)[genus.id0] <- gsub("\\.__\\.__\\.__\\.__",".D_2__.D_3__.D_4__.D_5__",colnames(genus)[genus.id0])

genus.id1 <- grep("\\.__\\.__\\.__",colnames(genus))
colnames(genus)[genus.id1] <- gsub("\\.__\\.__\\.__",".D_3__.D_4__.D_5__",colnames(genus)[genus.id1])

genus.id2 <- grep("\\.__\\.__",colnames(genus))
colnames(genus)[genus.id2] <- gsub("\\.__\\.__",".D_4__.D_5__",colnames(genus)[genus.id2])

genus.id3 <- grep("\\.__",colnames(genus))
colnames(genus)[genus.id3] <- gsub("\\.__",".D_5__",colnames(genus)[genus.id3])


#species
species.id0 <- grep("\\.__\\.__\\.__\\.__\\.__\\.__",colnames(species))
colnames(species)[species.id0] <- gsub("\\.__\\.__\\.__\\.__\\.__\\.__",".D_1__.D_2__.D_3__.D_4__.D_5__.D_6__",colnames(species)[species.id0])

species.id0 <- grep("\\.__\\.__\\.__\\.__\\.__",colnames(species))
colnames(species)[species.id0] <- gsub("\\.__\\.__\\.__\\.__\\.__","D_2__.D_3__.D_4__.D_5__.D_6__",colnames(species)[species.id0])

species.id0 <- grep("\\.__\\.__\\.__\\.__",colnames(species))
colnames(species)[species.id0] <- gsub("\\.__\\.__\\.__\\.__",".D_3__.D_4__.D_5__.D_6__",colnames(species)[species.id0])

species.id1 <- grep("\\.__\\.__\\.__",colnames(species))
colnames(species)[species.id1] <- gsub("\\.__\\.__\\.__",".D_4__.D_5__.D_6__",colnames(species)[species.id1])

species.id2 <- grep("\\.__\\.__",colnames(species))
colnames(species)[species.id2] <- gsub("\\.__\\.__",".D_5__.D_6__",colnames(species)[species.id2])

species.id3 <- grep("\\.__",colnames(species))
colnames(species)[species.id3] <- gsub("\\.__",".D_6__",colnames(species)[species.id3])

colnames(kingdom)
kingdom <- kingdom[,-c(1,3:18),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
colnames(phylum)
phylum <- phylum[,-c(1,28:42),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
colnames(class)
class <- class[,-c(1,54:68),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
colnames(ord)
ord <- ord[,-c(1,110:124),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
colnames(family)
family <- family[,-c(1,203:217),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
colnames(genus)
genus <- genus[,-c(1,487:501),drop=F]
genus.nor <- data.Normalization(genus,type = "n10", normalization = "row")

All.tax <- cbind.data.frame(kingdom.nor,phylum.nor,class.nor,ord.nor,family.nor,genus.nor)

length(colnames(All.tax)) ==
    length(unique(colnames(All.tax)))

All.tax[,colnames(All.tax)[duplicated(colnames(All.tax))]]
#####reorder
all.equal(rownames(simpson),rownames(obt))
all.equal(rownames(simpson),rownames(shannon))
all.equal(rownames(simpson),rownames(chao1.d))
all.equal(rownames(simpson),rownames(kingdom))
all.equal(rownames(simpson),rownames(phylum))
all.equal(rownames(simpson),rownames(class))
all.equal(rownames(simpson),rownames(ord))
all.equal(rownames(simpson),rownames(family))
all.equal(rownames(simpson),rownames(genus))

alpha <- cbind.data.frame(obt,chao1.d,simpson,shannon)

common.id <- intersect(rownames(alpha),rownames(metadata))

metadata <- metadata[pmatch(common.id,rownames(metadata)),]
alpha <- alpha[pmatch(common.id,rownames(alpha)),]

alpha <- cbind.data.frame(alpha,metadata)
colnames(alpha)
write.table(alpha,"alpha_diverisy.rawdata.txt",sep="\t",quote = F,col.names = T,row.names = T)
###alpha diversity distribution
a1 <- ggplot(alpha, aes(x=observed_otus)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")
a2 <- ggplot(alpha, aes(x=chao1)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")
a3 <- ggplot(alpha, aes(x=simpson)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")
a4 <- ggplot(alpha, aes(x=shannon)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")

shapiro.test(alpha[,1])$p.value
shapiro.test(alpha[,2])$p.value
shapiro.test(alpha[,3])$p.value
shapiro.test(alpha[,4])$p.value

grid.arrange(a1,a2,a3,a4,ncol=2,nrow=2)

#####alpha diversity test

colnames(alpha)
loc <- names(table(alpha$AnatomicLocation))
for(j in loc){
    wil.res <- c()
    t.res <- c()
    dat.tmp <- alpha[alpha$AnatomicLocation == j,]
        
        for(i in c(1:4)){
            p.t <- t.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
            p.v <- wilcox.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
            means <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],mean))
            stand.d <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],sd))
            t.res<-rbind(t.res,c(p.t,means,stand.d))
            iqr <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],quantile))
            wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
        }
        
        rownames(wil.res) <- rownames(t.res) <- colnames(alpha)[c(1:4)]
        colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
        colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")
        
        wil.res <- cbind(wil.res,dim(dat.tmp)[1],table(dat.tmp$depression)[1],table(dat.tmp$depression)[2])
        colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")
        
        t.res <- cbind(t.res,dim(dat.tmp)[1],table(dat.tmp$depression)[1],table(dat.tmp$depression)[2])
        colnames(t.res)[6:8] <- c("Total sample size","non-depression","depression")
        
        txt <- paste("wil.res.",j,"<- wil.res",sep="")
        eval((parse(text =txt)))
        txt <- paste("t.res",j,"<- t.res",sep="")
        eval((parse(text =txt)))
        
        txt <- paste(
        "write.table(wil.res.",j,",\"Ref65_dada2_alpha.diversity.wilcoxon.",j,".test.result.txt\",sep=\"\t\",quote = F,col.names = T,row.names = T)",sep = ""
        )
        eval((parse(text =txt)))
        
        txt <- paste(
            "write.table(t.res",j,",\"Ref65_dada2_alpha.diversity.t.test.",j,".result.txt\",sep=\"\t\",quote = F,col.names = T,row.names = T)",sep = ""
        )
        eval((parse(text =txt)))
}

library("kableExtra")

wil.res.Colon_Left[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (Colon_Left)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Right[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (Colon_Right)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Transverse[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (Colon_Transverse)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Ileum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (Ileum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Rectum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (Rectum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 


###########alpha divesity odd ratio
alpha.odd <- alpha[,c(1:4,13,14,12,5,11,15:17)]
colnames(alpha.odd)[5:12]<-c("Age","BMI","Sex","patient_ID","location","Disease","DiseaseActivity","Depression")
alpha.odd$DiseaseActivity <- as.numeric(factor(alpha.odd$DiseaseActivity,levels = c("No_Disease","Mild","Moderate","Severe")))
alpha.odd$Disease.activity <- 0
alpha.odd$Disease.activity[alpha.odd$DiseaseActivity >1] <- 1

odd.crude.res <- c()
odd.gee.res <- c()
odd.gee.loc.res <- c()
odd.gee.disease.res <- c()
library(survival)
library(geepack)
library(doBy)
for(i in 1:4){
    #crude
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    reg <- glm(Depression ~alpha, data=data, family=binomial)
    odd.crude.res <- rbind(odd.crude.res,
                           round(exp(cbind(coef(reg), confint(reg))),3)[2,])
    
    #adjust
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c("Depression","Age","BMI","Sex","location","patient_ID","Disease","DiseaseActivity","Disease.activity")])
    colnames(data)[1] <- "alpha"
    gee.adjust.location <- geeglm(Depression ~alpha+Age+BMI+Sex+location+Disease+DiseaseActivity, family=binomial, data=data,id=patient_ID, corstr = "independence")
    odd.gee.loc.res <- rbind(odd.gee.loc.res,
                         c(exp(c(esticon(gee.adjust.location, c(0,1,0,0,0,0,0,0,0,0,0))$estimate,esticon(gee.adjust.location,c(0,1,0,0,0,0,0,0,0,0,0))$lwr,esticon(gee.adjust.location, c(0,1,0,0,0,0,0,0,0,0,0))$upr)),paste("raw",colnames(alpha.odd)[i],sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    
    gee.adjust.disease <- geeglm(Depression ~alpha+Age+BMI+Sex+location+Disease.activity, family=binomial, data=data,id=patient_ID, corstr = "independence")
    odd.gee.disease.res <- rbind(odd.gee.disease.res,
                             c(exp(c(esticon(gee.adjust.disease, c(0,1,0,0,0,0,0,0,0,0))$estimate,esticon(gee.adjust.disease,c(0,1,0,0,0,0,0,0,0,0))$lwr,esticon(gee.adjust.disease, c(0,1,0,0,0,0,0,0,0,0))$upr)),paste("raw",colnames(alpha.odd)[i],sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    
    lev.loc <- unique(data$location)
    for(j in lev.loc){
        dat.e <- data[data$location == j,]
        colnames(dat.e)[1] <- "alpha"
        gee.adjust <- geeglm(Depression ~alpha+Age+BMI+Sex+Disease+DiseaseActivity, family=binomial, data=dat.e,  id=patient_ID, corstr = "independence")
        
        odd.gee.res <- rbind(odd.gee.res,
                             c(exp(c(esticon(gee.adjust, c(0,1,0,0,0,0,0))$estimate,esticon(gee.adjust,c(0,1,0,0,0,0,0))$lwr,esticon(gee.adjust, c(0,1,0,0,0,0,0))$upr)),paste("raw",colnames(alpha.odd)[i],j,sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    }
    
    
    ###seperated by median
    #crude
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg <- glm(Depression ~H.alpha, data=data, family=binomial)
    odd.crude.res <- rbind(odd.crude.res,
                           round(exp(cbind(coef(reg), confint(reg))),3)[2,])
    
    #adjust
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c("Depression","Age","BMI","Sex","location","patient_ID","Disease","DiseaseActivity","Disease.activity")])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    gee.adjust.location <- geeglm(Depression ~H.alpha+Age+BMI+Sex+location+Disease+DiseaseActivity, family=binomial, data=data,id=patient_ID, corstr = "independence")
    odd.gee.loc.res <- rbind(odd.gee.loc.res,
                             c(exp(c(esticon(gee.adjust.location, c(0,1,0,0,0,0,0,0,0,0,0))$estimate,esticon(gee.adjust.location,c(0,1,0,0,0,0,0,0,0,0,0))$lwr,esticon(gee.adjust.location, c(0,1,0,0,0,0,0,0,0,0,0))$upr)),paste("Higher",colnames(alpha.odd)[i],sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    
    gee.adjust.disease <- geeglm(Depression ~H.alpha+Age+BMI+Sex+location+Disease.activity, family=binomial, data=data,id=patient_ID, corstr = "independence")
    odd.gee.disease.res <- rbind(odd.gee.disease.res,
                                 c(exp(c(esticon(gee.adjust.disease, c(0,1,0,0,0,0,0,0,0,0))$estimate,esticon(gee.adjust.disease,c(0,1,0,0,0,0,0,0,0,0))$lwr,esticon(gee.adjust.disease, c(0,1,0,0,0,0,0,0,0,0))$upr)),paste("Higher",colnames(alpha.odd)[i],sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    
    
    lev.loc <- unique(data$location)
    for(j in lev.loc){
        dat.e <- data[data$location == j,]
        colnames(dat.e)[1] <- "alpha"
        med.alpha <- median(dat.e$alpha)
        dat.e$H.alpha <- 1
        dat.e$H.alpha[dat.e$alpha < med.alpha] <-0
        gee.adjust <- geeglm(Depression ~H.alpha+Age+BMI+Sex+Disease+DiseaseActivity, family=binomial, data=dat.e,  id=patient_ID, corstr = "independence")
        
        odd.gee.res <- rbind(odd.gee.res,
                             c(exp(c(esticon(gee.adjust, c(0,1,0,0,0,0,0))$estimate,esticon(gee.adjust,c(0,1,0,0,0,0,0))$lwr,esticon(gee.adjust, c(0,1,0,0,0,0,0))$upr)),paste("Higher",colnames(alpha.odd)[i],j,sep="_"),dim(data)[1],table(data$Depression)[1],table(data$Depression)[2]))
    }
}
rownames(odd.crude.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.gee.res) <- odd.gee.res[,4]
odd.gee.res<-odd.gee.res[,-4]
odd.crude.res
odd.crude.res <- cbind(odd.crude.res,dim(alpha)[1],table(alpha$depression)[1],table(alpha$depression)[2])
colnames(odd.crude.res)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.gee.res
colnames(odd.gee.res)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.gee.loc.res
rownames(odd.gee.loc.res) <- rownames(odd.gee.disease.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))
odd.gee.loc.res <- odd.gee.loc.res[,-4]
odd.gee.disease.res <- odd.gee.disease.res[,-4]
colnames(odd.gee.loc.res)[4:6] <- colnames(odd.gee.disease.res)[4:6] <- c("Total_sample_size","non_depression","depression")


write.table(odd.crude.res,"Ref65_dada_alpha.diversity.oR.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.res,"Ref65_dada_alpha.diversity.OR_GEE.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.loc.res,"Ref65_dada_alpha.diversity.OR_GEE.adjust.location.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.disease.res,"Ref65_dada_alpha.diversity.OR_GEE.adjust.location.disease.result.txt",sep="\t",quote = F,col.names = T,row.names = T)


###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum distribution
phylum.d <- cbind.data.frame(phylum.nor[pmatch(common.id,rownames(phylum.nor)),],metadata)
sort(apply(phylum.d[,c(1:26)],2,function(x){sum(x!=0)}),decreasing = T)
phy.p <- colnames(phylum.d)[
    which(
        apply(phylum.d[,c(1:26)],2,function(x){sum(x!=0)})/dim(phylum.d)[1] > 0.2
    )
]
phy.p
p1 <- ggplot(phylum.d, aes(x=D_0__Bacteria.D_1__Actinobacteria)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")
p1
p2 <- ggplot(phylum.d, aes(x=D_0__Bacteria.D_1__Bacteroidetes)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")

p3 <- ggplot(phylum.d, aes(x=D_0__Bacteria.D_1__Firmicutes)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")

p4 <- ggplot(phylum.d, aes(x=D_0__Bacteria.D_1__Proteobacteria)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")

p5 <- ggplot(phylum.d, aes(x=D_0__Bacteria.D_1__Cyanobacteria)) +
    geom_histogram(fill="black", alpha=0.5, position="identity")

grid.arrange(p1,p2,p3,p4,p5,ncol=2,nrow=3)

##normal distribution
nor.res <- c()
for(i in phy.p){
    nor.res <- c( nor.res,
                  shapiro.test(phylum.d[,i])$p.value
    )
}
names(nor.res) <- phy.p
nor.res[nor.res>0.05]
###all phylum are significant different with normal distribution

#####phylum test by wiltest and t.test
colnames(phylum.d)
loc <- names(table(phylum.d$AnatomicLocation))
for(j in loc){
    wil.res <- c()
    t.res <- c()
    dat.tmp <- phylum.d[phylum.d$AnatomicLocation == j,]
    
    for(i in c(1:26)){
        p.t <- t.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
        p.v <- wilcox.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
        means <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],mean))
        stand.d <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],sd))
        t.res<-rbind(t.res,c(p.t,means,stand.d))
        iqr <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],quantile))
        wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    }
    
    rownames(wil.res) <- rownames(t.res) <- colnames(phylum.d)[c(1:26)]
    colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
    colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")
    
    wil.res <- cbind(wil.res,dim(dat.tmp)[1],table(dat.tmp$depression)[1],table(dat.tmp$depression)[2])
    colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")
    
    txt <- paste("wil.res.",j,"<- wil.res",sep="")
    eval((parse(text =txt)))
    txt <- paste("t.res",j,"<- t.res",sep="")
    eval((parse(text =txt)))
    
    txt <- paste(
        "write.table(wil.res.",j,",\"Ref65_dada2_phylum.d.diversity.wilcoxon.",j,".test.result.txt\",sep=\"\t\",quote = F,col.names = T,row.names = T)",sep = ""
    )
    eval((parse(text =txt)))
}

library("kableExtra")

wil.res.Colon_Left[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%  filter(P.value <0.05) %>%
    kbl(caption = "phylum.d diversity (Colon_Left)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Right[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%  filter(P.value <0.05) %>%
    kbl(caption = "phylum.d diversity (Colon_Right)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Transverse[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "phylum.d diversity (Colon_Transverse)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Ileum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "phylum.d diversity (Ileum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Rectum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "phylum.d diversity (Rectum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 


###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus distribution
genus.d <- cbind.data.frame(genus.nor[pmatch(common.id,rownames(genus.nor)),],metadata)
colnames(genus.d)
sum(apply(genus.d,2,function(x){sum(x==0)})==465)

#####genus test
colnames(genus.d)
loc <- names(table(genus.d$AnatomicLocation))
for(j in loc){
    wil.res <- c()
    t.res <- c()
    dat.tmp <- genus.d[genus.d$AnatomicLocation == j,]
    
    for(i in c(1:485)){
        p.t <- t.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
        p.v <- wilcox.test(dat.tmp[,i] ~ dat.tmp[,"depression"])$p.value
        means <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],mean))
        stand.d <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],sd))
        t.res<-rbind(t.res,c(p.t,means,stand.d))
        iqr <- as.matrix(by(dat.tmp[,i],dat.tmp[,"depression"],quantile))
        wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    }
    
    rownames(wil.res) <- rownames(t.res) <- colnames(genus.d)[c(1:485)]
    colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
    colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")
    
    wil.res <- cbind(wil.res,dim(dat.tmp)[1],table(dat.tmp$depression)[1],table(dat.tmp$depression)[2])
    colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")
    
    txt <- paste("wil.res.",j,"<- wil.res",sep="")
    eval((parse(text =txt)))
    txt <- paste("t.res",j,"<- t.res",sep="")
    eval((parse(text =txt)))
    
    txt <- paste(
        "write.table(wil.res.",j,",\"Ref65_dada2_genus.d.diversity.wilcoxon.",j,".test.result.txt\",sep=\"\t\",quote = F,col.names = T,row.names = T)",sep = ""
    )
    eval((parse(text =txt)))
}

library("kableExtra")

wil.res.Colon_Left[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "genus.d diversity (Colon_Left)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Right[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "genus.d diversity (Colon_Right)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Colon_Transverse[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "genus.d diversity (Colon_Transverse)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Ileum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "genus.d diversity (Ileum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.Rectum[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>% filter(P.value <0.05) %>%
    kbl(caption = "genus.d diversity (Rectum)") %>%
    kable_classic(full_width = F, html_font = "Cambria") 


genus.m <- strsplit(id1,split='D_5__')
genus.m.g <-  sapply(genus.m, function(x){x[2]})
genus.m.g <- cbind(genus.m.g,"decrease")
genus.m.g[mark.1[,2]<mark.1[,3],2] <- "increase"

maker <- rbind(genus.m.g)

wil.res <- c()
t.res <- c()
for(i in 1:485){
    ##test for normal distribution.
    ##From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
  
    p.t <- t.test(genus.d[,i] ~ genus.d[,"depression"])$p.value
    p.v <- wilcox.test(genus.d[,i] ~ genus.d[,"depression"])$p.value
    means <- as.matrix(by(genus.d[,i],genus.d[,"depression"],mean))
    stand.d <- as.matrix(by(genus.d[,i],genus.d[,"depression"],sd))
    t.res<-rbind(t.res,c(p.normal,p.t,means,stand.d))
    iqr <- as.matrix(by(genus.d[,i],genus.d[,"depression"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(genus.d)[1:485]

colnames(t.res) <-  c("P.value_normal","P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")


colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(4) %>% filter(P.value <0.05) %>%
    kbl(caption = "Genus marker") %>%
    kable_classic(full_width = F, html_font = "Cambria") 


wil.res <- cbind(wil.res,dim(genus.d)[1],table(genus.d$depression)[1],table(genus.d$depression)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")

write.table(t.res,"Ref65_dada2_genus.d.diversity.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref65_dada2_genus.d.diversity.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
##consist
pmatch(increase$Var1,maker[maker[,2]=="increase",1])
pmatch(decrease$Var1,maker[maker[,2]=="decrease",1])
##conflict
pmatch(increase$Var1,maker[maker[,2]=="decrease",1])
pmatch(decrease$Var1,maker[maker[,2]=="increase",1])
####output lefseq
genus.d <- cbind.data.frame(genus,metadata)
colnames(genus.d)[1:10]
lefse.genus <- t(genus.d[,c(496,1:485)])
rownames(lefse.genus) <- rownames(lefse.genus) %>% gsub("D_0","k",.) %>% gsub("\\.D_1","|p",.) %>% gsub("\\.D_2","|c",.) %>% gsub("\\.D_3","|o",.) %>% gsub("\\.D_4","|f",.) %>% gsub("\\.D_5","|g",.) %>% gsub("\\.D_6","|s",.) 

lefse.genus[1,lefse.genus[1,]==1] <- "Depression"
lefse.genus[1,lefse.genus[1,]==0] <- "Non-Depression"
write.table(lefse.genus,"Ref65.lefse.genus.dada.txt",sep="\t",quote = F,col.names = F)


###compared with systematic review
library(metamicrobiomeR) 

taz.nor.com <- cbind(metadata[,c("AgeAtSampling","BMI","patient_number","AnatomicLocation","Gender","depression")],
                     All.tax[pmatch(common.id,rownames(All.tax)),])
genus.nor.com <- cbind(metadata[,c("AgeAtSampling","BMI","patient_number","AnatomicLocation","Gender","depression")],
                       genus.nor[pmatch(common.id,rownames(genus.nor)),])
taz.nor.com$depression <- gsub(1,"1_depression",taz.nor.com$depression )
taz.nor.com$depression <- gsub(0,"0_depression",taz.nor.com$depression )

genus.nor.com$depression <- gsub(1,"1_depression",genus.nor.com$depression )
genus.nor.com$depression <- gsub(0,"0_depression",genus.nor.com$depression )


colnames(taz.nor.com) <- colnames(taz.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.)
colnames(genus.nor.com) <- colnames(genus.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.) 

colnames(taz.nor.com)[1:6] <- c("Age","BMI","PatientID","AnatomicLocation","Gender","Depression") 
colnames(genus.nor.com)[1:6] <- c("Age","BMI","PatientID","AnatomicLocation","Gender","Depression") 

taxacom.tax<-taxa.compare(taxtab=taz.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI","AnatomicLocation","Gender"),longitudinal="yes",p.adjust.method="fdr", personid = "PatientID")

taxacom.genus<-taxa.compare(taxtab=genus.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Age","BMI","AnatomicLocation","Gender"),longitudinal="yes",p.adjust.method="fdr", personid = "PatientID")

write.table(taxacom.genus,"Ref65.metaR.genus.dada2.result.txt",sep="\t",quote = F,row.names = F,col.names = T)
write.table(taxacom.tax,"Ref65.metaR.tax.result.dada2.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.genus, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)

taxcomtab.show(taxcomtab=taxacom.tax, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
####output MMUPHin format
rownames(All.tax) <- paste("Ref65.",rownames(All.tax),sep="")
write.table(t(All.tax),"Ref65.abd.dada2.MMUPHin.txt",sep="\t",quote = F)
metadata.MMUPHin <- metadata[,c(1,2,3,7,8,9,10,11)]
colnames(metadata.MMUPHin) <- c("SampleID","Anxiety_score","depression_score","AnatomicLocation","Sex","Age","BMI","Depression")
rownames(metadata.MMUPHin) <- paste("Ref65.",rownames(metadata.MMUPHin),sep="")
metadata.MMUPHin$subjectID <- "Ref65"
write.table(metadata.MMUPHin,"Ref65.meta.dada2.MMUPHin.txt",sep="\t",quote = F)

