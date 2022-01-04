###Review 79 dada2
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
setwd("/Users/heinmintun/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref79/")

###marker in systematic review
decrease <- read.table("/Users/heinmintun/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/decrease.maker.list.txt",sep="\t")
increase <- read.table("/Users/heinmintun/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/increase.maker.list.txt",sep="\t")
##metadata
metadata <- read.table("/Users/heinmintun/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/DataDownload/Ref79/ANMAP.txt",sep="\t",row.names = 1,head=T,comment.char = "")
##depression cutoff
##in the paper: .
head(metadata)
metadata<-metadata[-1,]
metadata$thispapercatolog <- 0
metadata$thispapercatolog[metadata$Group  =="PRE"] <- 1
table(metadata$thispapercatolog)
metadata$sample <- substring(rownames(metadata),1,3)

######################
######################
##1.alpha diversity

#observe otu
obt <- read.table("dada2/observed_otus_vector.qza/12eee2b9-f796-4f25-86fa-92be1250b140/data/alpha-diversity.tsv")

#chao1
chao1.d <- read.table("dada2/chao1_alpha_diverisyt.qza/682853ee-88df-4dbe-87a3-2a9deb603a97/data/alpha-diversity.tsv")

#shannon
shannon <- read.table("dada2/shannon_vector.qza/43b6d600-0d17-4f84-892d-b291325df45c/data/alpha-diversity.tsv")

#simpson
simpson <- read.table("dada2/simpson_alpha_diverisyt.qza/74a4f0e0-7017-43e6-b278-2d32618181fb/data/alpha-diversity.tsv")

##2.taxonomy
kingdom <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-1.csv")
phylum <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-2.csv")
class <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-3.csv")
ord <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-4.csv")
family <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-5.csv")
genus <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-6.csv")
species <- read.csv("dada2/taxarare-bar-plots.qzv/555befc2-31f0-4dbc-8a1c-f026b3e9164b/data/level-7.csv")
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
kingdom <- kingdom[,-c(1,3:4),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
colnames(phylum)
phylum <- phylum[,-c(1,12:13),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
colnames(class)
class <- class[,-c(1,18:19),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
colnames(ord)
ord <- ord[,-c(1,20:21),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
colnames(family)
family <- family[,-c(1,44:45),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
colnames(genus)
genus <- genus[,-c(1,185:186),drop=F]
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

pid <- pmatch(rownames(chao1.d),rownames(metadata))
metadata <- metadata[pid,]

alpha <- cbind.data.frame(obt,chao1.d,simpson,shannon,metadata)
colnames(alpha)
write.table(alpha,"Ref79_aloha_diversity.txt",sep="\t",quote = F,col.names = T,row.names = T)
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
wil.res <- c()
t.res <- c()
colnames(alpha)
for(i in c(1:4)){
    p.t <- t.test(alpha[,i] ~ alpha[,"thispapercatolog"])$p.value
    p.v <- wilcox.test(alpha[,i] ~ alpha[,"thispapercatolog"])$p.value
    means <- as.matrix(by(alpha[,i],alpha[,"thispapercatolog"],mean))
    stand.d <- as.matrix(by(alpha[,i],alpha[,"thispapercatolog"],sd))
    t.res<-rbind(t.res,c(p.t,means,stand.d))
    iqr <- as.matrix(by(alpha[,i],alpha[,"thispapercatolog"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(alpha)[c(1:4)]

colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


wil.res[wil.res[,1]<0.05,]
t.res[t.res[,1]<0.05,]

library("kableExtra")

wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

t.res <- cbind(t.res,dim(alpha)[1],table(metadata$thispapercatolog)[1],table(metadata$thispapercatolog)[2])
colnames(t.res)[6:8] <-  c("Total sample size","non-depression","depression")

wil.res <- cbind(wil.res,dim(alpha)[1],table(metadata$thispapercatolog)[1],table(metadata$thispapercatolog)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")

write.table(t.res,"Ref79_dada2_alpha.diversity.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref79_dada2_alpha.diversity.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

###########alpha divesity odd ratio
alpha.odd <- alpha[,c(1:4,8,11,10)]
colnames(alpha.odd)[c(5:7)]<-c("Group","participant","Depression")
write.table(alpha.odd,"alpha_diverisy.rawdata.txt",quote = F,sep = "\t")
odd.crude.res <- c()
odd.gee.res <- c()

library(survival)
library(geepack)
library(doBy)
for(i in 1:4){
    ##crude
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    reg <- glm(Depression ~alpha, data=data, family=binomial)
    odd.crude.res <- rbind(odd.crude.res,
                           round(cbind(coef(reg), confint(reg)),3)[2,])
    
    ##adjust for age,sex,ibs,stress
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c(6:7)])
    colnames(data)[1] <- "alpha"
    gee.adjust <- geeglm(Depression ~alpha, family=binomial, data=data, 
                         id=participant, corstr = "independence")
    odd.gee.res <- rbind(odd.gee.res,
                         c(esticon(gee.adjust, c(0,1))$estimate,esticon(gee.adjust,c(0,1))$lwr,esticon(gee.adjust, c(0,1))$upr))

    ###seperated by median
    ##crude
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg <- glm(Depression ~H.alpha, data=data, family=binomial)
    odd.crude.res <- rbind(odd.crude.res,
                           #round(exp(cbind(coef(reg), confint(reg))),3)[2,])
                           exp(cbind(coef(reg), confint(reg)))[2,])
    
    ##adjust for age,sex,ibs,stress
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c(6:7)])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    gee.adjust <- geeglm(Depression ~H.alpha, family=binomial, data=data, 
                         id=participant, corstr = "independence")
    odd.gee.res <- rbind(odd.gee.res,
                         exp(c(esticon(gee.adjust, c(0,1))$estimate,esticon(gee.adjust,c(0,1))$lwr,esticon(gee.adjust, c(0,1))$upr)))
    
}
rownames(odd.crude.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.gee.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))


odd.crude.res
odd.crude.res <- cbind(odd.crude.res,dim(alpha)[1],table(alpha$thispapercatolog)[1],table(alpha$thispapercatolog)[2])
colnames(odd.crude.res)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.gee.res
odd.gee.res <- cbind(odd.gee.res,dim(alpha)[1],table(alpha$thispapercatolog)[1],table(alpha$thispapercatolog)[2])
colnames(odd.gee.res)[4:6] <- c("Total_sample_size","non_depression","depression")



write.table(odd.crude.res,"Ref79.dada_alpha.diversity.oR.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.res,"Ref79.dada_alpha.diversity.OR_adjust.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum distribution
phylum.d <- cbind.data.frame(phylum.nor,metadata)
sort(apply(phylum.d[,c(1:10)],2,function(x){sum(x!=0)}),decreasing = T)
phy.p <- colnames(phylum.d)[
    which(
        apply(phylum.d[,c(1:10)],2,function(x){sum(x!=0)})/dim(phylum.d)[1] > 0.2
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


#####phylum test by wiltest and t.test
wil.res <- c()
t.res <- c()

for(i in 1:10){
    p.t <- t.test(phylum.d[,i] ~ phylum.d[,"thispapercatolog"])$p.value
    p.v <- wilcox.test(phylum.d[,i] ~ phylum.d[,"thispapercatolog"])$p.value
    
    means <- as.matrix(by(phylum.d[,i],phylum.d[,"thispapercatolog"],mean))
    stand.d <- as.matrix(by(phylum.d[,i],phylum.d[,"thispapercatolog"],sd))
    t.res<-rbind(t.res,c(p.t,means,stand.d))
    
    iqr <- as.matrix(by(phylum.d[,i],phylum.d[,"thispapercatolog"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
    
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(phylum.d)[1:10]

colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


wil.res[wil.res[,1]<0.05,]
t.res[t.res[,1]<0.05,]

phylum.marker <- wil.res[wil.res[,1]<0.05,]

wil.res[,c(1,3,6)] %>% as.data.frame() %>%  filter(P.value <0.05) %>% round(2) %>%
    kbl(caption = "Phylum marker") %>%
    kable_classic(full_width = F, html_font = "Cambria")

wil.res <- cbind(wil.res,dim(phylum)[1],table(metadata$thispapercatolog)[1],table(metadata$thispapercatolog)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")


write.table(t.res,"Ref79.dada2_phylum.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref79.dada2_phylum.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus distribution
genus.d <- cbind.data.frame(genus.nor,metadata)
colnames(genus.d)
sum(apply(genus,2,function(x){sum(x!=0)})==25)
###genus transfer
names(sort(apply(genus.d[,c(1:dim(genus)[2])],2,function(x){sum(x!=0)}),decreasing = T))[1:20]
id <- which(apply(genus.d[,c(1:dim(genus)[2])],2,function(x){sum(x!=0)})==25)

df <- matrix(data=NA,nrow = 25,ncol = length(id))
p<-c()
for(i in 1:length(id)){
    n <- id[i]
    pt <- summary(powerTransform(genus[,n] ~ 1, genus))
    data.transfer <- genus[,n] ** pt$result[1]
    df[,i] <- data.transfer
    p1 <- shapiro.test(genus[,n])$p.value
    p2 <- shapiro.test(data.transfer)$p.value
    p <- rbind(p,c(p1,p2))
    
}
colnames(df) <- names(id)
genus.transf <- cbind.data.frame(df,metadata)
rownames(genus.transf) <- rownames(genus.d)

colnames(p) <- c("before transform","after transform")
rownames(p) <- names(id)
id.transf <- rownames(p[p[,2]>0.05,])

df <- as.data.frame(df)
b=1
for(i in 1:length(id)){
    dat <- cbind.data.frame(genus.d[,id[i]],df[,i])
    colnames(dat) <- c("Raw","Transfer")
    p1.t <- ggplot(dat, aes(x=Raw)) +
        geom_histogram(fill="black", alpha=0.5, position="identity" )
    p2.t <- ggplot(dat, aes(x=Transfer)) +
        geom_histogram(fill="black", alpha=0.5, position="identity" )
    
    txt<-paste("box.raw",b,"=p1.t",sep="")
    eval((parse(text =txt)))
    
    txt<-paste("box.transfer",b,"=p2.t",sep="")
    eval((parse(text =txt)))
    b=b+1
}
b
grid.arrange(box.raw1,box.raw2,ncol=2,nrow=4)
grid.arrange(box.transfer1,box.transfer2,ncol=2,nrow=4)

###the transform result is not significantly different with normal distribution, use t-test for these genus
#####genus test
wil.res <- c()
t.res <- c()

for(i in 1:183){
    ##test for normal distribution.
    ##From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
    p.normal <- shapiro.test(genus.d[,i])$p.value
    
    p.t <- t.test(genus.d[,i] ~ genus.d[,"thispapercatolog"])$p.value
    p.v <- wilcox.test(genus.d[,i] ~ genus.d[,"thispapercatolog"])$p.value
    means <- as.matrix(by(genus.d[,i],genus.d[,"thispapercatolog"],mean))
    stand.d <- as.matrix(by(genus.d[,i],genus.d[,"thispapercatolog"],sd))
    t.res<-rbind(t.res,c(p.normal,p.t,means,stand.d))
    iqr <- as.matrix(by(genus.d[,i],genus.d[,"thispapercatolog"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
    
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(genus.d)[1:183]

colnames(t.res) <-  c("P.value_normal","P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


mark.1 <- wil.res[wil.res[,2]<0.05 & wil.res[,1]<0.05,]
mark.2 <- t.res[t.res[,2]<0.05 & t.res[,1]>=0.05,]
mark.3 <- t.res.transf[t.res.transf[,2] <0.05,]
##no significant genus in transform profile

id1 <- rownames(mark.1)
id2 <- rownames(mark.2)
length(id1)
length(id2)


intersect(names(id1),names(id2))

wil.res[,c(1,3,6)] %>% as.data.frame() %>%  filter(P.value <0.05) %>% round(4) %>%
    kbl(caption = "Phylum marker") %>%
    kable_classic(full_width = F, html_font = "Cambria")

t.res %>% as.data.frame() %>%  filter(P.value <0.05) %>%
    kbl(caption = "Genus marker") %>%
    kable_classic(full_width = F, html_font = "Cambria")

wil.res <- cbind(wil.res,dim(genus)[1],table(metadata$thispapercatolog)[1],table(metadata$thispapercatolog)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")


write.table(t.res,"Ref79_dada2.genus.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref79_dada2.genus.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

genus.m <- strsplit(id1,split='D_5__')
genus.m.g <-  sapply(genus.m, function(x){x[2]})
genus.m.g <- cbind(genus.m.g,"decrease")
genus.m.g[mark.1[,2]<mark.1[,3],2] <- "increase"

maker <- rbind(genus.m.g)
##consist
pmatch(increase$Var1,maker[maker[,2]=="increase",1])
pmatch(decrease$Var1,maker[maker[,2]=="decrease",1])
##conflict
pmatch(increase$Var1,maker[maker[,2]=="decrease",1])
pmatch(decrease$Var1,maker[maker[,2]=="increase",1])
####output lefseq
genus.d <- cbind.data.frame(genus,metadata)
colnames(genus.d)[1:10]
lefse.genus <- t(genus.d[,c(189,1:183)])
rownames(lefse.genus) <- rownames(lefse.genus) %>% gsub("D_0","k",.) %>% gsub("\\.D_1","|p",.) %>% gsub("\\.D_2","|c",.) %>% gsub("\\.D_3","|o",.) %>% gsub("\\.D_4","|f",.) %>% gsub("\\.D_5","|g",.) %>% gsub("\\.D_6","|s",.) 

lefse.genus[1,lefse.genus[1,]==1] <- "Depression"
lefse.genus[1,lefse.genus[1,]==0] <- "Non-Depression"
write.table(lefse.genus,"Ref79.lefse.genus.dada.txt",sep="\t",quote = F,col.names = F)


###compared with systematic review
library(metamicrobiomeR) 

taz.nor.com <- cbind(metadata[,c("thispapercatolog","sample")],All.tax)
genus.nor.com <- cbind(metadata[,c("thispapercatolog","sample")],genus.nor)
taz.nor.com$thispapercatolog <- gsub(1,"1_depression",taz.nor.com$thispapercatolog )
taz.nor.com$thispapercatolog <- gsub(0,"0_depression",taz.nor.com$thispapercatolog )

genus.nor.com$thispapercatolog <- gsub(1,"1_depression",genus.nor.com$thispapercatolog )
genus.nor.com$thispapercatolog <- gsub(0,"0_depression",genus.nor.com$thispapercatolog )


colnames(taz.nor.com) <- colnames(taz.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.)
colnames(genus.nor.com) <- colnames(genus.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.) 

colnames(taz.nor.com)[1:2] <- c("Depression","PatientID")
colnames(genus.nor.com)[1:2] <- c("Depression","PatientID")
taxacom.tax<-taxa.compare(taxtab=taz.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="yes",p.adjust.method="fdr",personid = "PatientID")

taxacom.genus<-taxa.compare(taxtab=taz.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=NULL,longitudinal="yes",p.adjust.method="fdr",personid = "PatientID")
write.table(taxacom.genus,"Ref79.metaR.genus.dada2.result.txt",sep="\t",quote = F,row.names = F,col.names = T)
write.table(taxacom.tax,"Ref79.metaR.tax.result.dada2.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.genus, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)

taxcomtab.show(taxcomtab=taxacom.tax, tax.lev="l6",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
####output MMUPHin format
rownames(All.tax) <- paste("Ref79.",rownames(All.tax),sep="")
write.table(t(All.tax),"Ref79.abd.dada2.MMUPHin.txt",sep="\t",quote = F)
metadata.MMUPHin <- metadata[,c(6,7)]
colnames(metadata.MMUPHin) <- c("Depression","PatientID")
rownames(metadata.MMUPHin) <- paste("Ref79.",rownames(metadata.MMUPHin),sep="")
metadata.MMUPHin$subjectID <- "Ref79"
write.table(metadata.MMUPHin,"Ref79.meta.dada2.MMUPHin.txt",sep="\t",quote = F)

