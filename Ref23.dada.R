###Review 23 dada
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
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref23/")

###marker in systematic review
decrease <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/decrease.maker.list.txt",sep="\t")
increase <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/increase.maker.list.txt",sep="\t")
##metadata
metadata <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/DataDownload/ref23/match_sample_metadata.mvtitle.txt",sep="\t",row.names = 1,head=T,
                       colClasses=c("character","character","numeric","numeric","character","character","character","character","character","numeric","numeric","numeric","numeric"))


######################
######################
##1.alpha diversity

#observe otu
obt <- read.table("dada2/observed_otus_vector.qza/b38e6dc5-57e5-47f8-a47f-4e1645eeda8f/data/alpha-diversity.tsv",check.names = T)

#chao1
chao1.d <- read.table("dada2/chao1_alpha_diverisyt.qza/0b6c6832-c584-45cc-b6e9-0652d5924de9/data/alpha-diversity.tsv",check.names = T)

#shannon
shannon <- read.table("dada2/shannon_vector.qza/8e001ae8-393e-49b3-af3f-c64db8ba3a32/data/alpha-diversity.tsv",check.names = T)

#simpson
simpson <- read.table("dada2/simpson_alpha_diverisyt.qza/6805d56e-ab1d-4f95-aec4-9c8215a99f7d/data/alpha-diversity.tsv",check.names = T)

##2.taxonomy
kingdom <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-1.csv",colClasses=c("index"="character"))
phylum <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-2.csv",colClasses=c("index"="character"))
class <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-3.csv",colClasses=c("index"="character"))
ord <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-4.csv",colClasses=c("index"="character"))
family <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-5.csv",colClasses=c("index"="character"))
genus <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-6.csv",colClasses=c("index"="character"))
species <- read.csv("dada2/taxarare-bar-plots.qzv/262d70ed-7b64-4eec-ac9b-3f3be9f777b8/data/level-7.csv",colClasses=c("index"="character"))
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
kingdom <- kingdom[,-c(1,4:17),drop=F]
kingdom.nor <- data.Normalization(kingdom,type = "n10",normalization = "row")
colnames(phylum)
phylum <- phylum[,-c(1,27:39),drop=F]
phylum.nor <- data.Normalization(phylum,type = "n10", normalization = "row")
colnames(class)
class <- class[,-c(1,41:53),drop=F]
class.nor <- data.Normalization(class,type = "n10", normalization = "row")
colnames(ord)
ord <- ord[,-c(1,91:103),drop=F]
ord.nor <- data.Normalization(ord,type = "n10", normalization = "row")
colnames(family)
family <- family[,-c(1,202:214),drop=F]
family.nor <- data.Normalization(family,type = "n10", normalization = "row")
colnames(genus)
genus <- genus[,-c(1,666:678),drop=F]
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

###all alpha diversity significant different with normal distribution 
#####alpha diversity test
wil.res <- c()
t.res <- c()
for(i in 1:4){
    p.t <- t.test(alpha[,i] ~ alpha[,"mental_illness_type_depression_num"])$p.value
    p.v <- wilcox.test(alpha[,i] ~ alpha[,"mental_illness_type_depression_num"])$p.value
    means <- as.matrix(by(alpha[,i],alpha[,"mental_illness_type_depression_num"],mean))
    stand.d <- as.matrix(by(alpha[,i],alpha[,"mental_illness_type_depression_num"],sd))
    t.res<-rbind(t.res,c(p.t,means,stand.d))
    iqr <- as.matrix(by(alpha[,i],alpha[,"mental_illness_type_depression_num"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(alpha)[1:4]

colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")

names(wil.res[wil.res<0.05])
names(t.res[t.res<0.05])

library("kableExtra")


wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

t.res <- cbind(t.res,dim(alpha)[1],table(alpha$mental_illness_type_depression)[1],table(alpha$mental_illness_type_depression)[2])
colnames(t.res)[6:8] <- c("Total sample size","non-depression","depression")

wil.res <- cbind(wil.res,dim(alpha)[1],table(alpha$mental_illness_type_depression)[1],table(alpha$mental_illness_type_depression)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")


write.table(t.res,"Ref23_dada_alpha.diversity.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref23_dada_alpha.diversity.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)


###########alpha divesity odd ratio
alpha.odd <- alpha[,c(1:4,6,7,10,13,16)]
colnames(alpha.odd)[5:9]<-c("Age","BMI","Sex","Depression","subclass")
odd.crude.res <- c()
odd.adjust.res <- c()
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
                           round(exp(cbind(coef(reg), confint(reg))),3)[2,])
    
    
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c("Depression","Age","BMI","Sex","subclass")])
    colnames(data)[1] <- "alpha"
    reg.adjust <- clogit(Depression ~alpha+Age+BMI+Sex+strata(subclass) , data=data)
    #reg.adjust.2 <-clogit(Depression ~alpha+strata(subclass) , data=data)
    
    gee.adjust <- geeglm(Depression ~alpha+Age+BMI+Sex, family=binomial, data=data, 
                         id=subclass, corstr = "independence")
    
    odd.adjust.res <- rbind(odd.adjust.res,c(exp(coefficients(reg.adjust)[1]),exp(confint(reg.adjust)[1,])))
    odd.gee.res <- rbind(odd.gee.res,
                         exp(c(esticon(gee.adjust, c(0,1,0,0,0))$estimate,esticon(gee.adjust,c(0,1,0,0,0))$lwr,esticon(gee.adjust, c(0,1,0,0,0))$upr)))
    
    
    ###seperated by median
    ##crude
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg <- glm(Depression ~H.alpha, data=data, family=binomial)
    odd.crude.res <- rbind(odd.crude.res,
                           round(exp(cbind(coef(reg), confint(reg))),3)[2,])
    
    
    data <- cbind(alpha.odd[,i,drop=F],alpha.odd[,c("Depression","Age","BMI","Sex","subclass")])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg.adjust <- clogit(Depression ~H.alpha+Age+BMI+Sex+strata(subclass) , data=data)
    #reg.adjust.2 <-clogit(Depression ~alpha+strata(subclass) , data=data)
    
    gee.adjust <- geeglm(Depression ~H.alpha+Age+BMI+Sex, family=binomial, data=data, 
                         id=subclass, corstr = "independence")
    
    odd.adjust.res <- rbind(odd.adjust.res,c(exp(coefficients(reg.adjust)[1]),exp(confint(reg.adjust)[1,])))
    odd.gee.res <- rbind(odd.gee.res,
                         exp(c(esticon(gee.adjust, c(0,1,0,0,0))$estimate,esticon(gee.adjust,c(0,1,0,0,0))$lwr,esticon(gee.adjust, c(0,1,0,0,0))$upr)))
    
}
rownames(odd.crude.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.adjust.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.gee.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))

odd.crude.res
odd.crude.res <- cbind(odd.crude.res,dim(alpha)[1],table(alpha$mental_illness_type_depression)[1],table(alpha$mental_illness_type_depression)[2])
colnames(odd.crude.res)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.adjust.res
odd.adjust.res <- cbind(odd.adjust.res,dim(alpha)[1],table(alpha$mental_illness_type_depression)[1],table(alpha$mental_illness_type_depression)[2])
colnames(odd.adjust.res)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.gee.res
odd.gee.res <- cbind(odd.gee.res,dim(alpha)[1],table(alpha$mental_illness_type_depression)[1],table(alpha$mental_illness_type_depression)[2])
colnames(odd.gee.res)[1:3] <- c("","2.5 %","97.5 %")
colnames(odd.gee.res)[4:6] <- c("Total_sample_size","non_depression","depression")


write.table(odd.crude.res,"Ref23_dada_alpha.diversity.oR.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.adjust.res,"Ref23_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.res,"Ref23_dada_alpha.diversity.OR_GEE.result.txt",sep="\t",quote = F,col.names = T,row.names = T)


###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum ################################################################
###phylum distribution
phylum.d <- cbind.data.frame(phylum.nor,metadata)
sort(apply(phylum.d[,c(1:25)],2,function(x){sum(x!=0)}),decreasing = T)
phy.p <- colnames(phylum.d)[
    which(
        apply(phylum.d[,c(1:25)],2,function(x){sum(x!=0)})/dim(phylum.d)[1] > 0.2
    )
]

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
wil.res <- c()
t.res <- c()

for(i in 1:25){
    p.t <- t.test(phylum.d[,i] ~ phylum.d[,"mental_illness_type_depression_num"])$p.value
    p.v <- wilcox.test(phylum.d[,i] ~ phylum.d[,"mental_illness_type_depression_num"])$p.value
    
    means <- as.matrix(by(phylum.d[,i],phylum.d[,"mental_illness_type_depression_num"],mean))
    stand.d <- as.matrix(by(phylum.d[,i],phylum.d[,"mental_illness_type_depression_num"],sd))
    t.res<-rbind(t.res,c(p.t,means,stand.d))
    
    iqr <- as.matrix(by(phylum.d[,i],phylum.d[,"mental_illness_type_depression_num"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(phylum.d)[1:25]

colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(4) %>% filter(P.value <0.05) %>%
    kbl(caption = "phylum marker") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

t.res <- cbind(t.res,dim(phylum)[1])
colnames(t.res)[6] <- "Total sample size"

wil.res <- cbind(wil.res,dim(phylum)[1],table(phylum.d$mental_illness_type_depression)[1],table(phylum.d$mental_illness_type_depression)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")

write.table(t.res,"Ref23.dada_phylum.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref23.dada_phylum.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

b=1
id <- rownames(wil.res[wil.res[,1]<0.05,])
for(i in id){
    dat <- phylum.d[,c(i,"mental_illness_type_depression_num")]
    colnames(dat) <- c("phylum","Depression")
    p<-ggboxplot(dat,x="Depression", y="phylum",color="Depression",palette = "jco", add = "jitter", short.panel.labs = F,ylab=i,main=i,ggtheme=theme_classic())
    txt<-paste("box",b,"=p",sep="")
    eval((parse(text =txt)))
    b=b+1
}
grid.arrange(box1,box2,ncol=2,nrow=2)


#####phylum test by transfer data
###phylum transfer
sort(apply(phylum.d[,c(1:25)],2,function(x){sum(x!=0)}),decreasing = T)
id <- names(sort(apply(phylum.d[,c(1:25)],2,function(x){sum(x!=0)}),decreasing = T))

phylum.transf <- matrix(data=NA,nrow = 772,ncol = length(phy.p))
for(i in 1:length(phy.p)){
    n <- id[i]
    pt <- summary(powerTransform(phylum[,n] ~ 1, phylum, family = "bcnPower" ))
    data.transfer <- phylum[,n] ** pt$result[1]
    phylum.transf[,i] <- data.transfer
    
}
colnames(phylum.transf) <- id[1:length(phy.p)]
rownames(phylum.transf) <- rownames(phylum.d)
phylum.transf <- as.data.frame(phylum.transf)

b=1
nor.p <- c()
for(i in phy.p){
    dat <- phylum.transf[,i,drop=F]
    nor.p <- c(nor.p, shapiro.test(dat[,1])$p.value)
    txt <- paste("p <- ggplot(dat, aes(x=",i,"))",sep="" )
    eval((parse(text =txt)))
    
    p <- p +
        geom_histogram(fill="black", alpha=0.5, position="identity")
    
    txt<-paste("box",b,"=p",sep="")
    eval((parse(text =txt)))
    b=b+1
}
grid.arrange(box1,box2,box3,box4,box5,box6,box7,box8,ncol=2,nrow=4)
names(nor.p) <- phy.p
nor.p
###still significant different with normal distribution.

###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus ################################################################
###genus distribution
genus.d <- cbind.data.frame(genus.nor,metadata)
colnames(genus.d)
sum(apply(genus,2,function(x){sum(x!=0)})==48)
###genus transfer
# names(sort(apply(genus.d[,c(1:dim(genus)[2])],2,function(x){sum(x!=0)}),decreasing = T))[1:20]
# id <- which(apply(genus.d[,c(1:dim(genus)[2])],2,function(x){sum(x!=0)})==772)
# 
# df <- matrix(data=NA,nrow = 772,ncol = length(id))
# p<-c()
# for(i in 1:length(id)){
#     n <- id[i]
#     pt <- summary(powerTransform(genus[,n] ~ 1, genus))
#     data.transfer <- genus[,n] ** pt$result[1]
#     df[,i] <- data.transfer
#     p1 <- shapiro.test(genus[,n])$p.value
#     p2 <- shapiro.test(data.transfer)$p.value
#     p <- rbind(p,c(p1,p2))
#     
# }
# colnames(df) <- names(id)
# genus.transf <- cbind.data.frame(df,metadata)
# rownames(genus.transf) <- rownames(genus.d)
# 
# colnames(p) <- c("before transform","after transform")
# rownames(p) <- names(id)
# id.transf <- rownames(p[p[,2]>0.05,])
# 
# df <- as.data.frame(df)
# b=1
# for(i in 1:length(id)){
#     dat <- cbind.data.frame(genus.d[,id[i]],df[,i])
#     colnames(dat) <- c("Raw","Transfer")
#     p1.t <- ggplot(dat, aes(x=Raw)) +
#         geom_histogram(fill="black", alpha=0.5, position="identity" )
#     p2.t <- ggplot(dat, aes(x=Transfer)) +
#         geom_histogram(fill="black", alpha=0.5, position="identity" )
#     
#     txt<-paste("box.raw",b,"=p1.t",sep="")
#     eval((parse(text =txt)))
#     
#     txt<-paste("box.transfer",b,"=p2.t",sep="")
#     eval((parse(text =txt)))
#     b=b+1
# }
# b
# grid.arrange(box.raw1,box.raw2,box.raw3,box.raw4,box.raw5,box.raw6,box.raw7,ncol=2,nrow=4)
# grid.arrange(box.transfer1,box.transfer2,box.transfer3,box.transfer4,box.transfer5,box.transfer6,box.transfer7,ncol=2,nrow=4)

###the transform result is not significantly different with normal distribution, use t-test for these genus
#####genus test
wil.res <- c()
t.res <- c()

for(i in 1:664){
    ##test for normal distribution.
    ##From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
    p.normal <- shapiro.test(genus.d[,i])$p.value
    
    p.t <- t.test(genus.d[,i] ~ genus.d[,"mental_illness_type_depression_num"])$p.value
    p.v <- wilcox.test(genus.d[,i] ~ genus.d[,"mental_illness_type_depression_num"])$p.value
    means <- as.matrix(by(genus.d[,i],genus.d[,"mental_illness_type_depression_num"],mean))
    stand.d <- as.matrix(by(genus.d[,i],genus.d[,"mental_illness_type_depression_num"],sd))
    t.res<-rbind(t.res,c(p.normal,p.t,means,stand.d))
    iqr <- as.matrix(by(genus.d[,i],genus.d[,"mental_illness_type_depression_num"],quantile))
    wil.res<-wil.res<-rbind(wil.res,c(p.v,iqr[[1]][c(2,3,4)],iqr[[2]][c(2,3,4)]))
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(genus.d)[1:664]

colnames(t.res) <-  c("P.value_normal","P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")


colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")


wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(4) %>% filter(P.value <0.05) %>%
    kbl(caption = "phylum marker") %>%
    kable_classic(full_width = F, html_font = "Cambria") 


wil.res <- cbind(wil.res,dim(genus.d)[1],table(genus.d$mental_illness_type_depression)[1],table(genus.d$mental_illness_type_depression)[2])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")

mark.1 <- wil.res[wil.res[,2]<0.05 & wil.res[,1]<0.05,]
mark.2 <- t.res[t.res[,2]<0.05 & t.res[,1]>=0.05,]
mark.3 <- t.res.transf[t.res.transf[,2] <0.05,]
##no significant genus in transform profile

id1 <- rownames(mark.1)
id2 <- rownames(mark.2)
length(id1)
length(id2)


intersect(names(id1),names(id2))

write.table(t.res,"Ref23_dada.genus.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"Ref23_dada.genus.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

genus.m <- strsplit(id1,split='D_5__')
genus.m.g <-  sapply(genus.m, function(x){x[2]})
genus.m.g <- cbind(genus.m.g,"decrease")
genus.m.g[mark.1[,2]<mark.1[,3],2] <- "increase"

maker <- genus.m.g
##consist
pmatch(increase$Var1,maker[maker[,2]=="increase",1])
pmatch(decrease$Var1,maker[maker[,2]=="decrease",1])
##conflict
pmatch(increase$Var1,maker[maker[,2]=="decrease",1])
pmatch(decrease$Var1,maker[maker[,2]=="increase",1])
####output lefseq
genus.d <- cbind.data.frame(genus,metadata)
colnames(genus.d)[1:10]
lefse.genus <- t(genus.d[,c(673,1:664)])
rownames(lefse.genus) <- rownames(lefse.genus) %>% gsub("D_0","k",.) %>% gsub("\\.D_1","|p",.) %>% gsub("\\.D_2","|c",.) %>% gsub("\\.D_3","|o",.) %>% gsub("\\.D_4","|f",.) %>% gsub("\\.D_5","|g",.) %>% gsub("\\.D_6","|s",.) 

lefse.genus[1,lefse.genus[1,]==1] <- "Depression"
lefse.genus[1,lefse.genus[1,]==0] <- "Non-Depression"
write.table(lefse.genus,"lefse.genus.dada.txt",sep="\t",quote = F,col.names = F)

###compared with systematic review
library(metamicrobiomeR) 

taz.nor.com <- cbind(metadata[,c("sample","sex","age_years","bmi_corrected","mental_illness_type_depression_num")],All.tax)
genus.nor.com <- cbind(metadata[,c("sample","sex","age_years","bmi_corrected","mental_illness_type_depression_num")],genus.nor)
taz.nor.com$mental_illness_type_depression_num <- gsub(1,"1_depression",taz.nor.com$mental_illness_type_depression_num )
taz.nor.com$mental_illness_type_depression_num <- gsub(0,"0_depression",taz.nor.com$mental_illness_type_depression_num )

genus.nor.com$mental_illness_type_depression_num <- gsub(1,"1_depression",genus.nor.com$mental_illness_type_depression_num )
genus.nor.com$mental_illness_type_depression_num <- gsub(0,"0_depression",genus.nor.com$mental_illness_type_depression_num )


colnames(taz.nor.com) <- colnames(taz.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.)
colnames(taz.nor.com)[1:5] <- c("PatientID","Gender","Age","BMI","Depression") 
colnames(genus.nor.com) <- colnames(genus.nor.com) %>% gsub("D_0","k",.) %>% gsub("D_1","p",.) %>% gsub("D_2","c",.) %>% gsub("D_3","o",.) %>% gsub("D_4","f",.) %>% gsub("D_5","g",.) %>% gsub("D_6","s",.) 
colnames(genus.nor.com)[1:5] <- c("PatientID","Gender","Age","BMI","Depression") 
taxacom.tax<-taxa.compare(taxtab=taz.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Gender","Age","BMI"),longitudinal="no",p.adjust.method="fdr")

taxacom.genus<-taxa.compare(taxtab=genus.nor.com,propmed.rel="gamlss",comvar="Depression",adjustvar=c("Gender","Age","BMI"),longitudinal="no",p.adjust.method="fdr")

write.table(taxacom.genus,"Ref23.metaR.genus.dada.result.txt",sep="\t",quote = F,row.names = F,col.names = T)
write.table(taxacom.tax,"Ref23.metaR.tax.result.dada.txt",sep="\t",quote = F,row.names = F,col.names = T)

taxcomtab.show(taxcomtab=taxacom.genus, tax.lev="l5",tax.select="none",showvar="Depression1_depression",p.cutoff=1)

taxcomtab.show(taxcomtab=taxacom.tax, tax.lev="l3",tax.select="none",showvar="Depression1_depression",p.cutoff=0.05)
####output MMUPHin format
rownames(All.tax) <- paste("Ref23.",rownames(All.tax),sep="")
write.table(t(All.tax),"Ref23.abd.dada.MMUPHin.txt",sep="\t",quote = F)

metadata.MMUPHin <- metadata[,c(2,3,6,9,11)]
rownames(metadata.MMUPHin) <- paste("Ref23.",rownames(metadata.MMUPHin),sep="")
colnames(metadata.MMUPHin) <- c("Age","BMI","Sex","Depression","ps_group")
metadata.MMUPHin$subjectID <- "Ref23"
write.table(metadata.MMUPHin,"Ref23.meta.dada.MMUPHin.txt",sep="\t",quote = F)

