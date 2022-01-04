rm(list=ls())

ref1 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/Ref1_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)
ref23 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref23/Ref23_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)
ref27<- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref27/Ref27_dada_alpha.diversity.BDH.OR.result.txt",sep="\t",check.names = F)
ref40 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref40/Ref40_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)
ref60 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref60/Ref60_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)
ref77 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref77/Ref77_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)
ref79 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref79/Ref79.dada_alpha.diversity.OR_adjust.result.txt",sep="\t",check.names = F)

setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/alpha_diversity")
id2title <- read.table("../ID2Title_meta_analysis_study.txt",sep="\t",head=T)
id2title$id <- paste("ref",id2title$id,sep= "")

head(ref1)
head(ref23)
head(ref27)
head(ref40)
head(ref60)
head(ref77)
head(ref79)
all.equal(rownames(ref1) ,rownames(ref23), rownames(ref40),rownames(ref60),rownames(ref77),rownames(ref79),rownames(ref27))


library("tidyverse")
library("meta")
library("metafor")
file.id <-ls(pattern = "ref")

###obs_OTU higher
data.h <- c()
for(i in file.id){
    txt <- paste("dat.tmp <- ",i,sep="")
    eval((parse(text =txt)))
    colnames(dat.tmp)[1:3] <- c("obs_OTU","2.5 %","97.5 %")
    data.h <- rbind(data.h,dat.tmp[2,])
}
rownames(data.h) <- file.id

data.h$logOR <- log(data.h$obs_OTU)
data.h$se.logOR <- (log(data.h$`97.5 %`)- log(data.h$`2.5 %`))/(2*1.96)
data.h$slab <- id2title[pmatch(file.id,id2title$id),"title"]

obs_OTU_H.re <- rma(yi=logOR, sei=se.logOR, slab=slab, method="REML", data=data.h,measure ="OR" )

forest(obs_OTU_H.re, transf=exp, refline=1)

obs_OTU_H.re.metafor <- metagen(TE=logOR, seTE=se.logOR, studlab=slab, method.tau="REML",method.tau.ci="QP", sm="OR",data=data.h)

dev.off()
pdf("obs.adjust.model.pdf",width = 10,height = 5)
metafor::forest(obs_OTU_H.re.metafor,transf=exp, refline=1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Observed OTU")
dev.off()
###chao1 higher
n=4
data.h <- c()
for(i in file.id){
    txt <- paste("dat.tmp <- ",i,sep="")
    eval((parse(text =txt)))
    colnames(dat.tmp)[1:3] <- c("chao1","2.5 %","97.5 %")
    data.h <- rbind(data.h,dat.tmp[n,])
}
rownames(data.h) <- file.id

data.h$logOR <- log(data.h$chao1)
data.h$se.logOR <- (log(data.h$`97.5 %`)- log(data.h$`2.5 %`))/(2*1.96)
data.h$slab <- id2title[pmatch(file.id,id2title$id),"title"]
chao1_H.re <- rma(yi=logOR, sei=se.logOR, slab=slab, method="REML", data=data.h,measure ="OR" )

chao1.re.metafor <- metagen(TE=logOR, seTE=se.logOR, studlab=slab, method.tau="REML",method.tau.ci="QP", sm="OR",data=data.h)

pdf("chao1.adjust.model.pdf",width = 10,height = 5)
metafor::forest(chao1.re.metafor,transf=exp, refline=1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Chao1 index")
dev.off()
###simpson higher
n=6
data.h <- c()
for(i in file.id){
    txt <- paste("dat.tmp <- ",i,sep="")
    eval((parse(text =txt)))
    colnames(dat.tmp)[1:3] <- c("simpson","2.5 %","97.5 %")
    data.h <- rbind(data.h,dat.tmp[n,])
}
rownames(data.h) <- file.id

data.h$logOR <- log(data.h$simpson)
data.h$se.logOR <- (log(data.h$`97.5 %`)- log(data.h$`2.5 %`))/(2*1.96)
data.h$slab <- id2title[pmatch(file.id,id2title$id),"title"]
simpson_H.re <- rma(yi=logOR, sei=se.logOR, slab=slab, method="REML", data=data.h,measure ="OR" )
forest(simpson_H.re,transf = exp,refline = 1)

simpson.re.metafor <- metagen(TE=logOR, seTE=se.logOR, studlab=slab, method.tau="REML",method.tau.ci="QP", sm="OR",data=data.h)

pdf("simpson.adjust.model.pdf",width = 10,height = 5)
metafor::forest(simpson.re.metafor,transf=exp, refline=1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Simpson index")
dev.off()
###shannon higher
n=8
data.h <- c()
for(i in file.id){
    txt <- paste("dat.tmp <- ",i,sep="")
    eval((parse(text =txt)))
    colnames(dat.tmp)[1:3] <- c("shannon","2.5 %","97.5 %")
    data.h <- rbind(data.h,dat.tmp[n,])
}
rownames(data.h) <- file.id

data.h$logOR <- log(data.h$shannon)
data.h$se.logOR <- (log(data.h$`97.5 %`)- log(data.h$`2.5 %`))/(2*1.96)
data.h$slab <- id2title[pmatch(file.id,id2title$id),"title"]
shannon_H.re <- rma(yi=logOR, sei=se.logOR, slab=slab, method="REML", data=data.h,measure ="OR" )
shannon.re.metafor <- metagen(TE=logOR, seTE=se.logOR, studlab=slab, method.tau="REML",method.tau.ci="QP", sm="OR",data=data.h)

dev.off()
pdf("shannon.adjust.model.pdf",width = 10,height = 5)
metafor::forest(shannon.re.metafor,transf=exp, refline=1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Shannon index")
dev.off()


#No strong evidence of asymmetry (but note: small sample size)

leave1out(shannon_H.re, transf=exp)
leave1out(obs_OTU_H.re, transf=exp)
leave1out(chao1_H.re, transf=exp)
leave1out(simpson_H.re, transf=exp)

forest(obs_OTU_H.re, transf=exp, refline=1,header = "Observed OTU")
forest(chao1_H.re, transf=exp, refline=1,header = "Chao1 index")
forest(simpson_H.re, transf=exp, refline=1,header = "Simpson index")
forest(shannon_H.re, transf=exp, refline=1,header = "Shannon index")
funnel(shannon_H.re, atransf=exp)
funnel(simpson_H.re, atransf=exp)
funnel(chao1_H.re, atransf=exp)
funnel(obs_OTU_H.re, atransf=exp)
regtest(obs_OTU_H.re)
regtest(chao1_H.re)
regtest(simpson_H.re)
regtest(shannon_H.re)

