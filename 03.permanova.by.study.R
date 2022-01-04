###########silva138
###########silva138
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/")
rm(list=ls())

library(ade4)
library(vegan)
library(ggplot2)
library(viridis)
library(MMUPHin)
library(dplyr)
library(qiime2R)
library(boot)
library(erer)
library(ggpubr)
###ID and Title
id_title <- read.table("ID2Title_meta_analysis_study.txt",sep="\t",head=T)

##enterotype and permonova by study

source("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/tmp/script_final/enterotype.function.R")
####ref1
ref1 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/silva138/level-6.csv",header = T,check.names = F)
rownames(ref1) <- ref1$index
colnames(ref1)
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data

ref1.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/Ref1.meta.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref1.meta),paste("ref1.",rownames(ref1),sep="")) != c(1:dim(ref1)[1]))
rownames(ref1) <- paste("ref1.",rownames(ref1),sep="")
ref1 <- ref1[,-c(1,265:281)]
rowSums(ref1)

sum(pmatch(rownames(ref1.meta),paste("ref1.",labels(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref1.meta.rn <- ref1.meta
rownames(ref1.meta.rn) <- gsub("ref1.","",rownames(ref1.meta))

set.seed(1)
colnames(ref1.meta)

ref1.denose.adonis <- c()
ref1.boos.ci.adonis <- c()
ref1.f.ci.adonis <- c()
ref1.betadisper <- c()
ref1.betadisper.HSD <- c()
ref1.anosim <-c()
#depression
i <- 12
file <- ls(pattern=".dis$")

for(dist.file in file){
    rm.na.id <- !is.na(ref1.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref1.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref1.betadisper.HSD <- rbind(ref1.betadisper.HSD,ps.disper.HSD)
    
    ref1.betadisper <- rbind(ref1.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Depression <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Depression", y = "Distance_to_centroid",
                   fill = "Depression", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")+ylab(dist.file)
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref1.denose.adonis <- rbind(ref1.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref1.boos.ci.adonis <- rbind(ref1.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref1.f.ci.adonis <- rbind(ref1.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref1.anosim <- rbind(ref1.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}


rownames(ref1.denose.adonis) <- rownames(ref1.boos.ci.adonis) <- rownames(ref1.f.ci.adonis) <- rownames(ref1.betadisper) <- rownames(ref1.anosim) <- file

rownames(ref1.betadisper.HSD) <- paste(rownames(ref1.betadisper.HSD),file,sep="_")

list.adonis.ref1 <- list(adonis.output = ref1.denose.adonis, adonis.R.CI.bootrap = ref1.boos.ci.adonis, adonis.R.CI.F = ref1.f.ci.adonis, betadisper.output=ref1.betadisper,betadisper.HSD = ref1.betadisper.HSD,anosim = ref1.anosim)
save(z = list.adonis.ref1, file = "pernomova/ref1.output.RData",row.names = T)
ref1.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)

####ref23
####ref23
####ref23
####ref23
####ref23
####ref23

ref23 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/silva138/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref23) <- ref23$index
colnames(ref23)

ref23.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/Ref23.meta.dada.MMUPHin.txt",sep="\t",header = T,check.names = F)

####ref23
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref23.meta),paste("Ref23.",rownames(ref23),sep="")) != c(1:dim(ref23)[1]))
rownames(ref23) <- paste("Ref23.",rownames(ref23),sep="")
ref23 <- ref23[,-c(1,617:629)]
rowSums(ref23)

sum(pmatch(rownames(ref23.meta),paste("Ref23.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref23.meta.rn <- ref23.meta
rownames(ref23.meta.rn) <- gsub("Ref23.","",rownames(ref23.meta))

set.seed(1)
colnames(ref23.meta)

ref23.denose.adonis <- c()
ref23.boos.ci.adonis <- c()
ref23.f.ci.adonis <- c()
ref23.betadisper <- c()
ref23.betadisper.HSD <- c()
ref23.anosim <- c()
#depression
i <- 4

for(dist.file in file){
    rm.na.id <- !is.na(ref23.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref23.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref23.betadisper.HSD <- rbind(ref23.betadisper.HSD,ps.disper.HSD)
    ref23.betadisper <- rbind(ref23.betadisper,beta.disp$tab[1,])
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref23.denose.adonis <- rbind(ref23.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref23.boos.ci.adonis <- rbind(ref23.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref23.f.ci.adonis <- rbind(ref23.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref23.anosim <- rbind(ref23.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref23.denose.adonis) <- rownames(ref23.boos.ci.adonis) <- rownames(ref23.f.ci.adonis) <- rownames(ref23.betadisper) <- rownames(ref23.anosim) <-file

rownames(ref23.betadisper.HSD) <- paste(rownames(ref23.betadisper.HSD),file,sep="_")

list.adonis.ref23 <- list(adonis.output = ref23.denose.adonis, adonis.R.CI.bootrap = ref23.boos.ci.adonis, adonis.R.CI.F = ref23.f.ci.adonis, betadisper.output=ref23.betadisper,betadisper.HSD = ref23.betadisper.HSD,anosim = ref1.anosim)
save(z = list.adonis.ref23, file = "pernomova/ref23.output.RData",row.names = T)
ref23.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)

####ref27
####ref27
####ref27
####ref27
####ref27
####ref27


ref27 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref27) <- ref27$index
colnames(ref27)
ref27.meta <- ref27[,c(1,269:270)]
ref27 <- ref27[,-c(1,269:270)]
rowSums(ref27)

####ref27
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data

sum(pmatch(rownames(ref27.meta),label(unw.dis)) != c(1:length(labels(unw.dis))))

###BD vs H
set.seed(1)

ref27.denose.adonis <- c()
ref27.boos.ci.adonis <- c()
ref27.f.ci.adonis <- c()
ref27.betadisper <- c()
ref27.betadisper.HSD <- c()
ref27.anosim <-c()
#depression
i <- 4
for(dist.file in file){
    BDH.sample <- rownames(ref27.meta[ref27.meta$Time == 1,])
    ref27.meta.BDH <- ref27.meta[BDH.sample,]
    colnames(ref27.meta.BDH)
    ref27.meta.BDH$Depression <- 0
    ref27.meta.BDH$Depression[ref27.meta.BDH$Group == "BD"] <- 1
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[BDH.sample,BDH.sample]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref27.meta.BDH[,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref27.betadisper.HSD <- rbind(ref27.betadisper.HSD,ps.disper.HSD)
    ref27.betadisper <- rbind(ref27.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref27.denose.adonis <- rbind(ref27.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref27.boos.ci.adonis <- rbind(ref27.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref27.f.ci.adonis <- rbind(ref27.f.ci.adonis,out)
    
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref27.anosim <- rbind(ref27.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref27.denose.adonis) <- rownames(ref27.boos.ci.adonis) <- rownames(ref27.f.ci.adonis) <- rownames(ref27.betadisper) <- rownames(ref27.betadisper.HSD) <-  rownames(ref27.anosim)<-file

rownames(ref27.betadisper.HSD) <- paste(rownames(ref27.betadisper.HSD),file,sep="_")

list.adonis.ref27 <- list(adonis.output = ref27.denose.adonis, adonis.R.CI.bootrap = ref27.boos.ci.adonis, adonis.R.CI.F = ref27.f.ci.adonis, betadisper.output=ref27.betadisper,betadisper.HSD = ref27.betadisper.HSD,anosim = ref27.anosim)
save(z = list.adonis.ref27, file = "pernomova/ref27.output.RData",row.names = T)
ref27.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)

####ref40
####ref40
####ref40
####ref40
####ref40
ref40 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/silva138/level-6.csv",header = T,check.names = F)
rownames(ref40) <- ref40$index

ref40.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/Ref40.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref40.meta),paste("Ref40.",rownames(ref40),sep="")) != c(1:dim(ref40)[1]))
rownames(ref40) <-paste("Ref40.",rownames(ref40),sep="")
colnames(ref40)
ref40 <- ref40[,-c(1,245:299)]
rowSums(ref40)
####ref40
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref40.meta),paste("Ref40.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref40.meta.rn <- ref40.meta
rownames(ref40.meta.rn) <- gsub("Ref40.","",rownames(ref40.meta))

set.seed(1)
colnames(ref40.meta)

ref40.denose.adonis <- c()
ref40.boos.ci.adonis <- c()
ref40.f.ci.adonis <- c()
ref40.betadisper <- c()
ref40.betadisper.HSD <- c()
ref40.anosim <- c()
#depression
i <- 4

for(dist.file in file){
    rm.na.id <- !is.na(ref40.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref40.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref40.betadisper.HSD <- rbind(ref40.betadisper.HSD,ps.disper.HSD)
    ref40.betadisper <- rbind(ref40.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref40.denose.adonis <- rbind(ref40.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref40.boos.ci.adonis <- rbind(ref40.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref40.f.ci.adonis <- rbind(ref40.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref40.anosim <- rbind(ref40.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref40.denose.adonis) <- rownames(ref40.boos.ci.adonis) <- rownames(ref40.f.ci.adonis) <- rownames(ref40.betadisper) <-rownames(ref40.anosim) <- file

rownames(ref40.betadisper.HSD) <- paste(rownames(ref40.betadisper.HSD),file,sep="_")

list.adonis.ref40 <- list(adonis.output = ref40.denose.adonis, adonis.R.CI.bootrap = ref40.boos.ci.adonis, adonis.R.CI.F = ref40.f.ci.adonis, betadisper.output=ref40.betadisper,betadisper.HSD = ref40.betadisper.HSD,anosim = ref40.anosim)
save(z = list.adonis.ref40, file = "pernomova/ref40.output.RData",row.names = T)
ref40.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)

####ref60
####ref60
####ref60
####ref60
####ref60
ref60 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/silva138/level-6.csv",header = T,check.names = F)
rownames(ref60) <- ref60$index
colnames(ref60)
ref60 <- ref60[,-c(1,294:315)]
rowSums(ref60)

ref60.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/ref60.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref60.meta),paste("Ref60.",rownames(ref60),sep="")) != c(1:dim(ref60)[1]))
rownames(ref60)<- paste("Ref60.",rownames(ref60),sep="")
####ref60
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref60.meta),paste("Ref60.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref60.meta.rn <- ref60.meta
rownames(ref60.meta.rn) <- gsub("Ref60.","",rownames(ref60.meta))

set.seed(1)
colnames(ref60.meta)

ref60.denose.adonis <- c()
ref60.boos.ci.adonis <- c()
ref60.f.ci.adonis <- c()
ref60.betadisper <- c()
ref60.betadisper.HSD <- c()
ref60.anosim <- c()
#depression
i <- 5

for(dist.file in file){
    rm.na.id <- !is.na(ref60.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref60.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref60.betadisper.HSD <- rbind(ref60.betadisper.HSD,ps.disper.HSD)
    
    ref60.betadisper <- rbind(ref60.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref60.denose.adonis <- rbind(ref60.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref60.boos.ci.adonis <- rbind(ref60.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref60.f.ci.adonis <- rbind(ref60.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref60.anosim <- rbind(ref60.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref60.denose.adonis) <- rownames(ref60.boos.ci.adonis) <- rownames(ref60.f.ci.adonis) <- rownames(ref60.betadisper) <-rownames(ref60.anosim) <- file

rownames(ref60.betadisper.HSD) <- paste(rownames(ref60.betadisper.HSD),file,sep="_")

list.adonis.ref60 <- list(adonis.output = ref60.denose.adonis, adonis.R.CI.bootrap = ref60.boos.ci.adonis, adonis.R.CI.F = ref60.f.ci.adonis, betadisper.output=ref60.betadisper,betadisper.HSD = ref60.betadisper.HSD,anosim = ref60.anosim)
save(z = list.adonis.ref60, file = "pernomova/ref60.output.RData",row.names = T)
ref60.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)

####ref65
####ref65
####ref65
####ref65
####ref65
ref65 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/silva138/level-6.csv",header = T,check.names = F)
rownames(ref65) <- ref65$index
colnames(ref65)
ref65 <- ref65[,-c(1,474:488)]
rowSums(ref65)

ref65.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/ref65.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref65.meta),paste("Ref65.",rownames(ref65),sep="")) != c(1:dim(ref65)[1]))
pid <- pmatch(rownames(ref65.meta),paste("Ref65.",rownames(ref65),sep=""))
ref65 <- ref65[pid,]
sum(pmatch(rownames(ref65.meta),paste("Ref65.",rownames(ref65),sep="")) != c(1:dim(ref65)[1]))
rownames(ref65) <- paste("Ref65.",rownames(ref65),sep="")

####ref65
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref65.meta),paste("Ref65.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref65.meta.rn <- ref65.meta
rownames(ref65.meta.rn) <- gsub("Ref65.","",rownames(ref65.meta))

set.seed(1)
colnames(ref65.meta)

ref65.denose.adonis <- c()
ref65.boos.ci.adonis <- c()
ref65.f.ci.adonis <- c()
ref65.betadisper <- c()
ref65.betadisper.HSD <- c()
ref65.anosim <- c()
#depression
i <- 8

for(dist.file in file){
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    sample.id <- intersect(labels(unw.dis),rownames(ref65.meta.rn))
    ref65.meta.rn<-ref65.meta.rn[sample.id,]
    dat.sub <- dat[sample.id,sample.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref65.meta.rn[sample.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref65.betadisper.HSD <- rbind(ref65.betadisper.HSD,ps.disper.HSD)
    ref65.betadisper <- rbind(ref65.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref65.denose.adonis <- rbind(ref65.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref65.boos.ci.adonis <- rbind(ref65.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref65.f.ci.adonis <- rbind(ref65.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref65.anosim <- rbind(ref65.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref65.denose.adonis) <- rownames(ref65.boos.ci.adonis) <- rownames(ref65.f.ci.adonis) <- rownames(ref65.betadisper) <-rownames(ref65.anosim) <- file

rownames(ref65.betadisper.HSD) <- paste(rownames(ref65.betadisper.HSD),file,sep="_")

list.adonis.ref65 <- list(adonis.output = ref65.denose.adonis, adonis.R.CI.bootrap = ref65.boos.ci.adonis, adonis.R.CI.F = ref65.f.ci.adonis, betadisper.output=ref65.betadisper,betadisper.HSD = ref65.betadisper.HSD ,anosim = ref65.anosim)
save(z = list.adonis.ref65, file = "pernomova/ref65.output.RData",row.names = T)
ref65.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)
####ref77
####ref77
####ref77
####ref77
####ref77
ref77 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/sliva138/level-6.csv",header = T,check.names = F)
rownames(ref77) <- ref77$index
colnames(ref77)
ref77 <- ref77[,-c(1,300:303)]
rowSums(ref77)

ref77.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/ref77.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref77.meta),paste("Ref77.",rownames(ref77),sep="")) != c(1:dim(ref77)[1]))
rownames(ref77) <- paste("Ref77.",rownames(ref77),sep="")

####ref77
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref77.meta),paste("Ref77.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref77.meta.rn <- ref77.meta
rownames(ref77.meta.rn) <- gsub("Ref77.","",rownames(ref77.meta))

set.seed(1)
colnames(ref77.meta)

ref77.denose.adonis <- c()
ref77.boos.ci.adonis <- c()
ref77.f.ci.adonis <- c()
ref77.betadisper <- c()
ref77.betadisper.HSD <- c()
ref77.anosim <- c()
#depression
i <- 5

for(dist.file in file){
    rm.na.id <- !is.na(ref77.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref77.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref77.betadisper.HSD <- rbind(ref77.betadisper.HSD,ps.disper.HSD)
    ref77.betadisper <- rbind(ref77.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref77.denose.adonis <- rbind(ref77.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref77.boos.ci.adonis <- rbind(ref77.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref77.f.ci.adonis <- rbind(ref77.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref77.anosim <- rbind(ref77.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref77.denose.adonis) <- rownames(ref77.boos.ci.adonis) <- rownames(ref77.f.ci.adonis) <- rownames(ref77.betadisper) <- rownames(ref77.anosim) <- file

rownames(ref77.betadisper.HSD) <- paste(rownames(ref77.betadisper.HSD),file,sep="_")

list.adonis.ref77 <- list(adonis.output = ref77.denose.adonis, adonis.R.CI.bootrap = ref77.boos.ci.adonis, adonis.R.CI.F = ref77.f.ci.adonis, betadisper.output=ref77.betadisper,betadisper.HSD = ref77.betadisper.HSD,anosim = ref77.anosim)
save(z = list.adonis.ref77, file = "pernomova/ref77.output.RData",row.names = T)
ref77.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)
####ref79
####ref79
####ref79
####ref79
####ref79
ref79 <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/silva138/level-6.csv",header = T,check.names = F)
rownames(ref79) <- ref79$index
colnames(ref79)
ref79 <- ref79[,-c(1,174:175)]
rowSums(ref79)

ref79.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/ref79.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)

sum(pmatch(rownames(ref79.meta),paste("Ref79.",rownames(ref79),sep="")) != c(1:dim(ref79)[1]))
rownames(ref79) <- paste("Ref79.",rownames(ref79),sep="")

####ref79
unw.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/dist.qza/unweighted_unifrac_distance_matrix.qza")
unw.dis <- unw.dis.qza$data
wuni.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/dist.qza/weighted_unifrac_distance_matrix.qza")
wuni.dis <- wuni.dis.qza$data

jsd.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/dist.qza/jaccard_distance_matrix.qza")
jsd.dis <- jsd.dis.qza$data

bray.dis.qza <- read_qza("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/dada2/dist.qza/bray_curtis_distance_matrix.qza")
bray.dis <- bray.dis.qza$data


sum(pmatch(rownames(ref79.meta),paste("Ref79.",label(unw.dis),sep="")) != c(1:length(labels(unw.dis))))
ref79.meta.rn <- ref79.meta
rownames(ref79.meta.rn) <- gsub("Ref79.","",rownames(ref79.meta))

set.seed(1)
colnames(ref79.meta)

ref79.denose.adonis <- c()
ref79.boos.ci.adonis <- c()
ref79.f.ci.adonis <- c()
ref79.betadisper <- c()
ref79.betadisper.HSD <- c()
ref79.anosim <- c()
#depression
i <- 1

for(dist.file in file){
    rm.na.id <- !is.na(ref79.meta.rn[,i])
    txt <- paste("dat <- as.matrix(",dist.file,")",sep="")
    eval((parse(text =txt)))
    dat.sub <- dat[rm.na.id,rm.na.id]
    dat.sub <- as.dist(dat.sub)
    metadate.rm.na <- ref79.meta.rn[rm.na.id,i,drop=F]
    ps.disper <- betadisper(dat.sub, metadate.rm.na$Depression)
    beta.disp <- vegan::permutest(ps.disper)
    ps.disper.HSD <- TukeyHSD(ps.disper)$group
    
    ref79.betadisper.HSD <- rbind(ref79.betadisper.HSD,ps.disper.HSD)
    ref79.betadisper <- rbind(ref79.betadisper,beta.disp$tab[1,])
    
    ####betadisper plot
    df <- data.frame(Distance_to_centroid=ps.disper$distances,Group=ps.disper$group)
    groups <- ps.disper$group
    
    my_comparisons <- list(c("0", "1"))
    df$Group <- factor(df$Group,levels = c(1,0))
    p <- ggboxplot(df, x = "Group", y = "Distance_to_centroid",
                   fill = "Group", palette = "jco")+ 
        stat_compare_means(comparisons = my_comparisons)+xlab("Depression")
    txt <- paste("p.",dist.file," <- p",sep="")
    eval((parse(text =txt)))
    
    fit_adonis_before.denozed.adjust <- adonis2(dat.sub ~ ., data = metadate.rm.na)
    
    adonis.out <- function (data,dist,indices){
        dist.matrix <- as.matrix(dist)
        dist.matrix.sub <- dist.matrix[indices,indices]
        dist.matrix.sub <- as.dist(dist.matrix.sub)
        colnames(data) <- "fac"
        data <- data[indices,,drop=F]
        data[,1] <- as.factor(data[,1] )
        return(adonis2(dist.matrix.sub ~ fac, data = data)$R2[1])
    }
    
    adonis.b.out <- boot(data=metadate.rm.na, dist =dat.sub ,statistic=adonis.out, 
                         R=1000)
    
    ref79.denose.adonis <- rbind(ref79.denose.adonis,cbind(as.matrix(fit_adonis_before.denozed.adjust)[1,,drop=F]))
    ref79.boos.ci.adonis <- rbind(ref79.boos.ci.adonis,c(adonis.b.out$t0,quantile(adonis.b.out$t,c(0.025,0.975))))
    
    f.value <- fit_adonis_before.denozed.adjust$F[1]
    df.1 <- fit_adonis_before.denozed.adjust$Df[1]
    df.2 <-  fit_adonis_before.denozed.adjust$Df[2]
    Delta <- unlist(MBESS::conf.limits.ncf(
        F.value = f.value, 
        df.1, df.2, 
        conf.level = 0.95)[c("Lower.Limit", "Upper.Limit")])
    out <- c(Delta / (Delta + df.1 + df.2 + 1))
    out[is.na(out)] <- 0
    ref79.f.ci.adonis <- rbind(ref79.f.ci.adonis,out)
    
    pathotype.anosim <- anosim(dat.sub, metadate.rm.na$Depression)
    ref79.anosim <- rbind(ref79.anosim,c(pathotype.anosim$statistic,pathotype.anosim$signif))
}

rownames(ref79.denose.adonis) <- rownames(ref79.boos.ci.adonis) <- rownames(ref79.f.ci.adonis) <- rownames(ref79.betadisper) <- rownames(ref79.anosim) <- file

rownames(ref79.betadisper.HSD) <- paste(rownames(ref79.betadisper.HSD),file,sep="_")

list.adonis.ref79 <- list(adonis.output = ref79.denose.adonis, adonis.R.CI.bootrap = ref79.boos.ci.adonis, adonis.R.CI.F = ref79.f.ci.adonis, betadisper.output=ref79.betadisper,betadisper.HSD = ref79.betadisper.HSD,anosim = ref79.anosim)

save(z = list.adonis.ref79, file = "pernomova/ref79.output.RData",row.names = T)

ref79.fig <- list(Bray = p.bray.dis, JSD = p.jsd.dis, WeightUnifrac = p.wuni.dis, Unweigh = p.unw.dis)


save(z = ref1.fig, file = "pernomova/ref1.fig.RData",row.names = T)
save(z = ref23.fig, file = "pernomova/ref23.fig.RData",row.names = T)
save(z = ref27.fig, file = "pernomova/ref27.fig.RData",row.names = T)
save(z = ref40.fig, file = "pernomova/ref40.fig.RData",row.names = T)
save(z = ref60.fig, file = "pernomova/ref60.fig.RData",row.names = T)
save(z = ref65.fig, file = "pernomova/ref65.fig.RData",row.names = T)
save(z = ref77.fig, file = "pernomova/ref77.fig.RData",row.names = T)
save(z = ref79.fig, file = "pernomova/ref79.fig.RData",row.names = T)

#######load the output data
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/")
load("pernomova/ref1.output.RData")
load("pernomova/ref23.output.RData")
load("pernomova/ref27.output.RData")
load("pernomova/ref40.output.RData")
load("pernomova/ref60.output.RData")
load("pernomova/ref65.output.RData")
load("pernomova/ref77.output.RData")
load("pernomova/ref79.output.RData")

adonis.list.file <- ls(pattern="list.adonis.ref")
adonis.estimat.bray <- c()
adonis.estimat.jsd <- c()
adonis.estimat.unw <- c()
adonis.estimat.wuni <- c()
for(i in adonis.list.file){
    txt <- paste("estimat.list <-", i,sep="")
    eval((parse(text =txt)))
    
    adonis.estimat.bray <- rbind(adonis.estimat.bray,
                                 c(estimat.list$adonis.output[1,3],
                                   estimat.list$adonis.R.CI.bootrap[1,1],
                                   estimat.list$adonis.R.CI.bootrap[1,2],
                                   estimat.list$adonis.R.CI.bootrap[1,3],
                                   estimat.list$adonis.output[1,5]))
    adonis.estimat.jsd <- rbind(adonis.estimat.jsd,
                                c(estimat.list$adonis.output[2,3],
                                  estimat.list$adonis.R.CI.bootrap[2,1],
                                  estimat.list$adonis.R.CI.bootrap[2,2],
                                  estimat.list$adonis.R.CI.bootrap[2,3],
                                  estimat.list$adonis.output[2,5]))
    adonis.estimat.unw <- rbind(adonis.estimat.unw,
                                c(estimat.list$adonis.output[3,3],
                                  estimat.list$adonis.R.CI.bootrap[3,1],
                                  estimat.list$adonis.R.CI.bootrap[3,2],
                                  estimat.list$adonis.R.CI.bootrap[3,3],
                                  estimat.list$adonis.output[3,5]))
    adonis.estimat.wuni <- rbind(adonis.estimat.wuni,
                                 c(estimat.list$adonis.output[4,3],
                                   estimat.list$adonis.R.CI.bootrap[4,1],
                                   estimat.list$adonis.R.CI.bootrap[4,2],
                                   estimat.list$adonis.R.CI.bootrap[4,3],
                                   estimat.list$adonis.output[4,5]))
    
}
colnames(adonis.estimat.bray) <-colnames(adonis.estimat.jsd) <-colnames(adonis.estimat.unw) <-colnames(adonis.estimat.wuni) <-c("R2","boostrap_r2","boostrap_2.5%","boostrap_97.5%","pvalue")
id_title <- read.table("ID2Title_meta_analysis_study.txt",sep="\t",head=T)
###bray distance
adonis.estimat.bray <- as.data.frame(adonis.estimat.bray)
adonis.estimat.bray$logR <- log(adonis.estimat.bray[,1])
adonis.estimat.bray$se.logR <- (log(adonis.estimat.bray[,4])-
                                    log(adonis.estimat.bray[,3]))/(2*1.96)
cbind(adonis.list.file,id_title)
adonis.estimat.bray$study <- id_title$title

library(metafor)
bray.fe <- rma(yi=logR, sei=se.logR, slab=study,
               method="FE", data=adonis.estimat.bray)
with(bray.fe, exp(c(b, ci.lb, ci.ub)))

bray.re <- rma(yi=logR, sei=se.logR, slab=study,
               method="REML", data=adonis.estimat.bray)
with(bray.re, exp(c(b, ci.lb, ci.ub)))
library("meta")
pdf("pernomova/bray.R2.model.pdf",width = 10,height = 5)
bray.re.metafor <- metagen(TE=logR, seTE=se.logR, studlab=study, method.tau="REML",method.tau.ci="QP", sm="OR",data=adonis.estimat.bray)
metafor::forest(bray.re.metafor,transf=exp, refline=0.1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Bray–Curtis dissimilarity")
dev.off()
###jsd
adonis.estimat.jsd <- as.data.frame(adonis.estimat.jsd)
adonis.estimat.jsd$logR <- log(adonis.estimat.jsd[,1])
adonis.estimat.jsd$se.logR <- (log(adonis.estimat.jsd[,4])-
                                   log(adonis.estimat.jsd[,3]))/(2*1.96)
adonis.estimat.jsd$study <- id_title$title
library(metafor)
jsd.fe <- rma(yi=logR, sei=se.logR, slab=study,
              method="FE", data=adonis.estimat.jsd)
with(jsd.fe, exp(c(b, ci.lb, ci.ub)))

jsd.re <- rma(yi=logR, sei=se.logR, slab=study,
              method="REML", data=adonis.estimat.jsd)
with(jsd.re, exp(c(b, ci.lb, ci.ub)))


pdf("pernomova/jsd.R2.model.pdf",width = 10,height = 5)
jsd.re.metafor <- metagen(TE=logR, seTE=se.logR, studlab=study, method.tau="REML",method.tau.ci="QP", sm="OR",data=adonis.estimat.jsd)
metafor::forest(jsd.re.metafor,transf=exp, refline=0.1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Jensen–Shannon divergence")
dev.off()
###unweight unfrac
adonis.estimat.unw <- as.data.frame(adonis.estimat.unw)
adonis.estimat.unw$logR <- log(adonis.estimat.unw[,1])
adonis.estimat.unw$se.logR <- (log(adonis.estimat.unw[,4])-
                                   log(adonis.estimat.unw[,3]))/(2*1.96)
adonis.estimat.unw$study <- id_title$title
library(metafor)
unw.fe <- rma(yi=logR, sei=se.logR, slab=study,
              method="FE", data=adonis.estimat.unw)
with(unw.fe, exp(c(b, ci.lb, ci.ub)))

unw.re <- rma(yi=logR, sei=se.logR, slab=study,
              method="REML", data=adonis.estimat.unw)
with(unw.re, exp(c(b, ci.lb, ci.ub)))

pdf("pernomova/unw.R2.model.pdf",width = 10,height = 5)
unw.re.metafor <- metagen(TE=logR, seTE=se.logR, studlab=study, method.tau="REML",method.tau.ci="QP", sm="OR",data=adonis.estimat.unw)
metafor::forest(unw.re.metafor,transf=exp, refline=0.1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Unweighted UniFrac distance")
dev.off()
####weight Unifrac
adonis.estimat.wuni <- as.data.frame(adonis.estimat.wuni)
adonis.estimat.wuni$logR <- log(adonis.estimat.wuni[,1])
adonis.estimat.wuni$se.logR <- (log(adonis.estimat.wuni[,4])-
                                    log(adonis.estimat.wuni[,3]))/(2*1.96)
adonis.estimat.wuni$study <- id_title$title
library(metafor)
wuni.fe <- rma(yi=logR, sei=se.logR, slab=study,
               method="FE", data=adonis.estimat.wuni)
with(wuni.fe, exp(c(b, ci.lb, ci.ub)))

wuni.re <- rma(yi=logR, sei=se.logR, slab=study,
               method="REML", data=adonis.estimat.wuni)
with(wuni.re, exp(c(b, ci.lb, ci.ub)))
pdf("pernomova/wuf.R2.model.pdf",width = 10,height = 5)
wuni.re.metafor <- metagen(TE=logR, seTE=se.logR, studlab=study, method.tau="REML",method.tau.ci="QP", sm="OR",data=adonis.estimat.wuni)
metafor::forest(wuni.re.metafor,transf=exp, refline=0.1,digits=4,  print.stat = TRUE,backtransf=TRUE,xlab="Weighted UniFrac distance")
dev.off()
###code from meta package
##sm == "PLOGIT":Logit transformation
#TE <- log((event + incr.event)/(n - event + incr.event))
#seTE <- sqrt(1/(event + incr.event) + 1/((n - event +                                                       incr.event)))
#transf.null.effect <- log(null.effect/(1 - null.effect))
##
# TE <- log(adonis.estimat.wuni[,1])
# seTE <- sqrt(1/event + 1/(n - event)) ####cannot calculated the seTE only by R^2
# m <- metagen(TE, seTE, studlab)

pdf("pernomova/adonise.forestplot.pdf",width = 10,height = 6)
par(mfrow=c(2,2),mar=c(2,2,2,2))
forest(bray.re, transf=exp, refline=0,xlab = "Bray–Curtis dissimilarity")
forest(jsd.re, transf=exp, refline=0,xlab = "Jensen–Shannon divergence")
forest(unw.re, transf=exp, refline=0,xlab = "Unweighted UniFrac distance")
forest(wuni.re, transf=exp, refline=0,xlab="Weighted UniFrac distance")

dev.off()

####Heatmap
#permonova heatmap

dat.R2 <- cbind(adonis.estimat.bray[,1],adonis.estimat.jsd[,1],
                adonis.estimat.unw[,1],adonis.estimat.wuni[,1] )
dat.p <- cbind(adonis.estimat.bray[,5],adonis.estimat.jsd[,5],
               adonis.estimat.unw[,5],adonis.estimat.wuni[,5] )
colnames(dat.R2) <- colnames(dat.p) <- c("Bray–Curtis dissimilarity","Jensen–Shannon divergence","Unweighted UniFrac distance","Weighted UniFrac distance")
rownames(dat.R2) <- rownames(dat.p) <- adonis.estimat.bray[,8]
dat.p.sig <-dat.p
dat.p.sig[dat.p<0.05] <- "*"
dat.p.sig[dat.p>0.05] <- ""
pdf("pernomova/permanova.heatmap.pdf")
pheatmap::pheatmap(mat = dat.R2,display_numbers = dat.p.sig)
dev.off()
#Beta-Dispersion heatmap

betadisper.estimat.bray <- c()
betadisper.estimat.jsd <- c()
betadisper.estimat.unw <- c()
betadisper.estimat.wuni <- c()
for(i in adonis.list.file){
    txt <- paste("estimat.list <-", i,sep="")
    eval((parse(text =txt)))
    
    betadisper.estimat.bray <- rbind(betadisper.estimat.bray,
                                     c(estimat.list$betadisper.HSD[1,1],
                                       estimat.list$betadisper.HSD[1,2],
                                       estimat.list$betadisper.HSD[1,3],
                                       estimat.list$betadisper.HSD[1,4]))
    betadisper.estimat.jsd <- rbind(betadisper.estimat.jsd,
                                    c(estimat.list$betadisper.HSD[2,1],
                                      estimat.list$betadisper.HSD[2,2],
                                      estimat.list$betadisper.HSD[2,3],
                                      estimat.list$betadisper.HSD[2,4]))
    betadisper.estimat.unw <- rbind(betadisper.estimat.unw,
                                    c(estimat.list$betadisper.HSD[3,1],
                                      estimat.list$betadisper.HSD[3,2],
                                      estimat.list$betadisper.HSD[3,3],
                                      estimat.list$betadisper.HSD[3,4]))
    betadisper.estimat.wuni <- rbind(betadisper.estimat.wuni,
                                     c(estimat.list$betadisper.HSD[4,1],
                                       estimat.list$betadisper.HSD[4,2],
                                       estimat.list$betadisper.HSD[4,3],
                                       estimat.list$betadisper.HSD[4,4]))
    
}
colnames(betadisper.estimat.bray) <-colnames(betadisper.estimat.jsd) <-colnames(betadisper.estimat.unw) <-colnames(betadisper.estimat.wuni) <-c("diff1-0","lwr","upr","pvalue")

##heatmap

dat.diff <- cbind(betadisper.estimat.bray[,1],betadisper.estimat.jsd[,1],
                  betadisper.estimat.unw[,1],betadisper.estimat.wuni[,1] )
dat.diff.p <- cbind(betadisper.estimat.bray[,4],betadisper.estimat.jsd[,4],
               betadisper.estimat.unw[,4],betadisper.estimat.wuni[,4] )
colnames(dat.diff) <- colnames(dat.p) <- c("Bray–Curtis dissimilarity","Jensen–Shannon divergence","Unweighted UniFrac distance","Weighted UniFrac distance")
rownames(dat.diff) <- rownames(dat.p) <- id_title$title

dat.p.sig <-dat.diff.p
dat.p.sig[dat.diff.p<0.05] <- "*"
dat.p.sig[dat.diff.p>0.05] <- ""
pdf("pernomova/betadisper.heatmap.pdf")
pheatmap::pheatmap(mat = dat.diff,display_numbers = dat.p.sig)
dev.off()

