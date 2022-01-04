###Review 27 dada
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
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/")

###marker in systematic review
decrease <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/decrease.maker.list.txt",sep="\t")
increase <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/increase.maker.list.txt",sep="\t")
##metadata
metadata <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/DataDownload/ref27/metadata.txt",sep="\t",row.names = 1,head=T)
metadata$group_time <- paste(metadata$Group,metadata$Time,sep="_")
###Among the 52 BD patients, fecal samples were re‐examined in 20 of them after quetiapine treatment. The HDRS‐17 and Montgomery‐Åsberg Depression Rating Scale (MADRS) scores were significantly decreased following quetiapine treatment in BD patients (P < 0.05)
###comparision bettween BD_1 vs H_1 and BD_1 vs BD_2


#######match id
#######match id
metadata.BD <- metadata[metadata$Group == "BD",]
metadata.BD$subject <- rownames(metadata.BD)
metadata.BD$subject <-  gsub("-2","",metadata.BD$subject )
metadata.BD$subject <-  gsub("-","",metadata.BD$subject )

id <- table(metadata.BD$subject)
match.id <- names(id[id==2])
length(match.id)
t1<- metadata.BD[metadata.BD$subject %in% match.id & metadata.BD$Time ==1,]
t2 <- metadata.BD[metadata.BD$subject %in% match.id & metadata.BD$Time ==2,]

t1 <- t1[order(t1$subject),]
t2 <- t2[order(t2$subject),]
sum(pmatch(t1$subject,t2$subject) == c(1:20))

id.t1 <- rownames(t1)
id.t2 <- rownames(t2)

######################
######################
##1.alpha diversity

#observe otu
obt <- read.table("dada2/not_cb/observed_otus_vector.qza/2a9ff6c1-c987-48c9-99ad-56b411920be8/data/alpha-diversity.tsv",check.names = T)

#chao1
chao1.d <- read.table("dada2/not_cb/chao1_alpha_diverisyt.qza/e711455b-9477-49e5-a31b-c7f8a2c611a5/data/alpha-diversity.tsv",check.names = T)

#shannon
shannon <- read.table("dada2/not_cb/shannon_vector.qza/d337d6a4-b381-4fdd-84d7-d28be10be458/data/alpha-diversity.tsv",check.names = T)

#simpson
simpson <- read.table("dada2/not_cb/simpson_alpha_diverisyt.qza/6511ddec-ab24-4f94-a8c7-2888c6bf8cd1/data/alpha-diversity.tsv",check.names = T)




pid <- pmatch(rownames(chao1.d),rownames(metadata))
metadata <- metadata[pid,]

alpha <- cbind.data.frame(obt,chao1.d,simpson,shannon,metadata)
colnames(alpha)
write.table(alpha,"Ref27_aloha_diversity.txt",sep="\t",quote = F,col.names = T,row.names = T)
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

###simpsom index significant different with normal distribution 
#####alpha diversity test
wil.res <- c()
t.res <- c()

wil.res.t <- c()
t.res.t <- c()
for(i in 1:4){
    #BD_1 vs H
    dat1 <-  alpha[alpha$group_time=="BD_1" | alpha$group_time=="H_1",]
    p.t <- t.test(dat1[,i] ~ dat1[,"group_time"])$p.value
    p.v <- wilcox.test(dat1[,i] ~ dat1[,"group_time"])$p.value
    means <- as.matrix(by(dat1[,i] , dat1[,"group_time"],mean))[2:1,]
    stand.d <- as.matrix(by(dat1[,i] , dat1[,"group_time"],sd))[2:1,]
    t.res<-rbind(t.res,c(p.t,means,stand.d))

    iqr <- as.matrix(by(dat1[,i],dat1[,"group_time"],quantile))
    wil.res<-rbind(wil.res,c(p.v,iqr[[2]][c(2,3,4)],iqr[[1]][c(2,3,4)]))
    
    #BD_1 vs BD_2
    dat2.1 <- alpha[pmatch(id.t1,rownames(alpha)),] 
    dat2.2 <- alpha[pmatch(id.t2,rownames(alpha)),] 
    p.t.2 <- t.test(dat2.1[,i] , dat2.2[,i],paired = T)$p.value
    p.v.2 <- wilcox.test(dat2.1[,i] , dat2.2[,i],paired = T)$p.value
    
    means1 <- mean(dat2.1[,i])
    means2 <- mean(dat2.2[,i])
    stand.d1 <- sd(dat2.1[,i])
    stand.d2 <- sd(dat2.2[,i])
    t.res.t<-rbind(t.res.t,c(p.t.2,means2,means1,stand.d2,stand.d1))
    
    iqr.1 <- quantile(dat2.1[,i])[c(2,3,4)]
    iqr.2 <- quantile(dat2.2[,i])[c(2,3,4)]
    ##BD_treat is less depression compared the BD_t1
    wil.res.t<-rbind(wil.res.t,c(p.v.2,iqr.2,iqr.1))
    
}

rownames(wil.res) <- rownames(t.res) <- colnames(alpha)[1:4]
rownames(wil.res.t) <- rownames(t.res.t) <- colnames(alpha)[1:4]

colnames(t.res) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")

colnames(t.res.t) <-  c("P.value","mean in non-depression","mean in depression","sd in non-depression","sd in depression")
colnames(wil.res.t) <-  c("P.value","25% IQR in non-depression","median in non-depression","75% IQR in non-depression","25% IQR in depression","median in depression","75% IQR in depression")

wil.res[wil.res[,1]<0.05,]
t.res[t.res[,1]<0.05,]

wil.res.t[wil.res.t[,1]<0.05,]
t.res.t[t.res.t[,1]<0.05,]

library("kableExtra")


wil.res[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (BD_T1 vs Health) wilcox test") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res.t[,c(1,3,6)] %>% as.data.frame() %>% round(2) %>%
    kbl(caption = "Alpha diversity (BD vs BD_Treat) ") %>%
    kable_classic(full_width = F, html_font = "Cambria") 

wil.res <- cbind(wil.res,dim(dat1)[1],table(dat1$group_time)[2],table(dat1$group_time)[1])
colnames(wil.res)[8:10] <- c("Total sample size","non-depression","depression")

wil.res.t <- cbind(wil.res.t,dim(dat2.2)[1]+dim(dat2.1)[1],dim(dat2.2)[1],dim(dat2.1)[1])
colnames(wil.res.t)[8:10] <- c("Total sample size","non-depression","depression")


t.res.t <- cbind(t.res.t,dim(dat2.2)[1]+dim(dat2.1)[1],dim(dat2.2)[1],dim(dat2.1)[1])
colnames(t.res.t)[6:8] <- c("Total sample size","non-depression","depression")

t.res <- cbind(t.res,dim(dat1)[1],table(dat1$group_time)[2],table(dat1$group_time)[1])
colnames(t.res)[6:8] <- c("Total sample size","non-depression","depression")


write.table(t.res,"ref27_dada_alpha.diversity.BD_H.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res,"ref27_dada_alpha.diversity.BD_H.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

write.table(t.res.t,"ref27_dada_alpha.diversity.BD1_BD2.t.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(wil.res.t,"ref27_dada_alpha.diversity.BD1_BD2.wilcoxon.test.result.txt",sep="\t",quote = F,col.names = T,row.names = T)

###########alpha divesity odd ratio
alpha.odd <- alpha
write.table(alpha.odd,"alpha_diverisy.rawdata.txt",quote = F,sep = "\t")
odd.crude.res_T<- c()
odd.crude.res.BD_H <- c()
odd.gee.res <- c()

library(survival)
library(geepack)
library(doBy)
for(i in 1:4){
    #BD_1 vs H
    dat1 <-  alpha.odd[alpha.odd$group_time=="BD_1" | alpha.odd$group_time=="H_1",]
    dat1$Depression <- 1
    dat1$Depression[dat1$Group=="H"] <- 0
    ##Group
    data <- cbind(dat1[,i,drop=F],dat1[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    reg <- glm(Depression ~alpha, data=data, family=binomial)
    odd.crude.res.BD_H <- rbind(odd.crude.res.BD_H,
                           #round(exp(cbind(coef(reg), confint(reg))),3)[2,])
                           c(cbind(coef(reg), confint(reg))[2,],"BDH_raw_logOR",colnames(alpha.odd)[i]))
    
    #BD_1 vs BD_2
    dat2.1 <- alpha[pmatch(id.t1,rownames(alpha)),] 
    dat2.2 <- alpha[pmatch(id.t2,rownames(alpha)),] 
    dat2 <- rbind(dat2.1,dat2.2)
    dat2$Depression <- 1
    dat2$Depression[dat2$Time==2] <- 0
    dat2$subject <- rownames(dat2)
    dat2$subject <-  gsub("-2","",dat2$subject )
    dat2$subject <-  gsub("-","",dat2$subject )
    data <- cbind(dat2[,i,drop=F],dat2[,c("Depression","subject"),drop=F])
    colnames(data)[1] <- "alpha"
    data$subject <- as.factor(data$subject)
    reg2 <- glm(Depression ~alpha, data=data, family=binomial)
    odd.crude.res_T <- rbind(odd.crude.res_T,
                           #round(exp(cbind(coef(reg), confint(reg))),3)[2,])
                           c(cbind(coef(reg2), confint(reg2))[2,],"T_raw_LogOR",colnames(alpha.odd)[i]))
    gee.adjust <- geeglm(Depression ~alpha, family=binomial, data=data, 
                         id=subject, corstr = "independence")
    odd.gee.res <- rbind(odd.gee.res,
                         #exp(c(esticon(gee.adjust, c(0,1))$estimate,esticon(gee.adjust,c(0,1))$lwr,esticon(gee.adjust, c(0,1))$upr)))
                         c(esticon(gee.adjust, c(0,1))$estimate,esticon(gee.adjust,c(0,1))$lwr,esticon(gee.adjust, c(0,1))$upr,"T_raw_LogOR",colnames(alpha.odd)[i]))
    
    
    ###seperated by median
    ##group
    data <- cbind(dat1[,i,drop=F],dat1[,"Depression",drop=F])
    colnames(data)[1] <- "alpha"
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg <- glm(Depression ~H.alpha, data=data, family=binomial)
    odd.crude.res.BD_H <- rbind(odd.crude.res.BD_H,
                           c(round(exp(cbind(coef(reg), confint(reg))),3)[2,],"BDH_higher_OR",colnames(alpha.odd)[i]))
    ###BD_1 vs BD_2
    data <- cbind(dat2[,i,drop=F],dat2[,c("Depression","subject"),drop=F])
    colnames(data)[1] <- "alpha"
    data$subject <- as.factor(data$subject)
    med.alpha <- median(data$alpha)
    data$H.alpha <- 1
    data$H.alpha[data$alpha < med.alpha] <-0
    reg2 <- glm(Depression ~H.alpha, data=data, family=binomial)
    odd.crude.res_T <- rbind(odd.crude.res_T,
                           c(round(exp(cbind(coef(reg2), confint(reg2))),3)[2,],"T_heigher_OR",colnames(alpha.odd)[i]))
    
    gee.adjust <- geeglm(Depression ~H.alpha, family=binomial, data=data, 
                         id=subject, corstr = "independence")
    odd.gee.res <- rbind(odd.gee.res,
                         c(exp(c(esticon(gee.adjust, c(0,1))$estimate,esticon(gee.adjust,c(0,1))$lwr,esticon(gee.adjust, c(0,1))$upr)),"T_heigher_OR",colnames(alpha.odd)[i]))
    
}
rownames(odd.crude.res_T) <- c(paste(c("raw_t","Higher_t"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.crude.res.BD_H) <- c(paste(c("raw_BD_H","Higher_BD_H"),rep(colnames(alpha)[1:4],each=2)))
rownames(odd.gee.res) <- c(paste(c("raw","Higher"),rep(colnames(alpha)[1:4],each=2)))

odd.crude.res.BD_H
odd.crude.res.BD_H <- cbind(odd.crude.res.BD_H[,1:3],dim(dat1)[1],table(dat1$group_time)[2],table(dat1$group_time)[1])
colnames(odd.crude.res.BD_H)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.crude.res_T
odd.crude.res_T <- cbind(odd.crude.res_T[,1:3],dim(dat2.2)[1]+dim(dat2.1)[1],dim(dat2.2)[1],dim(dat2.1)[1])
colnames(odd.crude.res_T)[4:6] <- c("Total_sample_size","non_depression","depression")

odd.gee.res
odd.gee.res <- cbind(odd.gee.res[,1:3],dim(dat2.2)[1]+dim(dat2.1)[1],dim(dat2.2)[1],dim(dat2.1)[1])
colnames(odd.gee.res)[4:6] <- c("Total_sample_size","non_depression","depression")


write.table(odd.crude.res_T,"Ref27_dada_alpha.diversity.T1T2.OR.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.crude.res.BD_H,"Ref27_dada_alpha.diversity.BDH.OR.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
write.table(odd.gee.res,"Ref27_dada_alpha.diversity.OR_adjust.result.txt",sep="\t",quote = F,col.names = T,row.names = T)
