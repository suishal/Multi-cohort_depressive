###meta-analysis of median
rm(list=ls())
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/alpha_diversity")
ref1 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/Ref1_alpha.diversity.wilcoxon.test.result.txt",sep="\t")
ref23 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref23/Ref23_dada_alpha.diversity.wilcoxon.test.result.txt",sep="\t")
ref27_BD_H <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref27/ref27_dada_alpha.diversity.BD_H.wilcoxon.test.result.txt",sep="\t")
ref27_BD_BD_treat <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref27/ref27_dada_alpha.diversity.BD1_BD2.wilcoxon.test.result.txt",sep="\t")
ref40 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref40/Ref40_dada2_alpha.diversity.wilcoxon.test.result.txt",sep="\t")
ref60 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref60/Ref60_dada2_alpha.diversity.wilcoxon.test.result.txt",sep="\t")
ref65_left <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/Ref65_dada2_alpha.diversity.wilcoxon.Colon_Left.test.result.txt",sep="\t")
ref65_right <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/Ref65_dada2_alpha.diversity.wilcoxon.Colon_Right.test.result.txt",sep="\t")
ref65_Transverse <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/Ref65_dada2_alpha.diversity.wilcoxon.Colon_Transverse.test.result.txt",sep="\t")
ref65_Ileum <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/Ref65_dada2_alpha.diversity.wilcoxon.Ileum.test.result.txt",sep="\t")
ref65_Rectum <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/Ref65_dada2_alpha.diversity.wilcoxon.Rectum.test.result.txt",sep="\t")
ref77 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref77/Ref77_dada2_alpha.diversity.wilcoxon.test.result.txt",sep="\t")
ref79 <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref79/Ref79_dada2_alpha.diversity.wilcoxon.test.result.txt",sep="\t")

#shannono
file <- ls()
observed_otus <-c()
chao1 <-c()
simpson <-c()
shannon <-c()
for(i in file){
    txt <- paste("observed_otus <- rbind(observed_otus,",i,"[1,])",sep="")
    eval((parse(text =txt)))
    
    txt <- paste("chao1 <- rbind(chao1,",i,"[2,])",sep="")
    eval((parse(text =txt)))
    
    txt <- paste("simpson <- rbind(simpson,",i,"[3,])",sep="")
    eval((parse(text =txt)))
    
    txt <- paste("shannon <- rbind(shannon,",i,"[4,])",sep="")
    eval((parse(text =txt)))
}

rownames(observed_otus) <- rownames(chao1) <- rownames(simpson) <- rownames(shannon) <- file

library(metamedian)

alpha <- c("observed_otus","chao1","simpson","shannon")

##observed_otus
med.g1 <- observed_otus$median.in.non.depression
q1.g1 <-  observed_otus$X25..IQR.in.non.depression 
q3.g1 <- observed_otus$X75..IQR.in.non.depression 
n.g1 <- observed_otus$non.depression

med.g2 <- observed_otus$median.in.depression
q1.g2 <-  observed_otus$X25..IQR.in.depression
q3.g2 <- observed_otus$X75..IQR.in.depression
n.g2 <- observed_otus$depression

qe.obs <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                         q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.obs)

effect.size <- qe.obs$yi
QE.lb <- qe.obs$yi - qnorm(0.975) * sqrt(qe.obs$vi)
QE.ub <-  qe.obs$yi + qnorm(0.975) * sqrt(qe.obs$vi)


QE <- paste0(as.character(sprintf("%.2f", round(effect.size, 2))), 
                                  " [", as.character(sprintf("%.2f", round(QE.lb, 2))), ", ", 
                                  as.character(sprintf("%.2f", round(QE.ub, 2))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.2f", round(qe.obs$b, 2))), 
                         " [", 
                         as.character(sprintf("%.2f", round(qe.obs$ci.lb, 2))), 
                         ", ",  
                         as.character(sprintf("%.2f", round(qe.obs$ci.ub, 2))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Observed OTU",rownames(observed_otus),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.obs$b)
lb.vals.medians <- c(NA, QE.lb, qe.obs$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.obs$ci.ub)

library("forestplot")                     
p1 <- forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, length(effect.size)), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "observed_otus", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
                      hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.c <- cbind.data.frame(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table <- cb.table.c


for(i in 1:13){
    med.g1 <- observed_otus$median.in.non.depression[-i]
    q1.g1 <-  observed_otus$X25..IQR.in.non.depression[-i]
    q3.g1 <- observed_otus$X75..IQR.in.non.depression[-i] 
    n.g1 <- observed_otus$non.depression[-i]
    
    med.g2 <- observed_otus$median.in.depression[-i]
    q1.g2 <-  observed_otus$X25..IQR.in.depression[-i]
    q3.g2 <- observed_otus$X75..IQR.in.depression[-i]
    n.g2 <- observed_otus$depression[-i]
    
    qe.obs <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                 q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
    
    plot(qe.obs)
    
}

##chao1
med.g1 <- chao1$median.in.non.depression
q1.g1 <-  chao1$X25..IQR.in.non.depression 
q3.g1 <- chao1$X75..IQR.in.non.depression 
n.g1 <- chao1$non.depression

med.g2 <- chao1$median.in.depression
q1.g2 <-  chao1$X25..IQR.in.depression
q3.g2 <- chao1$X75..IQR.in.depression
n.g2 <- chao1$depression

qe.chao1 <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
             q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.chao1)


effect.size <- qe.chao1$yi
QE.lb <- qe.chao1$yi - qnorm(0.975) * sqrt(qe.chao1$vi)
QE.ub <-  qe.chao1$yi + qnorm(0.975) * sqrt(qe.chao1$vi)


QE <- paste0(as.character(sprintf("%.2f", round(effect.size, 2))), 
             " [", as.character(sprintf("%.2f", round(QE.lb, 2))), ", ", 
             as.character(sprintf("%.2f", round(QE.ub, 2))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.2f", round(qe.chao1$b, 2))), 
                         " [", 
                         as.character(sprintf("%.2f", round(qe.chao1$ci.lb, 2))), 
                         ", ",  
                         as.character(sprintf("%.2f", round(qe.chao1$ci.ub, 2))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Chao1 index",rownames(chao1),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.chao1$b)
lb.vals.medians <- c(NA, QE.lb, qe.chao1$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.chao1$ci.ub)

par(mfcol=c(2,2))
forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, length(effect.size)), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "chao1", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.c <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table <- rbind.data.frame(cb.table,"",cb.table.c)

for(i in 1:13){
    med.g1 <- chao1$median.in.non.depression[-i]
    q1.g1 <-  chao1$X25..IQR.in.non.depression[-i]
    q3.g1 <- chao1$X75..IQR.in.non.depression[-i] 
    n.g1 <- chao1$non.depression[-i]
    
    med.g2 <- chao1$median.in.depression[-i]
    q1.g2 <-  chao1$X25..IQR.in.depression[-i]
    q3.g2 <- chao1$X75..IQR.in.depression[-i]
    n.g2 <- chao1$depression[-i]
    
    qe.chao1 <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                   q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
    
    plot(qe.chao1)
    
}

##shannon
med.g1 <- shannon$median.in.non.depression
q1.g1 <-  shannon$X25..IQR.in.non.depression
q3.g1 <- shannon$X75..IQR.in.non.depression 
n.g1 <- shannon$non.depression

med.g2 <- shannon$median.in.depression
q1.g2 <-  shannon$X25..IQR.in.depression
q3.g2 <- shannon$X75..IQR.in.depression
n.g2 <- shannon$depression

qe.shannon <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                 q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.shannon)

effect.size <- qe.shannon$yi
QE.lb <- qe.shannon$yi - qnorm(0.975) * sqrt(qe.shannon$vi)
QE.ub <-  qe.shannon$yi + qnorm(0.975) * sqrt(qe.shannon$vi)


QE <- paste0(as.character(sprintf("%.4f", round(effect.size, 4))), 
             " [", as.character(sprintf("%.4f", round(QE.lb, 4))), ", ", 
             as.character(sprintf("%.4f", round(QE.ub, 4))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.4f", round(qe.shannon$b, 4))), 
                         " [", 
                         as.character(sprintf("%.4f", round(qe.shannon$ci.lb, 4))), 
                         ", ",  
                         as.character(sprintf("%.4f", round(qe.shannon$ci.ub, 4))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Shannon index",rownames(shannon),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.shannon$b)
lb.vals.medians <- c(NA, QE.lb, qe.shannon$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.shannon$ci.ub)


forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, length(effect.size)), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "shannon", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.c <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table <- rbind.data.frame(cb.table,"",cb.table.c)


##simpson
med.g1 <- simpson$median.in.non.depression
q1.g1 <-  simpson$X25..IQR.in.non.depression 
q3.g1 <- simpson$X75..IQR.in.non.depression 
n.g1 <- simpson$non.depression

med.g2 <- simpson$median.in.depression
q1.g2 <-  simpson$X25..IQR.in.depression
q3.g2 <- simpson$X75..IQR.in.depression
n.g2 <- simpson$depression

qe.simpson <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                 q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.simpson)

effect.size <- qe.simpson$yi
QE.lb <- qe.simpson$yi - qnorm(0.975) * sqrt(qe.simpson$vi)
QE.ub <-  qe.simpson$yi + qnorm(0.975) * sqrt(qe.simpson$vi)


QE <- paste0(as.character(sprintf("%.4f", round(effect.size, 4))), 
             " [", as.character(sprintf("%.4f", round(QE.lb, 4))), ", ", 
             as.character(sprintf("%.4f", round(QE.ub, 4))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.4f", round(qe.simpson$b, 4))), 
                         " [", 
                         as.character(sprintf("%.4f", round(qe.simpson$ci.lb, 4))), 
                         ", ",  
                         as.character(sprintf("%.4f", round(qe.simpson$ci.ub, 4))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Simpson index",rownames(simpson),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.simpson$b)
lb.vals.medians <- c(NA, QE.lb, qe.simpson$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.simpson$ci.ub)


forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, length(effect.size)), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "simpson", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

####combind
cb.table.c <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table <- rbind.data.frame(cb.table,"",cb.table.c)
id2title <- read.table("../ID2Title_meta_analysis_study.txt",sep="\t",head=T)
id2title$id <- paste("ref",id2title$id,sep= "")

cb.table.rn <- cb.table
for(i in 1:dim(cb.table.rn)[1]){
    if(cb.table.rn$V1 %in% id2title$id){
        pid <- pmatch(cb.table.rn$V1[i],id2title$id)
        cb.table.rn$V1[i] <- id2title$title[pid]
    }
    else{
        id <- unlist(strsplit(cb.table.rn$V1[i],"_"))[1]
        pid <- pmatch(id,id2title$id)
        if(!is.na(pid)){
        cb.table.rn$V1[i] <- gsub(id,id2title$title[pid],cb.table.rn$V1[i])
        }
        else{
            
        }
    }
}
cb.table <- cb.table.rn
pdf("alpha_diversity_raw_median_obs.pdf",height=5)
forestplot(labeltext = cb.table[1:15,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 13), rep(TRUE, 1)), 
           mean = as.numeric(cb.table[1:15,3]), lower = as.numeric(cb.table[1:15,4]), 
           upper = as.numeric(cb.table[,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

dev.off()

pdf("alpha_diversity_raw_median_chao.pdf",height=5)
forestplot(labeltext = cb.table[17:31,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 13), rep(TRUE, 1)), 
           mean = as.numeric(cb.table[17:31,3]), lower = as.numeric(cb.table[17:31,4]), 
           upper = as.numeric(cb.table[,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

dev.off()
pdf("alpha_diversity_raw_median_shannon.pdf",height = 5)
forestplot(labeltext = cb.table[33:47,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 13), rep(TRUE, 1)), 
           mean = as.numeric(cb.table[33:47,3]), lower = as.numeric(cb.table[33:47,4]), 
           upper = as.numeric(cb.table[33:63,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

dev.off()

pdf("alpha_diversity_raw_median_Simpson.pdf",height = 5)
forestplot(labeltext = cb.table[49:63,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 13), rep(TRUE, 1)), 
           mean = as.numeric(cb.table[49:63,3]), lower = as.numeric(cb.table[49:63,4]), 
           upper = as.numeric(cb.table[49:63,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `15` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))
dev.off()
###remove BD_H and biopsies 
##observed_otus
id <- c(3,7:11)

med.g1 <- observed_otus$median.in.non.depression[-id]
q1.g1 <-  observed_otus$X25..IQR.in.non.depression[-id]
q3.g1 <- observed_otus$X75..IQR.in.non.depression[-id] 
n.g1 <- observed_otus$non.depression[-id]

med.g2 <- observed_otus$median.in.depression[-id]
q1.g2 <-  observed_otus$X25..IQR.in.depression[-id]
q3.g2 <- observed_otus$X75..IQR.in.depression[-id]
n.g2 <- observed_otus$depression[-id]

qe.obs <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
             q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.obs)

effect.size <- qe.obs$yi
QE.lb <- qe.obs$yi - qnorm(0.975) * sqrt(qe.obs$vi)
QE.ub <-  qe.obs$yi + qnorm(0.975) * sqrt(qe.obs$vi)


QE <- paste0(as.character(sprintf("%.2f", round(effect.size, 2))), 
             " [", as.character(sprintf("%.2f", round(QE.lb, 2))), ", ", 
             as.character(sprintf("%.2f", round(QE.ub, 2))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.2f", round(qe.obs$b, 2))), 
                         " [", 
                         as.character(sprintf("%.2f", round(qe.obs$ci.lb, 2))), 
                         ", ",  
                         as.character(sprintf("%.2f", round(qe.obs$ci.ub, 2))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Observed OTU",rownames(observed_otus[-id,]),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.obs$b)
lb.vals.medians <- c(NA, QE.lb, qe.obs$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.obs$ci.ub)

library("forestplot")                     

forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "Observed OTU", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `9` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.rm <- cbind.data.frame(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table.rem <- cb.table.rm


##chao1
med.g1 <- chao1$median.in.non.depression[-id]
q1.g1 <-  chao1$X25..IQR.in.non.depression[-id] 
q3.g1 <- chao1$X75..IQR.in.non.depression[-id] 
n.g1 <- chao1$non.depression[-id]

med.g2 <- chao1$median.in.depression[-id]
q1.g2 <-  chao1$X25..IQR.in.depression[-id]
q3.g2 <- chao1$X75..IQR.in.depression[-id]
n.g2 <- chao1$depression[-id]

qe.chao1 <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
               q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.chao1)


effect.size <- qe.chao1$yi
QE.lb <- qe.chao1$yi - qnorm(0.975) * sqrt(qe.chao1$vi)
QE.ub <-  qe.chao1$yi + qnorm(0.975) * sqrt(qe.chao1$vi)


QE <- paste0(as.character(sprintf("%.2f", round(effect.size, 2))), 
             " [", as.character(sprintf("%.2f", round(QE.lb, 2))), ", ", 
             as.character(sprintf("%.2f", round(QE.ub, 2))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.2f", round(qe.chao1$b, 2))), 
                         " [", 
                         as.character(sprintf("%.2f", round(qe.chao1$ci.lb, 2))), 
                         ", ",  
                         as.character(sprintf("%.2f", round(qe.chao1$ci.ub, 2))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Chao1 index",rownames(chao1)[-id],"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.chao1$b)
lb.vals.medians <- c(NA, QE.lb, qe.chao1$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.chao1$ci.ub)

forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "chao1", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `9` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.rm <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table.rem <- rbind.data.frame(cb.table.rem,"",cb.table.rm)


###shannon
med.g1 <- shannon$median.in.non.depression[-id]
q1.g1 <-  shannon$X25..IQR.in.non.depression[-id]
q3.g1 <- shannon$X75..IQR.in.non.depression[-id]
n.g1 <- shannon$non.depression[-id]

med.g2 <- shannon$median.in.depression[-id]
q1.g2 <-  shannon$X25..IQR.in.depression[-id]
q3.g2 <- shannon$X75..IQR.in.depression[-id]
n.g2 <- shannon$depression[-id]

qe.shannon <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                 q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.shannon)

effect.size <- qe.shannon$yi
QE.lb <- qe.shannon$yi - qnorm(0.975) * sqrt(qe.shannon$vi)
QE.ub <-  qe.shannon$yi + qnorm(0.975) * sqrt(qe.shannon$vi)


QE <- paste0(as.character(sprintf("%.4f", round(effect.size, 4))), 
             " [", as.character(sprintf("%.4f", round(QE.lb, 4))), ", ", 
             as.character(sprintf("%.4f", round(QE.ub, 4))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.4f", round(qe.shannon$b, 4))), 
                         " [", 
                         as.character(sprintf("%.4f", round(qe.shannon$ci.lb, 4))), 
                         ", ",  
                         as.character(sprintf("%.4f", round(qe.shannon$ci.ub, 4))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Shannon index",rownames(shannon[-id,]),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.shannon$b)
lb.vals.medians <- c(NA, QE.lb, qe.shannon$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.shannon$ci.ub)


forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "shannon", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `9` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

cb.table.rm <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table.rem <- rbind.data.frame(cb.table.rem,"",cb.table.rm)

###find out outline


library(dmetar)
##simpson
med.g1 <- simpson$median.in.non.depression[-id]
q1.g1 <-  simpson$X25..IQR.in.non.depression[-id] 
q3.g1 <- simpson$X75..IQR.in.non.depression[-id] 
n.g1 <- simpson$non.depression[-id]

med.g2 <- simpson$median.in.depression[-id]
q1.g2 <-  simpson$X25..IQR.in.depression[-id]
q3.g2 <- simpson$X75..IQR.in.depression[-id]
n.g2 <- simpson$depression[-id]

qe.simpson <- qe(q1.g1 = q1.g1, med.g1 = med.g1, q3.g1 = q3.g1,n.g1 = n.g1, 
                 q1.g2 = q1.g2, med.g2 = med.g2, q3.g2 = q3.g2, n.g2 = n.g2)
plot(qe.simpson)

effect.size <- qe.simpson$yi
QE.lb <- qe.simpson$yi - qnorm(0.975) * sqrt(qe.simpson$vi)
QE.ub <-  qe.simpson$yi + qnorm(0.975) * sqrt(qe.simpson$vi)


QE <- paste0(as.character(sprintf("%.4f", round(effect.size, 4))), 
             " [", as.character(sprintf("%.4f", round(QE.lb, 4))), ", ", 
             as.character(sprintf("%.4f", round(QE.ub, 4))), "]")

pooled.est.all <- paste0(as.character(sprintf("%.4f", round(qe.simpson$b, 4))), 
                         " [", 
                         as.character(sprintf("%.4f", round(qe.simpson$ci.lb, 4))), 
                         ", ",  
                         as.character(sprintf("%.4f", round(qe.simpson$ci.ub, 4))), 
                         "]")

tabletext.medians <- cbind(QE, pooled.est.all)
right.col.medians <- c("Difference of Median", QE, pooled.est.all) 
tabletext.medians <- cbind(c("Simpson index",rownames(simpson[-id,]),"RE Model"), right.col.medians)

mean.vals.medians <- c(NA, effect.size, qe.simpson$b)
lb.vals.medians <- c(NA, QE.lb, qe.simpson$ci.lb)
ub.vals.medians <- c(NA, QE.ub, qe.simpson$ci.ub)


forestplot(labeltext = tabletext.medians, graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = mean.vals.medians, lower = lb.vals.medians, 
           upper = ub.vals.medians, xlab = "simpson", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `9` = gpar(lwd = 1, col = "#000044")), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

####combind
cb.table.rm <- cbind(tabletext.medians,mean.vals.medians,lb.vals.medians,ub.vals.medians)
cb.table.rem <- rbind.data.frame(cb.table.rem,"",cb.table.rm)
id2title <- read.table("../ID2Title_meta_analysis_study.txt",sep="\t",head=T)
id2title$id <- paste("ref",id2title$id,sep= "_")
id2title$id[3] <- "ref27_BD_H"
id2title<-id2title[-6,]
cb.table.rem[c(2:8,12:18,22:28,32:38),"V1"] <- id2title$title

pdf("alpha_diversity_raw_median_obs_chao_rm65andBDH.pdf")
forestplot(labeltext = cb.table.rem[1:19,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 3),rep(FALSE, 7),rep(TRUE, 1)), 
           mean = as.numeric(cb.table.rem[1:19,3]), lower = as.numeric(cb.table.rem[1:19,4]), 
           upper = as.numeric(cb.table.rem[1:19,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `12` = gpar(lwd = 1, col = "#000044"),
                             `9` = gpar(lwd = 1, col = "#000044",lty = 2),
                             `19` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))

dev.off()

pdf("alpha_diversity_raw_median_shannon_rm65andBDH.pdf",height = 5)
forestplot(labeltext = cb.table.rem[21:29,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = as.numeric(cb.table.rem[21:29,3]), lower = as.numeric(cb.table.rem[21:29,4]), 
           upper = as.numeric(cb.table.rem[21:29,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                              `9` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))
dev.off()

pdf("alpha_diversity_raw_median_simpson_rm65andBDH.pdf",height = 5)
forestplot(labeltext = cb.table.rem[31:39,1:2], graph.pos = 2, 
           is.summary = c(TRUE, rep(FALSE, 7), rep(TRUE, 1)), 
           mean = as.numeric(cb.table.rem[31:39,3]), lower = as.numeric(cb.table.rem[31:39,4]), 
           upper = as.numeric(cb.table.rem[31:39,5]), xlab = "Difference of Median", vertices = TRUE, 
           line.margin = 0.7, new_page = FALSE, 
           align = c("l", "r"), boxsize = 0.25, 
           hrzl_lines = list(`2` = gpar(lwd = 1, col = "#000044"), 
                             `9` = gpar(lwd = 1, col = "#000044",lty = 2)), 
           txt_gp = fpTxtGp(label = gpar(cex = 1.1), ticks = gpar(cex = 1.1), 
                            xlab = gpar(cex = 1.2), title = gpar(cex = 1.2)))
dev.off()

find.outliers(qe.simpson)
find.outliers(qe.shannon)
