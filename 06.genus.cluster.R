rm(list = ls())

library(vegan)
library(Mcluster.sample)
library(DirichletMultinomial)
library(reshape2)



noise.removal <- function(data, percent=0.01, top=NULL){
  matrix <- data
  rmnoise.bac <- apply(matrix, 1, function(x){sum(x>0)/dim(matrix)[2]}) > percent 
  matrix_1 <- matrix[rmnoise.bac,]
  return(matrix_1)
}


ref1 <- read.csv("~/Project/Review/Re1/Analysis/silva_138/811b49ca-cf28-40d3-bfd9-1e21c7932c74/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref1) <- ref1$index
colnames(ref1)
ref1 <- ref1[,-c(1,265:281)]
rowSums(ref1)
rownames(ref1) <- paste("ref1.",rownames(ref1),sep="")

ref23 <- read.csv("~/Project/Review/Re23/silva138/8d7d3b22-ebb3-4b71-9dfe-46b748cbd369/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref23) <- ref23$index
colnames(ref23)
ref23 <- ref23[,-c(1,617:629)]
rowSums(ref23)
rownames(ref23) <- paste("Ref23.",rownames(ref23),sep="")

ref27 <- read.csv("~/Project/Review/Re27/Analysis/silva138/e9976a0d-74ff-4109-9d27-78cd8b1877be/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref27) <- ref27$index
colnames(ref27)
ref27 <- ref27[,-c(1,269:270)]
rowSums(ref27)
rownames(ref27) <- paste("Ref27.",rownames(ref27),sep="")

ref40 <- read.csv("~/Project/Review/Re40/silva138/8be7cf1e-0dcc-46b9-9807-996561cdff21/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref40) <- ref40$index
colnames(ref40)
ref40 <- ref40[,-c(1,245:299)]
rowSums(ref40)
rownames(ref40) <- paste("Ref40.",rownames(ref40),sep="")

ref60 <- read.csv("~/Project/Review/Re60/EGAD00001004449/silva138/db80eb79-a176-4510-a408-9356f418a56d/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref60) <- ref60$index
colnames(ref60)
ref60 <- ref60[,-c(1,294:315)]
rowSums(ref60)
rownames(ref60) <- paste("Ref60.",rownames(ref60),sep="")

ref65 <- read.csv("~/Project/Review/Re65/IBD_RAR_DATA/silva138/0c84a17e-1f23-4ff4-bcad-b77c5a148b71/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref65) <- ref65$index
colnames(ref65)
ref65 <- ref65[,-c(1,474:488)]
rowSums(ref65)
rownames(ref65) <- paste("Ref65.",rownames(ref65),sep="")

ref77 <- read.csv("~/Project/Review/Re77/study_12382_raw_100520-200346/trim/flash/silva138/05e4d5b1-beed-4c42-997c-31cabb99faab/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref77) <- ref77$index
colnames(ref77)
ref77 <- ref77[,-c(1,300:303)]
rowSums(ref77)
rownames(ref77) <- paste("Ref77.",rownames(ref77),sep="")

ref79 <- read.csv("~/Project/Review/Re79/silva138/7de3da12-8452-4b5f-bc1a-171e48297eaa/data/level-6.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(ref79) <- ref79$index
colnames(ref79)
ref79 <- ref79[,-c(1,174:175)]
rowSums(ref79)
rownames(ref79) <- paste("Ref79.",rownames(ref79),sep="")

file.id <- ls(pattern="^ref\\d*$")

dat.list <- list(ref1 = ref1,
                 ref23 = ref23,
                 ref27 = ref27,
                 ref40 = ref40,
                 ref60 = ref60,
                 ref65 = ref65,
                 ref77 = ref77,
                 ref79 = ref79)

rowSums(sapply(dat.list,dim))

dat.cout.cb <- c()
for(i in  file.id){
  txt <- paste("dat <- ",i,sep="")
  eval((parse(text =txt)))
  dat<-data.frame(dat)
  dat.cout.cb <- dplyr::bind_rows(dat.cout.cb, dat)    
  
}

dat.cout.cb <- t(dat.cout.cb)
na.num <- apply(dat.cout.cb,1,function(x){sum(is.na(x))})
sort(na.num)
dim(dat.cout.cb)
table(dat.cout.cb)
dat.cout.cb[is.na(dat.cout.cb)]=0

#write.table(dat.cout.cb,"~/Project/Review/meta-analysis/allstudy.1842.sample.count.genus.prof",sep="\t",quote = F)

dat.cout.cb.denoinzed <- noise.removal(dat.cout.cb)

set.seed(0)
fit <- lapply(1:3, dmn, count = t(dat.cout.cb.denoinzed), verbose = TRUE)

lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace

plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

best <- fit[[which.min(unlist(lplc))]]
cluster.sample <- apply(mixture(best), 1, which.max)
table(cluster.sample)

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Genus", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by cluster.sampleignment strength
    arrange(value) %>%
    mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = Genus, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

ref1.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref1.meta.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref1.meta),rownames(ref1)) != c(1:dim(ref1)[1]))

ref23.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref23.meta.dada.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref23.meta),rownames(ref23)) != c(1:dim(ref23)[1]))

ref27.meta <- as.data.frame(rownames(ref27))
rownames(ref27.meta) <- rownames(ref27)
ref27.meta$Depression <- 1
ref27.meta$Depression[grepl("H",ref27.meta[,1]) | 
                        grepl("\\-2",ref27.meta[,1])] <- 0
ref27.meta <- ref27.meta[,-1,drop=F]
ref27.meta$subjectID <- "Ref27"

ref40.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref40.meta.dada2.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref40.meta),rownames(ref40)) != c(1:dim(ref40)[1]))

ref60.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref60.meta.dada2.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref60.meta),rownames(ref60)) != c(1:dim(ref60)[1]))

ref65.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref65.meta.dada2.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref65.meta),rownames(ref65)) != c(1:dim(ref65)[1]))
sample.id <- intersect(rownames(ref65),rownames(ref65.meta))
ref65.meta <- ref65.meta[sample.id,]
ref65 <- ref65[sample.id,]
sum(pmatch(rownames(ref65.meta),rownames(ref65)) != c(1:dim(ref65)[1]))

ref77.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref77.meta.dada2.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref77.meta),rownames(ref77)) != c(1:dim(ref77)[1]))

ref79.meta <- read.table("~/Project/Review/meta-analysis/metadata/Ref79.meta.dada2.MMUPHin.txt",sep="\t",head=T)
sum(pmatch(rownames(ref79.meta),rownames(ref79)) != c(1:dim(ref79)[1]))

dat.meta.cb <- c()
file.meta <- ls(pattern="ref\\d*.meta")
for(i in  file.meta){
  txt <- paste("dat <- ",i,sep="")
  eval((parse(text =txt)))
  dat<-data.frame(dat)
  dat.meta.cb <- dplyr::bind_rows(dat.meta.cb, dat)    
  
}

dim(dat.cout.cb)
dim(dat.meta.cb)

pid <- pmatch(rownames(dat.meta.cb),names(cluster.sample))

dat.meta.cb$DMM <- cluster.sample[pid]
#write.table(dat.meta.cb,"~/Project/Review/meta-analysis/allstudy.1842.sample.metadata.addclust.prof",sep="\t",quote = F)

table(dat.meta.cb$DMM,dat.meta.cb$subjectID)
table(dat.meta.cb$DMM,dat.meta.cb$Depression)
fisher.test(t(table(dat.meta.cb$DMM,dat.meta.cb$Depression)),simulate.p.value=T)

dat.meta.cb$DMM <- paste("Cluster",dat.meta.cb$DMM,sep="")
dat.meta.cb$Sex <- gsub("^male","Male",dat.meta.cb$Sex)
dat.meta.cb$Sex <- gsub("^female","Female",dat.meta.cb$Sex)

crude.m <- glm(Depression~as.factor(DMM),data = dat.meta.cb,family = binomial)
adjust.m <- glm(Depression~as.factor(DMM)+Sex+Age+BMI,data = dat.meta.cb,family = binomial)

summary(crude.m)
summary(adjust.m)
#####remove ref65
dat.cout.cb <- c()
for(i in  file.id[-6]){
  txt <- paste("dat <- ",i,sep="")
  eval((parse(text =txt)))
  dat<-data.frame(dat)
  dat.cout.cb <- dplyr::bind_rows(dat.cout.cb, dat)    
  
}

dat.cout.cb <- t(dat.cout.cb)
na.num <- apply(dat.cout.cb,1,function(x){sum(is.na(x))})
sort(na.num)
dim(dat.cout.cb)
table(dat.cout.cb)
dat.cout.cb[is.na(dat.cout.cb)]=0
#write.table(dat.cout.cb,"~/Project/Review/meta-analysis/allstudy.1842.sample.count.genus.rm65.prof",sep="\t",quote = F)

dat.cout.cb.denoinzed <- noise.removal(dat.cout.cb)

set.seed(0)
fit <- lapply(1:10, dmn, count = t(dat.cout.cb.denoinzed), verbose = TRUE)


lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace

#pdf("DMM_10_Dirichlet_Component.pdf")
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
dev.off()

best <- fit[[which.min(unlist(lplc))]]

set.seed(0)
fit <- lapply(1:5, dmn, count = t(dat.cout.cb.denoinzed), verbose = TRUE)


lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace


plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

best <- fit[[which.min(unlist(lplc))]]

cluster.sample <- apply(mixture(best), 1, which.max)
table(cluster.sample)

# for (k in seq(ncol(fitted(best)))) {
#   d <- melt(fitted(best))
#   colnames(d) <- c("Genus", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange OTUs by cluster.sampleignment strength
#     arrange(value) %>%
#     mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
#     # Only show the most important drivers
#     #filter(abs(value) > quantile(abs(value), 0.95))     
#     top_n(10)
#   txt <- paste("d",k,"<-d",sep="")
#   eval((parse(text =txt)))
#   
#   p <- ggplot(d, aes(x = Genus, y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     labs(title = paste("Top drivers: community type", k)) + theme_classic()
#   print(p)
#   
#   txt <- paste("p",k,"<-p",sep="")
#   eval((parse(text =txt)))
#   
# }
# 
# #pdf("Top10.Contribution_of_tax_to_DMM.pdf",height = 40,width = 20)
# ggarrange(p1,p2,p3,p4,p5,nrow = 5,align="h")
# dev.off()

###top taxo
group.num <- max(as.numeric(cluster.sample))
fitted1 <- fitted(fit[[1]], scale = TRUE)
fitted0 <- fitted(best, scale = TRUE)
score <- rowSums(abs(fitted0 - as.vector(fitted1)))
score <- score[rev(order(score))]
cscore <- cumsum(score)/sum(score)
score.out <- cbind(score = score, cscore)[1:10, ]
#write.table(score.out, "Top10.Contribution_of_tax_to_DMM.txt", quote = F, sep = "\t")

##top10 plot
dat <- t(dat.cout.cb)/colSums(dat.cout.cb)
name <- names(score[1:10])
dat.pick <- dat[, pmatch(name, colnames(dat))]
coln.tax<- strsplit(colnames(dat.pick) ,".f__")
colnames(dat.pick) <- as.matrix(sapply(coln.tax, function(x){x[2]}))
colnames(dat.pick) <- gsub(".g__","|",colnames(dat.pick))
colnames(dat.pick)[5] <- gsub(".__","|Unclassified",colnames(dat.pick)[5]) 

my_comparisons <- list( c("1", "3"), c("2", "3"), c("3", "4"), c("3","5") )
col <- c("#4D9DE0","#E15554","#3BB273","#E1BC29","#7768AE")
for (i in 1:10) {
  dat.tmp <- cbind.data.frame(dat.pick[, i],cluster.sample)
  colnames(dat.tmp) <- c("taxo","cluster")
  p <- ggboxplot(dat.tmp, x = "cluster", y = "taxo",
            fill = "cluster", palette = col) +
    #stat_compare_means(method = "kruskal.test", label.y = max(dat.tmp$taxo)+0.03, label.x = 2)+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = "3") +
    ylab(colnames(dat.pick)[i])+xlab("Cluster")
  txt <- paste("p",i," <- p",sep="")
  eval((parse(text =txt)))
}

#pdf("top.10.boxplot.pdf", height = 7, width = 15)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol = 5,nrow = 2,common.legend = T)
dev.off()

###median and mean of each tax
dat.pick.median <- aggregate(dat.pick,by = list(cluster.sample),median)
min.or.max <- apply(dat.pick.median, 2,function(x){sum(which.min(x) == 3 |which.max(x) == 3 ) })

# NMDS plot
dat.mds <- metaMDS(dat, distance = "bray")
dat.mds$stress
sum(pmatch(rownames(dat.mds$points),names(cluster.sample)) != c(1:1362))
#pdf("DMM_NMDS.pdf")
plot(dat.mds$points, col = col[cluster.sample], pch = 19)
legend("topright", legend=paste("Cluster",1:group.num,sep="_"),
       col=col, pch=16, cex=0.8)
dev.off() 

par(mfrow=c(2,3))
cluster.1 <- dat.mds$points[cluster.sample==1,]
plot(cluster.1, col = col[1], pch = 19,xlim=c(-1,2),ylim=c(-1,2))
cluster.2 <- dat.mds$points[cluster.sample==2,]
plot(cluster.2, col = col[2], pch = 19,xlim=c(-1,2),ylim=c(-1,2))
cluster.3 <- dat.mds$points[cluster.sample==3,]
plot(cluster.3, col = col[3], pch = 19,xlim=c(-1,2),ylim=c(-1,2))
cluster.4 <- dat.mds$points[cluster.sample==4,]
plot(cluster.4, col = col[4], pch = 19,xlim=c(-1,2),ylim=c(-1,2))
cluster.5 <- dat.mds$points[cluster.sample==5,]
plot(cluster.5, col = col[5], pch = 19,xlim=c(-1,2),ylim=c(-1,2))

dat.meta.cb <- c()
file.meta <- ls(pattern="ref\\d*.meta")
for(i in  file.meta[-6]){
  txt <- paste("dat <- ",i,sep="")
  eval((parse(text =txt)))
  dat<-data.frame(dat)
  dat.meta.cb <- dplyr::bind_rows(dat.meta.cb, dat)    
  
}

dim(dat.cout.cb)
dim(dat.meta.cb)

pid <- pmatch(rownames(dat.meta.cb),names(cluster.sample))
dat.meta.cb$DMM <- cluster.sample[pid]

#write.table(dat.meta.cb,"~/Project/Review/meta-analysis/allstudy.1842.sample.metadata.rm65.prof",sep="\t",quote = F)

###PCA
# library(ade4)
# library(RColorBrewer)
# pca<- dudi.pca(df = dat, center = T, scale = FALSE, scannf = FALSE, nf = 5)
# 
# 
# bb=pca$c1[order(sqrt(pca$c1[,1]^2+pca$c1[,2]^2),decreasing=T),][1:3,]
# pca_eig <- (pca$eig)[1:2] / sum(pca$eig)
# sample_site <- data.frame({pca$li})[1:2]
# sample_site$names <- rownames(sample_site)
# names(sample_site)[1:2] <- c('PCA1', 'PCA2')
# sum(pmatch(rownames(sample_site$names),names(cluster.sample)) != c(1:1362))
# sample_site$cluster<-as.factor(cluster.sample)
# pca_plot <- ggplot(sample_site, aes(PCA1, PCA2,color=cluster)) +
#   theme_classic()+#去掉背景框
#   geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
#   geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
#   geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
#   scale_color_manual(values = brewer.pal(6,"Set2")) + #可在这里修改点的颜色
#   theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1),
#         panel.background = element_rect(color = 'black', fill = 'transparent'),
#         legend.title=element_blank()
#   )+
#   labs(x = paste('PCA1: ', round(100 * pca_eig[1], 2), '%'), y = paste('PCA2: ', round(100 * pca_eig[2], 2), '%'))
# 
# pca_plot
# 
# ##PCOA
# dat.dist<-vegdist(dat,method='euclidean')#基于euclidean距离
# pcoa<- dudi.pco(dat.dist, scan = FALSE,nf=3)
# pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
# 
# #提取样本点坐标（前两轴）
# sample_site <- data.frame({pcoa$li})[1:2]
# sample_site$names <- rownames(sample_site)
# names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
# 
# #以最终成绩作为分组
# sum(pmatch(rownames(sample_site$names),names(cluster.sample)) != c(1:1362))
# sample_site$cluster<-as.factor(cluster.sample)
# 
# library(ggplot2)
# 
# pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2,color=cluster)) +
#   theme_classic()+#去掉背景框
#   geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
#   geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
#   geom_point(size = 1.5)+  #可在这里修改点的透明度、大小
#   scale_color_manual(values = brewer.pal(6,"Set2")) + #可在这里修改点的颜色
#   theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
#         panel.background = element_rect(color = 'black', fill = 'transparent'), 
#         legend.title=element_blank()
#   )+
#   labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 
# 
# pcoa_plot

###glm
table(dat.meta.cb$DMM,dat.meta.cb$subjectID)
table(dat.meta.cb$DMM,dat.meta.cb$Depression)

fisher.test(t(table(dat.meta.cb$DMM,dat.meta.cb$Depression)),simulate.p.value = T)

dat.meta.cb$DMM <- paste("Cluster",dat.meta.cb$DMM,sep="")
dat.meta.cb$DMM <- factor(dat.meta.cb$DMM,levels=c("Cluster3","Cluster1","Cluster2","Cluster4","Cluster5"))
dat.meta.cb$Sex <- gsub("^male","Male",dat.meta.cb$Sex)
dat.meta.cb$Sex <- gsub("^female","Female",dat.meta.cb$Sex)

crude.m <- glm(Depression~DMM,data = dat.meta.cb,family = binomial)
adjust.m <- glm(Depression~DMM+Sex+Age+BMI,data = dat.meta.cb,family = binomial)

crude.s <- summary(crude.m)
adjust.s <- summary(adjust.m)
crude.odd <- round(exp(cbind(coef(crude.m), confint(crude.m))),2)
adjust.odd <-round(exp(cbind(coef(adjust.m), confint(adjust.m))),2)
crude.odd <- cbind(crude.odd,crude.s$coefficients[,"Pr(>|z|)"])
adjust.odd <- cbind(adjust.odd,adjust.s$coefficients[,"Pr(>|z|)"])
#write.table(crude.odd,"DMM_odd.crude.txt",sep="\t",quote = F)
#write.table(adjust.odd,"DMM_odd.adjust.txt",sep="\t",quote = F)

library(dummies)
dumm.DMM <- dummy(dat.meta.cb$DMM, sep = ".")
rownames(dumm.DMM) <- rownames(dat.meta.cb)
dumm.DMM <- cbind.data.frame(dumm.DMM,dat.meta.cb[,c("Depression","Sex","Age","BMI")])
crude.dumm.m <- glm(Depression~DMM.Cluster3,data = dumm.DMM,family = binomial)
summary(crude.dumm.m)
adjust.dumm.m <- glm(Depression~DMM.Cluster3+Sex+Age+BMI,data = dumm.DMM,family = binomial)
summary(adjust.dumm.m)

regress.out <- list(crude.modle = crude.s$coefficients, adjust.model = adjust.s$coefficients)
#save(z = regress.out, file = "DMM_depression.regression.RData",row.names = T)
write.table(crude.dumm.m)