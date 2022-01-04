rm(list = ls())

library("igraph")
library("network")
library("sna")
library("ggraph")
library("visNetwork")
library("threejs")
library("networkD3")
library("ndtv")
library("Hmisc")
library("data.table")
library("meta")

noise.removal <- function(data, percent=0.1, top=NULL){
    matrix <- data
    rmnoise.bac <- apply(matrix, 1, function(x){sum(x>0)/dim(matrix)[2]}) > percent 
    matrix_1 <- matrix[rmnoise.bac,]
    return(matrix_1)
}

source("~/Project/Review/meta-analysis/taxname.reformat.R")
genus <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.count.genus.prof",sep="\t",check.names = F)
metadata <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.metadata.prof",sep="\t",check.names = F)

genus.denoinzed <- t(noise.removal(genus))

pid <- pmatch(rownames(metadata),rownames(genus.denoinzed))
sum(is.na(pid))
genus.denoinzed <- as.data.frame(genus.denoinzed[pid,])
ref27.treat.id <- metadata$subjectID == "Ref27" & grepl("-2",rownames(metadata)) & !grepl("H",rownames(metadata))
metadata <- metadata[!ref27.treat.id,]
genus.denoinzed <- genus.denoinzed[!ref27.treat.id,]

corr<- rcorr(as.matrix(genus.denoinzed),type="spearman")
dim(genus.denoinzed)

flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
    )
}
corr.wideformat <- flattenCorrMatrix(corr$r, corr$P)
corr.wideformat.all.sample.correlation <- corr.wideformat
###depression 
pid <- pmatch(rownames(metadata[metadata$Depression==1,]),rownames(genus.denoinzed))
genus.denoinzed.dep <- genus.denoinzed[pid,]
med.genus.denoinzed.dep <- apply(genus.denoinzed.dep, 2, median)
corr.dep<- rcorr(as.matrix(genus.denoinzed.dep),type="spearman")
corr.dep.wideformat <- flattenCorrMatrix(corr.dep$r, corr.dep$P)
colnames(corr.dep.wideformat) <- paste("dep",colnames(corr.dep.wideformat),sep = "_")
hist(corr.dep.wideformat$dep_p,breaks = 50)
###control
pid <- pmatch(rownames(metadata[metadata$Depression==0,]),rownames(genus.denoinzed))
genus.denoinzed.ctl  <- genus.denoinzed[pid,]
med.genus.denoinzed.ctl <- apply(genus.denoinzed.ctl, 2, median)
corr.ctl<- rcorr(as.matrix(genus.denoinzed.ctl),type="spearman")
corr.ctl.wideformat <- flattenCorrMatrix(corr.ctl$r, corr.ctl$P)
colnames(corr.ctl.wideformat) <- paste("ctl",colnames(corr.ctl.wideformat),sep = "_")
hist(corr.ctl.wideformat$ctl_p,breaks = 50)

all.equal(corr.wideformat$row,corr.dep.wideformat$dep_row)
all.equal(corr.wideformat$column,corr.dep.wideformat$dep_column)

all.equal(corr.ctl.wideformat$ctl_row,corr.dep.wideformat$dep_row)
all.equal(corr.ctl.wideformat$ctl_column,corr.dep.wideformat$dep_column)

cor.data <- cbind(corr.wideformat,corr.dep.wideformat,corr.ctl.wideformat)

conflict.cor <- cor.data[(!is.na(cor.data$ctl_p) & (!is.na(cor.data$dep_p)) 
                          & cor.data$ctl_p <0.05 & cor.data$dep_p <0.05) & 
                             (cor.data$ctl_cor > 0 & cor.data$dep_cor < 0),]
dim(conflict.cor)
consist.cor <- cor.data[(!is.na(cor.data$ctl_p) & (!is.na(cor.data$dep_p)) 
                         & cor.data$ctl_p <0.05 & cor.data$dep_p <0.05) & 
                            (cor.data$ctl_cor > 0 & cor.data$dep_cor >0) | 
                            (cor.data$ctl_cor < 0 & cor.data$dep_cor <0),]
dim(consist.cor)

###correlation meta-analysis
study.list <- unique(metadata$subjectID)

for(i in study.list){
    pid <- pmatch(rownames(metadata[metadata$subjectID==i,]),rownames(genus.denoinzed))
    genus.sep  <- genus.denoinzed[pid,]
    
    corr.matrix <- rcorr(as.matrix(genus.sep ),type="spearman")
    corr.wideformat <- flattenCorrMatrix(corr.matrix$r, corr.matrix$P)
    colnames(corr.wideformat)[3] <- "correlation" 
    corr.wideformat$study <- i
    txt <- paste("corr.",i," <- corr.wideformat",sep="")
    eval((parse(text =txt)))
    
    ##depression
    pid <- pmatch(rownames(metadata[metadata$subjectID==i & metadata$Depression==1,]),rownames(genus.denoinzed))
    genus.sep.dep  <- genus.denoinzed[pid,]
    
    corr.matrix <- rcorr(as.matrix(genus.sep.dep ),type="spearman")
    corr.wideformat.dep <- flattenCorrMatrix(corr.matrix$r, corr.matrix$P)
    colnames(corr.wideformat.dep)[3] <- "correlation" 
    corr.wideformat.dep$study <- i
    txt <- paste("corr.dep",i," <- corr.wideformat.dep",sep="")
    eval((parse(text =txt)))
    
    ##control
    pid <- pmatch(rownames(metadata[metadata$subjectID==i & metadata$Depression==0,]),rownames(genus.denoinzed))
    genus.sep.ctl  <- genus.denoinzed[pid,]
    
    corr.matrix <- rcorr(as.matrix(genus.sep.ctl ),type="spearman")
    corr.wideformat.ctl <- flattenCorrMatrix(corr.matrix$r, corr.matrix$P)
    colnames(corr.wideformat.ctl)[3] <- "correlation" 
    corr.wideformat.ctl$study <- i
    txt <- paste("corr.ctl",i," <- corr.wideformat.ctl",sep="")
    eval((parse(text =txt)))
    
    
}

stufy.n <-table(metadata$subjectID)

par(mfcol=c(3,3))
cor.id <- ls(pattern="corr.Ref\\d*$")
length.cor <- c()
for(i in cor.id){
    txt <- paste("cor.fild <- ",i ,sep="")
    eval((parse(text =txt)))
    hist(cor.fild[,4],breaks = 50)
    length.cor <- c(length.cor,dim(cor.fild)[1])
}


all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref23[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref27[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref40[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref60[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref65[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref77[,1:2]),check.attributes=F)
all.equal(as.data.table(corr.Ref1[,1:2]),as.data.table(corr.Ref79[,1:2]),check.attributes=F)

meta.cor.genus.fixed <- c()
meta.cor.genus.random <- c()
for(i in 1:dim(corr.Ref1)[1]){
    dat.sep.genus <- c()
    for(j in cor.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        m.cor <- metacor(cor = correlation, 
                         n = sample.size,
                         studlab = study,
                         data = dat.sep.genus,
                         title = i,sm="ZCOR",
                         backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.cor.genus.random <- rbind.data.frame(meta.cor.genus.random,
                                           cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.cor.genus.fixed <- rbind.data.frame(meta.cor.genus.fixed,
                                                 cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
                                           
    }
    else{}
    
}

###depressive
meta.cor.genus.dep.fixed <- c()
meta.cor.genus.dep.random <- c()
cor.dep.id <-  ls(pattern="corr.depRef\\d*$")
for(i in 1:dim(corr.depRef1)[1]){
    dat.sep.genus <- c()
    for(j in cor.dep.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    rm.id <- which(dat.sep.genus$correlation == 1)
    if(length(rm.id)){
        dat.sep.genus <- dat.sep.genus[-rm.id,]
    }
    
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        m.cor <- metacor(cor = correlation, 
                         n = sample.size,
                         studlab = study,
                         data = dat.sep.genus,
                         title = i,sm="ZCOR",
                         backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.cor.genus.dep.random <- rbind.data.frame(meta.cor.genus.dep.random,
                                                  cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.cor.genus.dep.fixed <- rbind.data.frame(meta.cor.genus.dep.fixed,
                                                 cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
        
    }
    else{}
    
}

###control
meta.cor.genus.ctl.fixed <- c()
meta.cor.genus.ctl.random <- c()
cor.ctl.id <-  ls(pattern="corr.ctlRef\\d*$")
for(i in 1:dim(corr.ctlRef1)[1]){
    dat.sep.genus <- c()
    for(j in cor.ctl.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    rm.id <- which(dat.sep.genus$correlation == 1)
    if(length(rm.id)){
        dat.sep.genus <- dat.sep.genus[-rm.id,]
    }
    
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        
        m.cor <- metacor(cor = correlation, 
                         n = sample.size,
                         studlab = study,
                         data = dat.sep.genus,
                         title = i,sm="ZCOR",
                         backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.cor.genus.ctl.random <- rbind.data.frame(meta.cor.genus.ctl.random,
                                                      cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.cor.genus.ctl.fixed <- rbind.data.frame(meta.cor.genus.ctl.fixed,
                                                     cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
        
    }
    else{}
    
}


dim(meta.cor.genus.random)
colnames(meta.cor.genus.random)[1:2] <- c("from","to")
meta.cor.genus.random$p.adjust <- p.adjust(meta.cor.genus.random$p,method = "BH")
genus.cor.pick <- meta.cor.genus.random[meta.cor.genus.random$p.adjust <0.05 & (abs(meta.cor.genus.random$TE) >0.4),]
dim(genus.cor.pick)
write.table(genus.cor.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.random.txt",sep="\t",quote = F)

dim(meta.cor.genus.dep.random)
colnames(meta.cor.genus.dep.random)[1:2] <- c("from","to")
meta.cor.genus.dep.random$p.adjust <- p.adjust(meta.cor.genus.dep.random$p,method = "BH")
genus.cor.dep.pick <- meta.cor.genus.dep.random[meta.cor.genus.dep.random$p.adjust <0.05 & (abs(meta.cor.genus.dep.random$TE) >0.4),]
dim(genus.cor.dep.pick)
write.table(genus.cor.dep.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.dep.random.txt",sep="\t",quote = F)

dim(meta.cor.genus.ctl.random)
colnames(meta.cor.genus.ctl.random)[1:2] <- c("from","to")
meta.cor.genus.ctl.random$p.adjust <- p.adjust(meta.cor.genus.ctl.random$p,method = "BH")
genus.cor.ctl.pick <- meta.cor.genus.ctl.random[meta.cor.genus.ctl.random$p.adjust <0.05 & meta.cor.genus.ctl.random$TE != "" & (abs(meta.cor.genus.ctl.random$TE) >0.4),]
dim(genus.cor.ctl.pick)
write.table(genus.cor.ctl.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.ctl.random.txt",sep="\t",quote = F)

hist(meta.cor.genus.dep.random$TE)
hist(meta.cor.genus.ctl.random$TE,add=T)

hist(meta.cor.genus.dep.random$p.adjust)
hist(meta.cor.genus.ctl.random$p.adjust)
###graphic data
mark.list <- read.table("corr/crude.and.adjust.model.marker.list.txt",head=T,row.names = 1,check.names = F)
mark.list$x <- gsub("k__b","d__B",mark.list$x)
mark.list$genus <-  extract.genus(mark.list$x)[,"genus"]
mark.list <- mark.list[(!grepl("s__",mark.list$x)) & grepl("g__",mark.list$x) ,]

node.meta <- as.data.frame(table(c(genus.cor.pick$from,genus.cor.pick$to)))
colnames(node.meta) <- c("id","size")
node.meta <- node.meta[order(node.meta$size,decreasing = T),]

node.meta.rename<- extract.genus(node.meta$id)
node.meta$nodelabel <- node.meta.rename[pmatch(node.meta$id , node.meta.rename$name_of_tax),"genus"]
node.meta$phlym <- extract.genus.phlym(node.meta$id)[,"phlym.name"]
node.meta$id <-  extract.genus(node.meta$id)[,"genus"]
node.meta$enrich <- "N"
node.meta$enrich[node.meta$nodelabel %in%  mark.list$genus] <-"Y"

link.meta <- genus.cor.pick
link.meta$dir <- 1
link.meta$dir[link.meta$TE <0] <- 0
link.meta$weight <- abs(link.meta$TE)
colnames(link.meta)[c(1:2)] <- c("from","to")
link.meta$from <- extract.genus(link.meta$from)[,"genus"]
link.meta$to <- extract.genus(link.meta$to)[,"genus"]

write.table(node.meta,"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.node.txt",sep="\t",quote = F,row.names = F)
write.table(link.meta[,c(1:3,17,18,19)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.link.txt",sep="\t",quote = F,row.names = F)


net.meta <- graph_from_data_frame(d=link.meta, vertices=node.meta, directed=F)
net.meta
net.meta <- simplify(net.meta, remove.multiple = F, remove.loops = T)
plot(net.meta,vertex.label.font=2)
l <- layout_in_circle(net.meta)
plot(net.meta, layout=l,edge.curved=.1,vertex.label.font=2)
l <- layout_with_fr(net.meta)
plot(net.meta, layout = l, vertex.label = "")
l <- layout_with_kk(net.meta)
plot(net.meta, layout=l, vertex.label.font=2)


##dep

node.dep.meta <- as.data.frame(table(c(genus.cor.dep.pick$from,genus.cor.dep.pick$to)))
colnames(node.dep.meta) <- c("id","size")
node.dep.meta <- node.dep.meta[order(node.dep.meta$size,decreasing = T),]

node.dep.meta.rename<- extract.genus(node.dep.meta$id)
node.dep.meta$nodelabel <- node.dep.meta.rename[pmatch(node.dep.meta$id , node.dep.meta.rename$name_of_tax),"genus"]
node.dep.meta$phlym <- extract.genus.phlym(node.dep.meta$id)[,"phlym.name"]

link.dep.meta <- genus.cor.dep.pick
link.dep.meta$dir <- 1
link.dep.meta$dir[link.dep.meta$TE <0] <- 0
link.dep.meta$weight <- abs(link.dep.meta$TE)
colnames(link.dep.meta)[c(1:2)] <- c("from","to")
link.dep.meta$from <- extract.genus(link.dep.meta$from)[,"genus"]
link.dep.meta$to <- extract.genus(link.dep.meta$to)[,"genus"]
node.dep.meta$enrich <- "NA"
node.dep.meta$enrich[node.dep.meta$nodelabel %in%  mark.list$genus] <-"Y"

write.table(node.dep.meta,"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.dep.node.txt",sep="\t",quote = F,row.names = F)


###ctl
node.ctl.meta <- as.data.frame(table(c(genus.cor.ctl.pick$from,genus.cor.ctl.pick$to)))
colnames(node.ctl.meta) <- c("id","size")
node.ctl.meta <- node.ctl.meta[order(node.ctl.meta$size,decreasing = T),]

node.ctl.meta.rename<- extract.genus(node.ctl.meta$id)
node.ctl.meta$nodelabel <- node.ctl.meta.rename[pmatch(node.ctl.meta$id , node.ctl.meta.rename$name_of_tax),"genus"]
node.ctl.meta$phlym <- extract.genus.phlym(node.ctl.meta$id)[,"phlym.name"]
node.ctl.meta$enrich <- "NA"
node.ctl.meta$enrich[node.ctl.meta$id %in%  mark.list$x] <-"Y"

link.ctl.meta <- genus.cor.ctl.pick
link.ctl.meta$dir <- 1
link.ctl.meta$dir[link.ctl.meta$TE <0] <- 0
link.ctl.meta$weight <- abs(link.ctl.meta$TE)
colnames(link.ctl.meta)[c(1:2)] <- c("from","to")
link.ctl.meta$from <- extract.genus(link.ctl.meta$from)[,"genus"]
link.ctl.meta$to <- extract.genus(link.ctl.meta$to)[,"genus"]

write.table(node.ctl.meta,"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.ctl.node.txt",sep="\t",quote = F,row.names = F)


#########remove ref65
#########remove ref65
#########remove ref65
#########remove ref65
##overall
cor.id <- cor.id[-6]
stufy.n <- stufy.n[-6]
meta.rm65.cor.genus.fixed <- c()
meta.rm65.cor.genus.random <- c()
for(i in 1:dim(corr.Ref1)[1]){
    dat.sep.genus <- c()
    for(j in cor.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        m.cor <- metacor(cor = correlation, 
                              n = sample.size,
                              studlab = study,
                              data = dat.sep.genus,
                              title = i,sm="ZCOR",
                              backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.rm65.cor.genus.random <- rbind.data.frame(meta.rm65.cor.genus.random,
                                                       cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.rm65.cor.genus.fixed <- rbind.data.frame(meta.rm65.cor.genus.fixed,
                                                      cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
        
    }
    else{}
    
}


###depressive
meta.rm65.cor.genus.dep.fixed <- c()
meta.rm65.cor.genus.dep.random <- c()
cor.dep.id <-  ls(pattern="corr.depRef\\d*$")
cor.dep.id <- cor.dep.id[-6]

for(i in 1:dim(corr.depRef1)[1]){
    dat.sep.genus <- c()
    for(j in cor.dep.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    rm.id <- which(dat.sep.genus$correlation == 1)
    if(length(rm.id)){
        dat.sep.genus <- dat.sep.genus[-rm.id,]
    }
    
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        m.cor <- metacor(cor = correlation, 
                         n = sample.size,
                         studlab = study,
                         data = dat.sep.genus,
                         title = i,sm="ZCOR",
                         backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.rm65.cor.genus.dep.random <- rbind.data.frame(meta.rm65.cor.genus.dep.random,
                                                           cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.rm65.cor.genus.dep.fixed <- rbind.data.frame(meta.rm65.cor.genus.dep.fixed,
                                                          cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
        
    }
    else{}
    
}

###control
meta.rm65.cor.genus.ctl.fixed <- c()
meta.rm65.cor.genus.ctl.random <- c()
cor.ctl.id <-  ls(pattern="corr.ctlRef\\d*$")
cor.ctl.id <-cor.ctl.id[-6]
for(i in 1:dim(corr.ctlRef1)[1]){
    dat.sep.genus <- c()
    for(j in cor.ctl.id){
        txt <- paste("cor.fild <- ",j ,sep="")
        eval((parse(text =txt)))
        dat.sep.genus <- rbind.data.frame(dat.sep.genus,cor.fild[i,])
    }
    dat.sep.genus$sample.size <- stufy.n
    rm.id <- which(dat.sep.genus$correlation == 1)
    if(length(rm.id)){
        dat.sep.genus <- dat.sep.genus[-rm.id,]
    }
    
    if(sum(is.na(dat.sep.genus$p))<4){
        dat.sep.genus <- as.data.frame(dat.sep.genus)
        dat.sep.genus$sample.size <- as.numeric(dat.sep.genus$sample.size)
        
        m.cor <- metacor(cor = correlation, 
                         n = sample.size,
                         studlab = study,
                         data = dat.sep.genus,
                         title = i,sm="ZCOR",
                         backtransf=TRUE)
        Q.p <- m.cor$pval.Q
        m.cor <- summary(m.cor)
        meta.rm65.cor.genus.ctl.random <- rbind.data.frame(meta.rm65.cor.genus.ctl.random,
                                                           cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$random),m.cor$tau2,m.cor$I2))
        meta.rm65.cor.genus.ctl.fixed <- rbind.data.frame(meta.rm65.cor.genus.ctl.fixed,
                                                          cbind.data.frame(dat.sep.genus$row[1],dat.sep.genus$column[1],as.data.frame(m.cor$fixed),m.cor$tau2,m.cor$I2))
        
    }
    else{}
    
}

################filter

dim(meta.rm65.cor.genus.random)
colnames(meta.rm65.cor.genus.random)[1:2] <- c("from","to")
meta.rm65.cor.genus.random$p.adjust <- p.adjust(meta.rm65.cor.genus.random$p,method = "BH")
genus.cor.pick <- meta.rm65.cor.genus.random[meta.rm65.cor.genus.random$p.adjust <0.05 & (abs(meta.rm65.cor.genus.random$TE) >0.4),]
dim(genus.cor.pick)
write.table(genus.cor.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.rm65.random.txt",sep="\t",quote = F)


dim(meta.rm65.cor.genus.dep.random)
colnames(meta.rm65.cor.genus.dep.random)[1:2] <- c("from","to")
meta.rm65.cor.genus.dep.random$p.adjust <- p.adjust(meta.rm65.cor.genus.dep.random$p,method = "BH")
genus.cor.dep.pick <- meta.rm65.cor.genus.dep.random[meta.rm65.cor.genus.dep.random$p.adjust <0.05 & (abs(meta.rm65.cor.genus.dep.random$TE) >0.4),]
dim(genus.cor.dep.pick)
write.table(genus.cor.dep.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.rm65.dep.random.txt",sep="\t",quote = F)

dim(meta.rm65.cor.genus.ctl.random)
colnames(meta.rm65.cor.genus.ctl.random)[1:2] <- c("from","to")
meta.rm65.cor.genus.ctl.random$p.adjust <- p.adjust(meta.rm65.cor.genus.ctl.random$p,method = "BH")
genus.cor.ctl.pick <- meta.rm65.cor.genus.ctl.random[meta.rm65.cor.genus.ctl.random$p.adjust <0.05 & meta.rm65.cor.genus.ctl.random$TE != "" & (abs(meta.rm65.cor.genus.ctl.random$TE) >0.4),]
dim(genus.cor.ctl.pick)
write.table(genus.cor.ctl.pick[,1:3],"~/Project/Review/meta-analysis/corr/genus.correlation.meta.rm65.ctl.random.txt",sep="\t",quote = F)


###graphic data
mark.list <- read.table("corr/crude.and.adjust.model.marker.list.txt",head=T,row.names = 1,check.names = F)
mark.list$x <- gsub("k__b","d__B",mark.list$x)
mark.list$genus <-  extract.genus(mark.list$x)[,"genus"]

node.meta.rm65 <- as.data.frame(table(c(genus.cor.pick$from,genus.cor.pick$to)))
colnames(node.meta.rm65) <- c("id","size")
node.meta.rm65 <- node.meta.rm65[order(node.meta.rm65$size,decreasing = T),]

node.meta.rm65.rename<- extract.genus(node.meta.rm65$id)
node.meta.rm65$nodelabel <- node.meta.rm65.rename[pmatch(node.meta.rm65$id , node.meta.rm65.rename$name_of_tax),"genus"]
node.meta.rm65$phlym <- extract.genus.phlym(node.meta.rm65$id)[,"phlym.name"]
node.meta.rm65$id <-  extract.genus(node.meta.rm65$id)[,"genus"]
node.meta.rm65$enrich <- "NA"
node.meta.rm65$enrich[node.meta.rm65$nodelabel %in%  mark.list$genus] <-"Y"
table(node.meta.rm65$enrich)
link.meta.rm65 <- genus.cor.pick
link.meta.rm65$dir <- 1
link.meta.rm65$dir[link.meta.rm65$TE <0] <- 0
link.meta.rm65$weight <- abs(link.meta.rm65$TE)
colnames(link.meta.rm65)[c(1:2)] <- c("from","to")
link.meta.rm65$from <- extract.genus(link.meta.rm65$from)[,"genus"]
link.meta.rm65$to <- extract.genus(link.meta.rm65$to)[,"genus"]

write.table(node.meta.rm65,"~/Project/Review/meta-analysis/corr/genus.cor.rm65.sample.node.txt",sep="\t",quote = F,row.names = F)
write.table(link.meta.rm65[,c(1:3,17,18,19)],"~/Project/Review/meta-analysis/corr/genus.cor.rm65.sample.link.txt",sep="\t",quote = F,row.names = F)


net.meta.rm65 <- graph_from_data_frame(d=link.meta.rm65, vertices=node.meta.rm65, directed=F)
net.meta.rm65
net.meta.rm65 <- simplify(net.meta.rm65, remove.multiple = F, remove.loops = T)
plot(net.meta.rm65,vertex.label.font=2)
l <- layout_in_circle(net.meta.rm65)
plot(net.meta.rm65, layout=l,edge.curved=.1,vertex.label.font=2)
l <- layout_with_fr(net.meta.rm65)
plot(net.meta.rm65, layout = l, vertex.label = "")
l <- layout_with_kk(net.meta.rm65)
plot(net.meta.rm65, layout=l, vertex.label.font=2)


##dep

node.dep.meta.rm65 <- as.data.frame(table(c(genus.cor.dep.pick$from,genus.cor.dep.pick$to)))
colnames(node.dep.meta.rm65) <- c("id","size")
node.dep.meta.rm65 <- node.dep.meta.rm65[order(node.dep.meta.rm65$size,decreasing = T),]

node.dep.meta.rm65.rename<- extract.genus(node.dep.meta.rm65$id)
node.dep.meta.rm65$nodelabel <- node.dep.meta.rm65.rename[pmatch(node.dep.meta.rm65$id , node.dep.meta.rm65.rename$name_of_tax),"genus"]
node.dep.meta.rm65$phlym <- extract.genus.phlym(node.dep.meta.rm65$id)[,"phlym.name"]

link.dep.meta.rm65 <- genus.cor.dep.pick
link.dep.meta.rm65$dir <- 1
link.dep.meta.rm65$dir[link.dep.meta.rm65$TE <0] <- 0
link.dep.meta.rm65$weight <- abs(link.dep.meta.rm65$TE)
colnames(link.dep.meta.rm65)[c(1:2)] <- c("from","to")
link.dep.meta.rm65$from <- extract.genus(link.dep.meta.rm65$from)[,"genus"]
link.dep.meta.rm65$to <- extract.genus(link.dep.meta.rm65$to)[,"genus"]
node.dep.meta.rm65$enrich <- "NA"
node.dep.meta.rm65$enrich[node.dep.meta.rm65$nodelabel %in%  mark.list$genus] <-"Y"
table(node.dep.meta.rm65$enrich)
write.table(node.dep.meta.rm65,"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.dep.rm65.node.txt",sep="\t",quote = F,row.names = F)
write.table(link.dep.meta.rm65[,c(1:3,17,18,19)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.dep.rm65.link.txt",sep="\t",quote = F,row.names = F)

###ctl
node.ctl.meta.rm65 <- as.data.frame(table(c(genus.cor.ctl.pick$from,genus.cor.ctl.pick$to)))
colnames(node.ctl.meta.rm65) <- c("id","size")
node.ctl.meta.rm65 <- node.ctl.meta.rm65[order(node.ctl.meta.rm65$size,decreasing = T),]

node.ctl.meta.rm65.rename<- extract.genus(node.ctl.meta.rm65$id)
node.ctl.meta.rm65$nodelabel <- node.ctl.meta.rm65.rename[pmatch(node.ctl.meta.rm65$id , node.ctl.meta.rm65.rename$name_of_tax),"genus"]
node.ctl.meta.rm65$phlym <- extract.genus.phlym(node.ctl.meta.rm65$id)[,"phlym.name"]
node.ctl.meta.rm65$enrich <- "NA"
node.ctl.meta.rm65$enrich[node.ctl.meta.rm65$id %in%  mark.list$x] <-"Y"

link.ctl.meta.rm65 <- genus.cor.ctl.pick
link.ctl.meta.rm65$dir <- 1
link.ctl.meta.rm65$dir[link.ctl.meta.rm65$TE <0] <- 0
link.ctl.meta.rm65$weight <- abs(link.ctl.meta.rm65$TE)
colnames(link.ctl.meta.rm65)[c(1:2)] <- c("from","to")
link.ctl.meta.rm65$from <- extract.genus(link.ctl.meta.rm65$from)[,"genus"]
link.ctl.meta.rm65$to <- extract.genus(link.ctl.meta.rm65$to)[,"genus"]

write.table(node.ctl.meta.rm65,"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.ctl.rm65.node.txt",sep="\t",quote = F,row.names = F)
write.table(link.ctl.meta.rm65[,c(1:3,17,18,19)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.ctl.rm65.link.txt",sep="\t",quote = F,row.names = F)

###compared all sample and remove 65 
link.meta$label <- paste(link.meta$from,link.meta$to,sep="_vs_")
link.meta.rm65$label <- paste(link.meta.rm65$from,link.meta.rm65$to,sep="_vs_")

robust.link.ctl <- intersect(link.meta$label,link.meta.rm65$label)
robust.link.ctl <- data.frame(robust.link.ctl)
robust.link.ctl$dir.all <- link.meta[pmatch(robust.link.ctl[,1],link.meta$label),"dir"]
robust.link.ctl$dir.rm <- link.meta.rm65[pmatch(robust.link.ctl[,1],link.meta.rm65$label),"dir"]

robust.link.ctl$weight.all <- link.meta[pmatch(robust.link.ctl[,1],link.meta$label),"weight"]
robust.link.ctl$weight.rm <- link.meta.rm65[pmatch(robust.link.ctl[,1],link.meta.rm65$label),"weight"]
sum(robust.link.ctl$dir.all == robust.link.ctl$dir.rm) == dim(robust.link.ctl)[1]
hist(abs(robust.link.ctl$weight.all - robust.link.ctl$weight.rm))

link.meta$robust <- "N"
link.meta$robust[link.meta$label %in% robust.link.ctl$robust.link.ctl] <- "Y"
write.table(link.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.link.txt",sep="\t",quote = F,row.names = F)
write.csv(link.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.link.csv",quote = F,row.names = F)


##control
link.ctl.meta$label <- paste(link.ctl.meta$from,link.ctl.meta$to,sep="_vs_")
link.ctl.meta.rm65$label <- paste(link.ctl.meta.rm65$from,link.ctl.meta.rm65$to,sep="_vs_")

robust.link.ctl <- intersect(link.ctl.meta$label,link.ctl.meta.rm65$label)
robust.link.ctl <- data.frame(robust.link.ctl)
robust.link.ctl$dir.all <- link.ctl.meta[pmatch(robust.link.ctl[,1],link.ctl.meta$label),"dir"]
robust.link.ctl$dir.rm <- link.ctl.meta.rm65[pmatch(robust.link.ctl[,1],link.ctl.meta.rm65$label),"dir"]

robust.link.ctl$weight.all <- link.ctl.meta[pmatch(robust.link.ctl[,1],link.ctl.meta$label),"weight"]
robust.link.ctl$weight.rm <- link.ctl.meta.rm65[pmatch(robust.link.ctl[,1],link.ctl.meta.rm65$label),"weight"]
sum(robust.link.ctl$dir.all == robust.link.ctl$dir.rm) == dim(robust.link.ctl)[1]
hist(abs(robust.link.ctl$weight.all - robust.link.ctl$weight.rm))

link.ctl.meta$robust <- "N"
link.ctl.meta$robust[link.ctl.meta$label %in% robust.link.ctl$robust.link.ctl] <- "Y"
write.table(link.ctl.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.ctl.link.txt",sep="\t",quote = F,row.names = F)
write.csv(link.ctl.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.ctl.link.csv",quote = F,row.names = F)

##########depression
link.dep.meta$label <- paste(link.dep.meta$from,link.dep.meta$to,sep="_vs_")
link.dep.meta.rm65$label <- paste(link.dep.meta.rm65$from,link.dep.meta.rm65$to,sep="_vs_")

robust.link.dep <- intersect(link.dep.meta$label,link.dep.meta.rm65$label)
robust.link.dep <- data.frame(robust.link.dep)
robust.link.dep$dir.all <- link.dep.meta[pmatch(robust.link.dep[,1],link.dep.meta$label),"dir"]
robust.link.dep$dir.rm <- link.dep.meta.rm65[pmatch(robust.link.dep[,1],link.dep.meta.rm65$label),"dir"]

robust.link.dep$weight.all <- link.dep.meta[pmatch(robust.link.dep[,1],link.dep.meta$label),"weight"]
robust.link.dep$weight.rm <- link.dep.meta.rm65[pmatch(robust.link.dep[,1],link.dep.meta.rm65$label),"weight"]
sum(robust.link.dep$dir.all == robust.link.dep$dir.rm) == dim(robust.link.dep)[1]
hist(abs(robust.link.dep$weight.all - robust.link.dep$weight.rm))

link.dep.meta$robust <- "N"
link.dep.meta$robust[link.dep.meta$label %in% robust.link.dep$robust.link.dep] <- "Y"
write.table(link.dep.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.dep.link.txt",sep="\t",quote = F,row.names = F)
write.csv(link.dep.meta[,c(1:3,17,18,19,21)],"~/Project/Review/meta-analysis/corr/genus.cor.all.sample.dep.link.csv",quote = F,row.names = F)

###compare ctl vs dep
link.dep.meta.p <- link.dep.meta[link.dep.meta$robust=="Y",]
link.ctl.meta.p <- link.ctl.meta[link.ctl.meta$robust=="Y",]

#common
common.link <- as.data.frame(intersect(link.dep.meta.p$label,link.ctl.meta.p$label))
common.link$dep <- link.dep.meta.p[pmatch(common.link[,1],link.dep.meta.p$label),"dir"]
common.link$ctl <- link.ctl.meta.p[pmatch(common.link[,1],link.ctl.meta.p$label),"dir"]

sum(common.link$dep ==common.link$ctl ) == dim(common.link)[1]

#novel
novel.dep <-link.dep.meta.p[link.dep.meta.p$label %in% common.link[,1],]
novel.ctl <-link.ctl.meta.p[link.ctl.meta.p$label %in% common.link[,1],]

###strength
d.dep <- density(abs(meta.cor.genus.dep.random$TE))
plot(d.dep,col = rgb(33, 152, 164,max=255 ))

d.ctl <- density(abs(meta.cor.genus.ctl.random$TE))
lines(d.ctl,col = rgb(255, 183, 3,max=255 ),add=T)


dat <- rbind.data.frame(cbind(abs(meta.cor.genus.dep.random[meta.cor.genus.dep.random$p.adjust < 0.05,"TE"]),"Depression"),
                  cbind(abs(meta.cor.genus.ctl.random[meta.cor.genus.ctl.random$p.adjust < 0.05,"TE"]),"Control"))

wilcox.test(as.numeric(dat$V1)~dat$V2)

pdf("corr/correlation.density.pdf")
ggplot(dat, aes(x = V1, fill = V2)) + geom_density(alpha = 0.5)
dev.off()