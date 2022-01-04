###########function

rm(list = ls())
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/")

pathway.compare.nocomvar <- function (pathtab, mapfile, sampleid = "sampleid", pathsum = "rel", 
                                      stat.med = "gamlss", transform = "none", comvar, 
                                      adjustvar, personid = "personid", longitudinal = "yes", 
                                      p.adjust.method = "fdr", percent.filter = 0.05, relabund.filter = 5e-05, 
                                      pooldata = FALSE) 
{
    pathlev <- paste("l", 1:length(pathtab), sep = "")
    estilist <- list()
    for (j in 1:length(pathlev)) {
        print(j)
        samlist <- rownames(pathtab[[j]])
        pathtab[[j]] <- as.data.frame(lapply(pathtab[[j]], as.character))
        pathtab[[j]] <- as.data.frame(lapply(pathtab[[j]], as.numeric))
        rownames(pathtab[[j]]) <- samlist
        pathlist <- colnames(pathtab[[j]])
        pathtest <- apply(pathtab[[j]], 2, function(x) {
            length(x[!is.na(x) & x > 0])
        })
        pathget <- pathtest[pathtest >= percent.filter * (nrow(pathtab[[j]]))]
        pathtests <- apply(pathtab[[j]], 2, function(x) {
            sum(x, na.rm = T)/sum(pathtab[[j]])
        })
        pathgets <- pathtests[pathtests > relabund.filter]
        pathname <- names(pathget)[names(pathget) %in% names(pathgets)]
        mapfile[, sampleid] <- tolower(mapfile[, sampleid])
        if (pathsum == "rel") {
            pathrel <- as.data.frame(t(apply(pathtab[[j]], 1, 
                                             function(x) x/sum(x)))[, pathname])
            pathrel[, sampleid] <- tolower(rownames(pathrel))
            pathdat <- merge(mapfile, pathrel, by = sampleid)
            if (stat.med == "gamlss" & transform != "none") {
                stop("gamlss with beta zero-inflated family should only be used for relative abundance without transformation")
            }
            if (stat.med == "lm" & transform == "asin.sqrt") {
                asintransform <- function(p) {
                    asin(sqrt(p))
                }
                pathdat[, pathname] <- apply(pathdat[, pathname], 
                                             2, asintransform)
            }
            if (stat.med == "lm" & transform == "logit") {
                logittransform <- function(p) {
                    log(p/(1 - p))
                }
                pathdat[, pathname] <- apply(pathdat[, pathname], 
                                             2, logittransform)
            }
        }
        if (pathsum == "log") {
            pathlog <- log2(pathtab[[j]][, pathname] + 1)
            pathlog[, sampleid] <- tolower(rownames(pathlog))
            pathdat <- merge(mapfile, pathlog, by = sampleid)
        }
        pathdat[, comvar] <- gdata::drop.levels(pathdat[, comvar], 
                                                reorder = FALSE)
        if (longitudinal == "yes") {
            pathdat$personid <- as.factor(pathdat[, personid])
        }
        estisum <- NULL
        for (i in 1:length(pathname)) {
            if (stat.med == "lm" & (pathsum == "rel" | 
                                    pathsum == "log")) {
                if (longitudinal == "yes") {
                    fitsum <- try(summary(lme4::lmer(stats::as.formula(paste(pathname[i], 
                                                                             paste(c(comvar, adjustvar, "(1|personid)"), 
                                                                                   collapse = "+"), sep = "~")), 
                                                     data = pathdat)))
                }
                if (longitudinal == "no") {
                    fitsum <- try(summary(stats::glm(stats::as.formula(paste(pathname[i], 
                                                                             paste(c(comvar, adjustvar), collapse = "+"), 
                                                                             sep = "~")), data = pathdat, family = "gaussian")))
                }
                if (class(fitsum) == "try-error") {
                    warning("Error in model fit, NA introduced.\n")
                    fitcoefw <- NULL
                    estisum <- plyr::rbind.fill(estisum, fitcoefw)
                }
                if (class(fitsum) != "try-error") {
                    if (length(which(rownames(fitsum$coefficients) != 
                                     "(Intercept)")) > 1) {
                        fitcoef <- as.data.frame(fitsum$coefficients[rownames(fitsum$coefficients) != 
                                                                         "(Intercept)", ])
                        if (longitudinal == "yes") {
                            fitcoef[, "Pr(>|t|)"] <- 1.96 * stats::pnorm(-abs(fitcoef[, 
                                                                                      "Estimate"]/fitcoef[, "Std. Error"]))
                        }
                        fitcoef[, "varname"] <- rownames(fitcoef)
                        fitcoef[, "id"] <- pathname[i]
                        fitcoefw <- stats::reshape(fitcoef, idvar = "id", 
                                                   timevar = "varname", direction = "wide")
                    }
                    if (length(which(rownames(fitsum$coefficients) != 
                                     "(Intercept)")) == 1) {
                        fitcoef <- as.data.frame(matrix(fitsum$coefficients[rownames(fitsum$coefficients) != 
                                                                                "(Intercept)", ], ncol = ncol(fitsum$coefficients)))
                        rownames(fitcoef) <- rownames(fitsum$coefficients)[rownames(fitsum$coefficients) != 
                                                                               "(Intercept)"]
                        colnames(fitcoef) <- colnames(fitsum$coefficients)
                        if (longitudinal == "yes") {
                            fitcoef[, "Pr(>|t|)"] <- 1.96 * stats::pnorm(-abs(fitcoef[, 
                                                                                      "Estimate"]/fitcoef[, "Std. Error"]))
                        }
                        fitcoef[, "varname"] <- rownames(fitcoef)
                        fitcoef[, "id"] <- pathname[i]
                        fitcoefw <- stats::reshape(fitcoef, idvar = "id", 
                                                   timevar = "varname", direction = "wide")
                    }
                    if (length(which(rownames(fitsum$coefficients) != 
                                     "(Intercept)")) == 0) {
                        fitcoefw <- NULL
                    }
                    estisum <- plyr::rbind.fill(estisum, fitcoefw)
                }
            }
            if (stat.med == "gamlss" & (pathsum == "log")) {
                stop("gamlss with beta zero-inflated family should only be used for relative abundance data")
            }
            if (stat.med == "gamlss" & (pathsum == "rel")) {
                if (longitudinal == "yes") {
                    testdat <- pathdat[, c(pathname[i], comvar, 
                                           adjustvar, "personid")]
                    testdat[, pathname[i]][testdat[, pathname[i]] == 
                                               1] <- 0.9999
                    testdat <- stats::na.omit(testdat)
                    fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i], 
                                                                                 paste(c(comvar, adjustvar, "random(personid)"), 
                                                                                       collapse = "+"), sep = "~")), 
                                                         family = BEZI, data = testdat, trace = FALSE), 
                                          save = TRUE))
                }
                if (longitudinal == "no") {
                    testdat <- pathdat[, c(pathname[i], comvar, 
                                           adjustvar)]
                    testdat[, pathname[i]][testdat[, pathname[i]] == 
                                               1] <- 0.9999
                    testdat <- stats::na.omit(testdat)
                    fitsum <- try(summary(gamlss::gamlss(stats::as.formula(paste(pathname[i], 
                                                                                 paste(c(comvar, adjustvar), collapse = "+"), 
                                                                                 sep = "~")), family = BEZI, data = testdat, 
                                                         trace = FALSE), save = TRUE))
                }
                if (class(fitsum) == "try-error") {
                    warning("Error in model fit, NA introduced.\n")
                    fitcoefw <- NULL
                    estisum <- plyr::rbind.fill(estisum, fitcoefw)
                }
                if (class(fitsum) != "try-error") {
                    if (length(which(rownames(fitsum$coef.table) != 
                                     "(Intercept)")) > 1) {
                        fitcoef <- as.data.frame(fitsum$coef.table[rownames(fitsum$coef.table) != 
                                                                       "(Intercept)", ])
                        fitcoef[, "varname"] <- rownames(fitcoef)
                        fitcoef[, "id"] <- pathname[i]
                        fitcoefw <- stats::reshape(fitcoef, idvar = "id", 
                                                   timevar = "varname", direction = "wide")
                    }
                    if (length(which(rownames(fitsum$coef.table) != 
                                     "(Intercept)")) == 1) {
                        fitcoef <- as.data.frame(matrix(fitsum$coef.table[rownames(fitsum$coef.table) != 
                                                                              "(Intercept)", ], ncol = ncol(fitsum$coef.table)))
                        rownames(fitcoef) <- rownames(fitsum$coef.table)[rownames(fitsum$coef.table) != 
                                                                             "(Intercept)"]
                        colnames(fitcoef) <- colnames(fitsum$coef.table)
                        fitcoef[, "varname"] <- rownames(fitcoef)
                        fitcoef[, "id"] <- pathname[i]
                        fitcoefw <- stats::reshape(fitcoef, idvar = "id", 
                                                   timevar = "varname", direction = "wide")
                    }
                    if (length(which(rownames(fitsum$coef.table) != 
                                     "(Intercept)")) == 0) {
                        fitcoefw <- NULL
                    }
                    estisum <- plyr::rbind.fill(estisum, fitcoefw)
                }
            }
        }
        if(length(adjustvar)){
            estisum[, sub(".*\\.", "pval.adjust.", colnames(estisum)[grep("Pr(>|t|)", 
                                                                          colnames(estisum))])] <- apply(estisum[, colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]], 2, stats::p.adjust, method = p.adjust.method)
        }
        else{
            estisum[, sub(".*\\.", "pval.adjust.", colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))])]<- stats::p.adjust(estisum[, colnames(estisum)[grep("Pr(>|t|)",colnames(estisum))]], method = p.adjust.method)
        }
        
        
        estilist[[j]] <- estisum
    }
    names(estilist) <- pathlev
    return(estilist)
}


#pathway.info <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/picrust2/metacyc_pathways_info.txt",sep="\t")
####reference 1
ref1.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
ref1.pathway[1:4,1:4]
colSums(ref1.pathway)
ref1.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref1/Ref1.meta.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref1.meta),paste("ref1.",colnames(ref1.pathway),sep="")) != c(1:dim(ref1.pathway)[2]))


####
ref1.pathway.meta <- ref1.meta[,c("Sex","Age","Depression")]
ref1.pathway.meta$Depression <- gsub(1,"1_Depression",ref1.pathway.meta$Depression )
ref1.pathway.meta$Depression <- gsub(0,"0_Depression",ref1.pathway.meta$Depression )

colnames(ref1.pathway.meta)[1:3] <- c("Gender","Age","Depression")
ref1.pathway.meta$sample.id <- rownames(ref1.pathway.meta)

library(metamicrobiomeR) 


colnames(ref1.pathway) <- paste("Ref1.",colnames(ref1.pathway),sep="")
ref1.pathway <- as.data.frame(t(ref1.pathway))
path1.crude<-pathway.compare.nocomvar(pathtab=list(ref1.pathway),sampleid = "sample.id",
                              mapfile=ref1.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                              comvar="Depression",adjustvar=NULL, longitudinal="no",
                              p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path1.adjust<-pathway.compare(pathtab=list(ref1.pathway),sampleid = "sample.id",
                       mapfile=ref1.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                       comvar="Depression",adjustvar=c("Gender","Age"), longitudinal="no",
                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path1.crude,"metaR.ref1.pathway.crude.txt",quote = F,sep="\t")
write.table(path1.adjust,"metaR.ref1.pathway.adjust.txt",quote = F,sep="\t")

######reference 23
ref23.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref23/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref23.pathway)
ref23.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref23/Ref23.meta.dada.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref23.meta),paste("Ref23.",colnames(ref23.pathway),sep="")) != c(1:dim(ref23.pathway)[2]))

####

ref23.pathway.meta <- ref23.meta[,c("Sex","Age","BMI","Depression")]

ref23.pathway.meta$Depression <- gsub(1,"1_Depression",ref23.pathway.meta$Depression )
ref23.pathway.meta$Depression <- gsub(0,"0_Depression",ref23.pathway.meta$Depression )

colnames(ref23.pathway.meta)[1:4] <- c("Gender","Age","BMI","Depression")

colnames(ref23.pathway) <- paste("ref23.",colnames(ref23.pathway),sep="")
ref23.pathway <- as.data.frame(t(ref23.pathway))
ref23.pathway.meta$sample.id <- rownames(ref23.pathway.meta)
path23.crude<-pathway.compare.nocomvar(pathtab=list(ref23.pathway),sampleid = "sample.id",
                                      mapfile=ref23.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                      comvar="Depression",adjustvar=NULL, longitudinal="no",
                                      p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path23.adjust<-pathway.compare(pathtab=list(ref23.pathway),sampleid = "sample.id",
                              mapfile=ref23.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                              comvar="Depression",adjustvar=c("Gender","Age","BMI"), longitudinal="no",
                              p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path23.crude,"metaR.ref23.pathway.crude.txt",quote = F,sep="\t")
write.table(path23.adjust,"metaR.ref23.pathway.adjust.txt",quote = F,sep="\t")


#######reference 27
ref27.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref27/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref27.pathway)
kingdom <- read.csv("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref27/dada2/silva138/level-1.csv",header = T,check.names = F,colClasses=c("index"="character"))
rownames(kingdom) <- kingdom$index
ref27.meta <- kingdom[,c(1,5:6)]
sum(pmatch(rownames(ref27.meta),colnames(ref27.pathway)) != c(1:dim(ref27.pathway)[2]))
ref27.pathway.meta <- ref27.meta[!(ref27.meta$Group == "BD" & ref27.meta$Time ==2),]
ref27.pathway <- ref27.pathway[,rownames(ref27.pathway.meta)]

colnames(ref27.pathway.meta)[1] <- "sample.id"


ref27.pathway.meta$Depression <- "1_Depression"
ref27.pathway.meta$Depression[ref27.pathway.meta$Group == "H"] <- "0_Depression"

ref27.pathway <- as.data.frame(t(ref27.pathway))
ref27.pathway.meta$sample.id <- rownames(ref27.pathway.meta)
path27.crude<-pathway.compare.nocomvar(pathtab=list(ref27.pathway),sampleid = "sample.id",
                                       mapfile=ref27.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path27.adjust<-path27.crude

write.table(path27.crude,"metaR.ref27.pathway.crude.txt",quote = F,sep="\t")
write.table(path27.adjust,"metaR.ref27.pathway.adjust.txt",quote = F,sep="\t")


#####reference 40
ref40.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref40/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref40.pathway)
ref40.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref40/Ref40.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref40.meta),paste("Ref40.",colnames(ref40.pathway),sep="")) != c(1:dim(ref40.pathway)[2]))

####
ref40.pathway.meta <- ref40.meta[,c("Age","BMI","Depression")]

ref40.pathway.meta$Depression <- gsub(1,"1_Depression",ref40.pathway.meta$Depression )
ref40.pathway.meta$Depression <- gsub(0,"0_Depression",ref40.pathway.meta$Depression )
colnames(ref40.pathway.meta)[1:3] <- c("Age","BMI","Depression")


colnames(ref40.pathway) <- paste("ref40.",colnames(ref40.pathway),sep="")
ref40.pathway <- as.data.frame(t(ref40.pathway))
ref40.pathway.meta$sample.id <- rownames(ref40.pathway.meta)
path40.crude<-pathway.compare.nocomvar(pathtab=list(ref40.pathway),sampleid = "sample.id",
                                       mapfile=ref40.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path40.adjust<-pathway.compare(pathtab=list(ref40.pathway),sampleid = "sample.id",
                               mapfile=ref40.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                               comvar="Depression",adjustvar=c("Age","BMI"), longitudinal="no",
                               p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path40.crude,"metaR.ref40.pathway.crude.txt",quote = F,sep="\t")
write.table(path40.adjust,"metaR.ref40.pathway.adjust.txt",quote = F,sep="\t")

#####reference 60
ref60.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref60/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref60.pathway)
ref60.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref60/Ref60.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref60.meta),paste("Ref60.",colnames(ref60.pathway),sep="")) != c(1:dim(ref60.pathway)[2]))
pid <- pmatch(paste("Ref60.",colnames(ref60.pathway),sep=""),rownames(ref60.meta))
ref60.meta <- ref60.meta[pid,]

ref60.pathway.meta <- ref60.meta[,c("Sex","Age","BMI","Depression")]
ref60.pathway.meta$Depression <- gsub(1,"1_Depression",ref60.pathway.meta$Depression )
ref60.pathway.meta$Depression <- gsub(0,"0_Depression",ref60.pathway.meta$Depression )

colnames(ref60.pathway.meta)[1:4] <- c("Gender","Age","BMI","Depression")

colnames(ref60.pathway) <- paste("ref60.",colnames(ref60.pathway),sep="")
ref60.pathway <- as.data.frame(t(ref60.pathway))
ref60.pathway.meta$sample.id <- rownames(ref60.pathway.meta)
path60.crude<-pathway.compare.nocomvar(pathtab=list(ref60.pathway),sampleid = "sample.id",
                                       mapfile=ref60.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path60.adjust<-pathway.compare(pathtab=list(ref60.pathway),sampleid = "sample.id",
                               mapfile=ref60.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                               comvar="Depression",adjustvar=c("Gender","Age","BMI"), longitudinal="no",
                               p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path60.crude,"metaR.ref60.pathway.crude.txt",quote = F,sep="\t")
write.table(path60.adjust,"metaR.ref60.pathway.adjust.txt",quote = F,sep="\t")


#######reference 65
ref65.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref65/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref65.pathway)
ref65.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref65/ref65.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref65.meta),paste("Ref65.",colnames(ref65.pathway),sep="")) != c(1:dim(ref65.pathway)[2]))
pid <- pmatch(rownames(ref65.meta),paste("Ref65.",colnames(ref65.pathway),sep=""))
ref65.pathway <- ref65.pathway[,pid]
dim(ref65.pathway)
dim(ref65.meta)
sum(pmatch(rownames(ref65.meta),paste("Ref65.",colnames(ref65.pathway),sep="")) != c(1:dim(ref65.pathway)[2]))
colnames(ref65.pathway) <- paste("Ref65.",colnames(ref65.pathway),sep="")

ref65.pathway.meta <- ref65.meta[,c("Sex","Age","BMI","Depression")]
ref65.pathway.meta$Depression <- gsub(1,"1_Depression",ref65.pathway.meta$Depression )
ref65.pathway.meta$Depression <- gsub(0,"0_Depression",ref65.pathway.meta$Depression )

colnames(ref65.pathway.meta)[1:4] <- c("Gender","Age","BMI","Depression")

ref65.pathway <- as.data.frame(t(ref65.pathway))
ref65.pathway.meta$sample.id <- rownames(ref65.pathway.meta)
path65.crude<-pathway.compare.nocomvar(pathtab=list(ref65.pathway),sampleid = "sample.id",
                                       mapfile=ref65.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path65.adjust<-pathway.compare(pathtab=list(ref65.pathway),sampleid = "sample.id",
                               mapfile=ref65.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                               comvar="Depression",adjustvar=c("Gender","Age","BMI"), longitudinal="no",
                               p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path65.crude,"metaR.ref65.pathway.crude.txt",quote = F,sep="\t")
write.table(path65.adjust,"metaR.ref65.pathway.adjust.txt",quote = F,sep="\t")


#####reference 77
ref77.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref77/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref77.pathway)
ref77.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref77/Ref77.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref77.meta),paste("Ref77.",colnames(ref77.pathway),sep="")) != c(1:dim(ref77.pathway)[2]))

ref77.pathway.meta <- cbind(ref77.meta[,c("Sex","Age","BMI","Depression")],ref77.meta)
ref77.pathway.meta$Depression <- gsub(1,"1_Depression",ref77.pathway.meta$Depression )
ref77.pathway.meta$Depression <- gsub(0,"0_Depression",ref77.pathway.meta$Depression )
colnames(ref77.pathway.meta)[1:4] <- c("Gender","Age","BMI","Depression")

colnames(ref77.pathway) <- paste("ref77.",colnames(ref77.pathway),sep="")
ref77.pathway <- as.data.frame(t(ref77.pathway))
ref77.pathway.meta$sample.id <- rownames(ref77.pathway.meta)
path77.crude<-pathway.compare.nocomvar(pathtab=list(ref77.pathway),sampleid = "sample.id",
                                       mapfile=ref77.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

path77.adjust<-pathway.compare(pathtab=list(ref77.pathway),sampleid = "sample.id",
                               mapfile=ref77.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                               comvar="Depression",adjustvar=c("Gender","Age","BMI"), longitudinal="no",
                               p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)

write.table(path77.crude,"metaR.ref77.pathway.crude.txt",quote = F,sep="\t")
write.table(path77.adjust,"metaR.ref77.pathway.adjust.txt",quote = F,sep="\t")


########reference 79
ref79.pathway <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/Ref79/dada2/picrust2/feature-table.biom.tsv",row.names = 1,head=T,check.names = F)
colSums(ref77.pathway)
ref79.meta <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/ref79/ref79.meta.dada2.MMUPHin.txt",sep="\t",header = T,check.names = F)
sum(pmatch(rownames(ref79.meta),paste("Ref79.",colnames(ref79.pathway),sep="")) != c(1:dim(ref79.pathway)[2]))
pid <- pmatch(rownames(ref79.meta),paste("Ref79.",colnames(ref79.pathway),sep="")) 
ref79.pathway <- ref79.pathway[,pid]

ref79.pathway.meta <- ref79.meta[,c("Depression"),drop=F]
ref79.pathway.meta$Depression <- gsub(1,"1_Depression",ref79.pathway.meta$Depression )
ref79.pathway.meta$Depression <- gsub(0,"0_Depression",ref79.pathway.meta$Depression )

colnames(ref79.pathway) <- paste("ref79.",colnames(ref79.pathway),sep="")
ref79.pathway <- as.data.frame(t(ref79.pathway))
ref79.pathway.meta$sample.id <- rownames(ref79.pathway.meta)


path79.crude<-pathway.compare.nocomvar(pathtab=list(ref79.pathway),sampleid = "sample.id",
                                       mapfile=ref79.pathway.meta,pathsum="rel", stat.med="gamlss",personid = "sample.id",
                                       comvar="Depression",adjustvar=NULL, longitudinal="no",
                                       p.adjust.method="fdr", percent.filter=0.05,relabund.filter=0.00005)
path79.adjust <- path79.crude


write.table(path79.crude,"metaR.ref79.pathway.crude.txt",quote = F,sep="\t")
write.table(path79.adjust,"metaR.ref79.pathway.adjust.txt",quote = F,sep="\t")

####meta analysis
###use all study

source("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/tmp/script_final/metaR.meta.niceplot.revised.R")
file.id <- ls(pattern="path\\d*.crude$")
dat.pathway.cb <- c()
for(i in file.id){
    txt <- paste("dat <- ",i,"$l1",sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat))
    )
    
    dat.pathway.cb <- dplyr::bind_rows(dat.pathway.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.pathway.cb,2,function(x){sum(!is.na(x))})


colnames(dat.pathway.cb)[colnames(dat.pathway.cb)=="i"]="study"

id_title <- read.table("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Final/metaanalysis/ID2Title_meta_analysis_study.txt",sep="\t",head=T)
studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.pathway.cb$study))
dat.pathway.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.pathway.cb$pop)

dat.pathway.cb$study <- dat.pathway.cb$pop
table(dat.pathway.cb$study)

meta.depress.crude.rr <- meta.taxa(taxcomdat=dat.pathway.cb,se.pattern="Std..Error.",
                                        summary.measure="RR", pool.var="id", studylab="study",
                                        backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")


metadat.random<-metatab.show(metatab=meta.depress.crude.rr$random,com.pooled.tab=dat.pathway.cb,sumvar="path",showvar="Depression1_Depression",p.cutoff.type="p", p.cutoff=0.05,display="data")
all.crude.mark <- metadat.random$taxsig$id
p.crude.all<-meta.niceplot.rev(metadat=metadat.random,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)",">=0.5"),leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")


####remove ref65
dat.pathway.rm65.cb <- c()
file.id <- file.id[-6]
for(i in file.id){
    txt <- paste("dat <- ",i,"$l1",sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat))
    )
    
    dat.pathway.rm65.cb <- dplyr::bind_rows(dat.pathway.rm65.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.pathway.rm65.cb,2,function(x){sum(!is.na(x))})

colnames(dat.pathway.rm65.cb)[colnames(dat.pathway.rm65.cb)=="i"]="study"
dat.pathway.rm65.cb$pop <- dat.pathway.rm65.cb$study
table(dat.pathway.rm65.cb$study)

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.pathway.rm65.cb$study))
dat.pathway.rm65.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.pathway.rm65.cb$pop)

dat.pathway.rm65.cb$study <- dat.pathway.rm65.cb$pop
table(dat.pathway.rm65.cb$study)

meta.depress.crude.rm65.rr <- meta.taxa(taxcomdat=dat.pathway.rm65.cb,se.pattern="Std..Error.",
                                        summary.measure="RR", pool.var="id", studylab="study",
                                        backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")


#pdf("function/metaR.pathway.rm65.occ0.5.p005.gamlss.crude.random.ef.pdf",width = 15,height = 10)
metadat.random<-metatab.show(metatab=meta.depress.crude.rm65.rr$random,com.pooled.tab=dat.pathway.rm65.cb,sumvar="path",showvar="Depression1_Depression",p.cutoff.type="p", p.cutoff=1,display="data")

all.crude.mark[!all.crude.mark %in% metadat.random$taxsig$id]

metadat.random$taxsig <- metadat.random$taxsig[metadat.random$taxsig$id %in% all.crude.mark,]
metadat.random$taxsig.all <- metadat.random$taxsig.all[metadat.random$taxsig.all$id %in% all.crude.mark,]


p.crude.rm65 <- meta.niceplot.rev(metadat=metadat.random,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)",">=0.5"),leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")
#dev.off()

############################
#########adjust model###################
############################
###use all study
file.id <- ls(pattern="path\\d*.adjust$")
colnames(path65.adjust$l1) <- gsub("GenderMale","Gendermale",colnames(path65.adjust$l1))
dat.pathway.cb <- c()
for(i in file.id){
    txt <- paste("dat <- ",i,"$l1",sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat)),
                   grep("Age",colnames(dat)),
                   grep("BMI",colnames(dat)),
                   grep("Gender",colnames(dat))
                   
                   
    )
    
    dat.pathway.cb <- dplyr::bind_rows(dat.pathway.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.pathway.cb,2,function(x){sum(!is.na(x))})

colnames(dat.pathway.cb)[colnames(dat.pathway.cb)=="i"]="study"
dat.pathway.cb$pop <- dat.pathway.cb$study
table(dat.pathway.cb$study)

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.pathway.cb$study))
dat.pathway.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.pathway.cb$pop)

dat.pathway.cb$study <- dat.pathway.cb$pop
table(dat.pathway.cb$study)

meta.depress.adjust.rr <- meta.taxa(taxcomdat=dat.pathway.cb,se.pattern="Std..Error.",
                                    summary.measure="RR", pool.var="id", studylab="study",
                                    backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")

#pdf("function/metaR.pathway.occ0.5.p005.gamlss.adjust.random.ef.pdf",width = 15,height = 10)
metadat.random<-metatab.show(metatab=meta.depress.adjust.rr$random,com.pooled.tab=dat.pathway.cb,sumvar="path",showvar="Depression1_Depression",p.cutoff.type="p", p.cutoff=1,display="data")

metadat.random$taxsig <- metadat.random$taxsig[metadat.random$taxsig$id %in% all.crude.mark,]
metadat.random$taxsig.all <- metadat.random$taxsig.all[metadat.random$taxsig.all$id %in% all.crude.mark,]


p.adjust.all <- meta.niceplot.rev(metadat=metadat.random,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)",">=0.5"),leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")


##remove ref65
colnames(path65.adjust$l1) <- gsub("GenderMale","Gendermale",colnames(path65.adjust$l1))
dat.pathway.rm65.cb <- c()
file.id<-file.id[-6]
for(i in file.id){
    txt <- paste("dat <- ",i,"$l1",sep="")
    eval((parse(text =txt)))
    coln.match <-c(1,grep("Depression",colnames(dat)),
                   grep("Age",colnames(dat)),
                   grep("BMI",colnames(dat)),
                   grep("Gender",colnames(dat))
                   
                   
    )
    
    dat.pathway.rm65.cb <- dplyr::bind_rows(dat.pathway.rm65.cb,cbind.data.frame(dat[,coln.match],i))
}
apply(dat.pathway.rm65.cb,2,function(x){sum(!is.na(x))})

colnames(dat.pathway.rm65.cb)[colnames(dat.pathway.rm65.cb)=="i"]="study"
dat.pathway.rm65.cb$pop <- dat.pathway.rm65.cb$study
table(dat.pathway.rm65.cb$study)

studyID <-  as.numeric(gsub(".*?([0-9]+).*", "\\1", dat.pathway.rm65.cb$study))
dat.pathway.rm65.cb$pop <- id_title[match(studyID,id_title$id),"title"]
table(dat.pathway.rm65.cb$pop)

dat.pathway.rm65.cb$study <- dat.pathway.rm65.cb$pop
table(dat.pathway.rm65.cb$study)

meta.depress.adjust.rm65.rr <- meta.taxa(taxcomdat=dat.pathway.rm65.cb,se.pattern="Std..Error.",
                                         summary.measure="RR", pool.var="id", studylab="study",
                                         backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")


#pdf("function/metaR.pathway.rm65.occ0.5.p005.gamlss.adjust.random.ef.pdf",width = 15,height = 10)
metadat.random<-metatab.show(metatab=meta.depress.adjust.rm65.rr$random,com.pooled.tab=dat.pathway.rm65.cb,sumvar="path",showvar="Depression1_Depression",p.cutoff.type="p", p.cutoff=1,display="data")

metadat.random$taxsig <- metadat.random$taxsig[metadat.random$taxsig$id %in% all.crude.mark,]
metadat.random$taxsig.all <- metadat.random$taxsig.all[metadat.random$taxsig.all$id %in% all.crude.mark,]

p.adjust.rm65 <-  meta.niceplot.rev(metadat=metadat.random,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)",">=0.5"),leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")
#dev.off()

sum(as.character(p.crude.all$forest$data$id) != as.character(p.crude.rm65$forest$data$id))
sum(as.character(p.crude.all$forest$data$id) != as.character(p.adjust.all$forest$data$id))
sum(as.character(p.crude.all$forest$data$id) != as.character(p.adjust.rm65$forest$data$id))

library(ggpubr)

ggarrange(p.crude.all$heatmap,p.crude.all$forest,p.crude.rm65$forest,
          p.adjust.all$forest,p.adjust.rm65$forest,
          ncol = 5,nrow = 1,align = "hv",heights=c(1,1),widths=c(1.8,2,2),legend="right",common.legend=T,legend.grob = get_legend(p.adjust.all))
dev.off()
write.table(all.crude.mark,"function/marker.list.txt",sep="\t",quote = F)
library(xlsx)


##add information
pathway.info <- read.csv("function/pathway_info.csv")
head(pathway.info)
sum(pmatch(metadat.random$taxsig$id,gsub("-",".",pathway.info$Pathway.ID)) != 1:64)

metadat.random$taxsig$id <- pathway.info[pmatch(metadat.random$taxsig$id,gsub("-",".",pathway.info$Pathway.ID)),"Pathway.information"]
metadat.random$taxsig.all$id <- pathway.info[match(metadat.random$taxsig.all$id,gsub("-",".",pathway.info$Pathway.ID)),"Pathway.information"]

p.adjust.rm65 <-  meta.niceplot.rev(metadat=metadat.random,sumtype="path",p="p",p.adjust="p.adjust",p.sig.heat="yes",heat.forest.width.ratio =c(1,1.3),est.break = c(-Inf, -0.5,-0.1,-0.05,0,0.05,0.1,0.5, Inf),est.break.label = c("<-0.5)", "[-0.5,-0.1)","[-0.1,-0.05)","[-0.05,0)","[0,0.05)","[0.05,0.1)", "[0.1,0.5)",">=0.5"),leg.key.size=0.8, leg.text.size=10, heat.text.x.size=10, forest.axis.text.y=8,forest.axis.text.x=10, point.ratio = c(4,2),line.ratio = c(2,1),neg.palette="Oranges" , pos.palette="PuBu")

pdf("function/metaR.allstudy.pathway.cb.pdf",width = 35,height = 15,onefile=FALSE)
ggarrange(p.crude.all$heatmap,p.crude.all$forest,p.crude.rm65$forest,
          p.adjust.all$forest,p.adjust.rm65$forest,
          ncol = 5,nrow = 1,align = "hv",heights=c(1,1),widths=c(1.8,2,2),legend="right",common.legend=T,legend.grob = get_legend(p.adjust.all)) 
dev.off()

###consistent
p.crude.mark <- p.crude.all$forest$data[p.crude.all$forest$data$p <0.05 & p.crude.all$forest$data$estimate >0,"taxa"]
p.crude.rm65.mark <- p.crude.rm65$forest$data[p.crude.rm65$forest$data$p <0.05 & p.crude.rm65$forest$data$estimate >0,"taxa"]

p.adjust.mark <- p.adjust.all$forest$data[p.adjust.all$forest$data$p <0.05 & p.adjust.all$forest$data$estimate >0,"taxa"]
p.adjust.rm65.mark <- p.adjust.rm65$forest$data[p.adjust.rm65$forest$data$p <0.05 & p.adjust.rm65$forest$data$estimate >0,"taxa"]

commen.path <- intersect(intersect(p.crude.mark,p.crude.rm65.mark),p.adjust.mark)
commen.path.info <- pathway.info[pmatch(commen.path,gsub("-",".",pathway.info$Pathway.ID)),]
intersect(commen.path.info$Pathway.information,p.adjust.rm65.mark)
