###Review 60
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
library(XML)
library("methods")
setwd("/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/Analysis/Re60")

sample1 <- xmlParse(file = "../../DataDownload/Ref60/EGAD00001004449/xmls/samples/EGAN00001892232.sample.xml")
rootnode <- xmlToDataFrame(sample1)
xmlSize(rootnode)
rootnode[1]
xmlToList(sample1)
SAMPLE_ATTRIBUTES

metadata.1 <- xmlToDataFrame(nodes = getNodeSet(sample1, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"))

fold.list <- list.files(path = "../../DataDownload/Ref60/EGAD00001004449/xmls/samples/")

metadata.t1 <- data.frame(Doubles=rep(as.double(NA),19) )
metadata.t2 <- data.frame(Doubles=rep(as.double(NA),6) )
for(i in fold.list){
    sample.file <- xmlParse(file = paste("../../DataDownload/Ref60/EGAD00001004449/xmls/samples/",i,sep=""))
    metadata.tmp <- xmlToDataFrame(nodes = getNodeSet(sample.file, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"))
    if(dim(metadata.tmp)[1]==19){
        metadata.t1 <- cbind.data.frame(metadata.t1,metadata.tmp)
    }
    if(dim(metadata.tmp)[1]==6){
        metadata.t2 <- cbind.data.frame(metadata.t2,metadata.tmp)
    }
    
}
dim(metadata.t1)
dim(metadata.t2)
n <-c(1, seq(2, 2066, 2))
metadata.d <- metadata.t1[,-n]
rownames(metadata.d) <- metadata.t1[,2]
metadata.d <- t(metadata.d)
rownames(metadata.d) <- metadata.d[,"subject_id"]
metadata.d[,"antidepressants"] <- gsub("FALSE",0,metadata.d[,"antidepressants"])
metadata.d[,"antidepressants"] <- gsub("TRUE",1,metadata.d[,"antidepressants"])

metadata.d[,"depression"] <- gsub("FALSE",0,metadata.d[,"depression"])
metadata.d[,"depression"] <- gsub("TRUE",1,metadata.d[,"depression"])
write.table(metadata.d,"metadata_FGFP.txt",sep="\t",quote = F)

##propensity score
metadata.d <- as.data.frame(metadata)

metadata.d[,"depression"] <- as.numeric(metadata.d[,"depression"])
metadata.d[,"age"] <- as.numeric(metadata.d[,"age"])
metadata.d[,"antidepressants"] <- as.numeric(metadata.d[,"antidepressants"])
metadata.d[,"bmi"] <- as.numeric(metadata.d[,"bmi"])
metadata.d[,"bss"] <- as.numeric(metadata.d[,"bss"])
m_ps <- glm(depression ~ bmi + bss + age + gender,
            family = binomial(), data = metadata.d)
summary(m_ps)
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     depression = m_ps$model$depression)
head(prs_df)

labs <- paste("Depression:", c("Yes", "No"))
prs_df %>% 
    mutate(depression = ifelse(depression == 1, labs[1], labs[2])) %>%
    ggplot(aes(x = pr_score)) +
    geom_histogram(color = "white") +
    facet_wrap(~depression) +
    xlab("Probability of depression") +
    theme_bw()
library(MatchIt)
mod_match <- matchit(depression ~ bmi + bss + age + gender,
                     method = "nearest", data = metadata.d)
dta_m <- match.data(mod_match)
dim(dta_m)
write.table(dta_m,"match_sample.metadata.txt",sep="\t",quote = F)

