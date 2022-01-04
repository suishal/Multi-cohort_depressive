rm(list = ls())
setwd("D:/OneDrive - connect.hku.hk/project/systematic_review_and_meta_analysis/paper/")
library(grDevices)
color <- c(rgb(255,183,3,maxColorValue=255),rgb(33,158,188,maxColorValue=255))
library(xlsx)
human.method.all<-read.xlsx("Stable1.characteristic of 29 studies.xlsx",sheetIndex = 1,header = T,rowIndex=1:30)

human.method <- human.method.all[human.method.all$Extracted.data.from != "\\",]
##################################
##################################
##################################
##summary based on the unique study
hist(human.method$Year)

table(human.method$Study.design)
table(paste(human.method$Study.design,human.method$Study.design2))

table(human.method$Country)
human.method$questionnaire.abb
question <- c()
for(j in 1:dim(human.method)[1]){
  if(!is.na(human.method[j,"questionnaire.abb"])){
    question.t<-unlist(strsplit(human.method[j,"questionnaire.abb"],split='; '))
    question<-c(question,paste(question.t,human.method[j,"ID"],sep=":"))
  }
}
question.m<-data.frame(t(unlist(sapply(question, function(x){unlist(strsplit(x,split = ":"))}, simplify = TRUE))))
colnames(question.m) <- c("question","ID")
library(ggplot2)
question.table <- data.frame(table(question.m[,"question"]))
question.table$Var1 <- factor(question.table$Var1, levels = question.table[order(question.table$Freq,decreasing = T),"Var1"])
pie.question <- ggplot(question.table, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)
pie <- pie.question + coord_polar("y", start=0)

question <- c()
for(j in 1:dim(human.method)[1]){
  if(!is.na(human.method[j,"questionnaire.abb"])){
    question.t<-unlist(strsplit(human.method[j,"questionnaire.abb"],split='; '))
    question<-c(question,paste(question.t,human.method[j,"ID"],sep=":"))
  }
}
question.m<-data.frame(t(unlist(sapply(question, function(x){unlist(strsplit(x,split = ":"))}, simplify = TRUE))))
colnames(question.m) <- c("question","ID")
question.m$question <- as.character(question.m$question)
question.family <- sapply(question.m$question, function(x){strsplit(x,split = "_")},simplify = "array")
question.m$family <- sapply(question.family, "[", 1)
question.m$sub_family <- sapply(question.family, "[", 2)
question.m$sub_family[!(question.m$family %in% c("HAMDS","BDI","DASS","Interview","DSM"))]=0
question.m$sub_family[is.na(question.m$sub_family)]="Not mention"

question.table <- data.frame(table(question.m[,"question"]))
question.table$Var1 <- factor(question.table$Var1, levels = question.table[order(question.table$Freq,decreasing = T),"Var1"])
pie.question <- ggplot(question.table, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)
pie <- pie.question + coord_polar("y", start=0)
ggplot(aes(x=Var1,y=Freq),data=question.table)+geom_bar(stat="identity")+ coord_flip()

question.table <- data.frame(table(question.m[,"family"]))
question.table$Var1 <- factor(question.table$Var1, levels = question.table[order(question.table$Freq,decreasing = T),"Var1"])
pie.question <- ggplot(question.table, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)
pie <- pie.question + coord_polar("y", start=0)
ggplot(aes(x=Var1,y=Freq),data=question.table)+geom_bar(stat="identity")+coord_flip()

pdf("questionnaire.stat.pdf",width = 8,height = 4)
st <- table(question.m$family)
question.m$family <- factor(question.m$family,levels = names(st)[order(st,decreasing = T)])
ggplot(data=question.m,aes(x=family)) + geom_bar(aes(fill=as.factor(sub_family)),position="stack")+coord_flip()+scale_y_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10))+theme_classic()
dev.off()

#country summary
table(human.method$Country)

#country and questionnaire
question.m$country <- human.method[match(question.m$ID,human.method$ID),"Country"]
table(question.m$family,question.m$country)

continent<-data.frame(table(question.m$country))
continent$continent <- c("Australia","Europe","Europe","Asia","Europe","Europe","Europe","Europe","Asia","Europe","North America")

question.m$continent <- continent[match(question.m$country,continent$Var1),"continent"]
table(question.m$family,question.m$continent)

##sequence method
table(human.method$Sequence.method)


#############
#############
###Taxonomy summary
human.result <- read.xlsx("Stable2.result of 24 studies_rmformate.xlsx",sheetIndex = 1,header = T,rowIndex = 1:35)
colnames(human.result)
tax.increasd<-c()
increase.list<-colnames(human.result)[16:20]
for(i in increase.list){
  human.result[,i]<-as.character(human.result[,i])
  for(j in 1:dim(human.result)[1]){
    if(!is.na(human.result[j,i])){
      tax.t<-unlist(strsplit(human.result[j,i],split='; '))
      tax.increasd<-c(tax.increasd,paste(tax.t,human.result[j,"ID"],human.result[j,"sub_id"],sep=":"))
    }
  }
}

tax.increasd.m<-t(sapply(tax.increasd, function(x){unlist(strsplit(x,split = ":"))}, simplify = TRUE))
tax.increasd.m[order(tax.increasd.m[,1]),]
tax.increasd.m<-as.data.frame(tax.increasd.m)
colnames(tax.increasd.m) <- c("Tax","ref")
increased.tax.stat<-table(tax.increasd.m$Tax)


tax.decreasd<-c()
decrease.list<-colnames(human.result)[21:25]
for(i in decrease.list){
  human.result[,i]<-as.character(human.result[,i])
  for(j in 1:dim(human.result)[1]){
    if(!is.na(human.result[j,i])){
      tax.t<-unlist(strsplit(human.result[j,i],split='; '))
      tax.decreasd<-c(tax.decreasd,paste(tax.t,human.result[j,"ID"],human.result[j,"sub_id"],sep=":"))
    }
  }
}

tax.decreasd.m<-t(unlist(sapply(tax.decreasd, function(x){unlist(strsplit(x,split = ":"))}, simplify = TRUE)))

colnames(tax.decreasd.m)<-c("Tax","Ref")
tax.decreasd.m<-as.data.frame(tax.decreasd.m)
decreased.tax.stat<-table(tax.decreasd.m$Tax)

decrease.stat<-as.data.frame(table(tax.decreasd.m$Tax))
increase.stat<-as.data.frame(table(tax.increasd.m$Tax))

length(unique(c(decrease.stat$Var1,increase.stat$Var1)))

conflict.tax<-merge(decrease.stat,increase.stat,by = "Var1")
decrease.pure<-decrease.stat[!(decrease.stat$Var1 %in% conflict.tax$Var1),]
increase.pure<-increase.stat[!(increase.stat$Var1 %in% conflict.tax$Var1),]

dim(conflict.tax)
dim(decrease.pure)
dim(increase.pure)
sum(decrease.pure$Freq>1) + sum(increase.pure$Freq>1)

length(unique(c(as.character(decrease.pure$Var1),as.character(increase.pure$Var1))))
length(c(as.character(decrease.pure$Var1),as.character(increase.pure$Var1)))

dim(conflict.tax)[1]/length(unique(c(decrease.stat$Var1,increase.stat$Var1)))
dim(conflict.tax)

library(ggplot2)

library(tidyverse)
decrease.pure.m3 <- decrease.pure[decrease.pure$Freq >=3,]
p.d <- ggplot(data=decrease.pure.m3, aes(x=reorder(Var1,Freq), y=Freq)) +
  geom_bar(stat="identity", fill=color[1])+
  geom_text(aes(label=Freq), color="white", size=3, hjust=1.6)+
  theme_classic()+
  labs(title="Decrease Taxonomy", 
       x="Taxonomy", y = "Related reference")
p.d <- p.d + coord_flip()

increase.pure.m3 <- increase.pure[increase.pure$Freq >=3,]
p.i <- ggplot(data=increase.pure.m3, aes(x=reorder(Var1,Freq), y=Freq)) +
  geom_bar(stat="identity", fill=color[2])+
  geom_text(aes(label=Freq), color="white", size=3, hjust=1.6)+
  theme_classic()+
  labs(title="Increase Taxonomy", 
       x="Taxonomy", y = "Related reference")
p.i <- p.i + coord_flip()

colnames(conflict.tax) <- c("Tax","Decrease in MDD","Increase in MDD")
conflict.tax <- conflict.tax[conflict.tax$`Decrease in MDD` >3 | conflict.tax$`Increase in MDD` >3,]
conflict.tax <- conflict.tax[order(rowSums(conflict.tax[,2:3])),]
library(reshape2)
conflict.tax.m<-melt(conflict.tax)
conflict.tax.m$Tax <- factor(conflict.tax.m$Tax,levels =conflict.tax$Tax)
p.conflict.tax <- ggplot(data=conflict.tax.m, aes(x=Tax, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=value), hjust=1.6, color="white",
            position = position_dodge(0.9), size=3)+
  scale_fill_brewer(palette="Paired")+
  theme_classic() +
  labs(title="Conflict Taxonomy", 
       x="Taxonomy", y = "Related reference")+scale_fill_manual(values=color)
p.conflict.tax.enrich <- p.conflict.tax <- p.conflict.tax +coord_flip()
library(ggpubr)
pdf("Taxonomy.result.summary.pdf",width = 7,height = 8)
ggarrange(ggarrange(p.d,p.i,ncol = 2,labels = c("B", "C")),
          p.conflict.tax,nrow = 2,labels = "D",heights  = c(1.5,3),legend="bottom")
dev.off()

#alpha diversity
colnames(human.result)
alpha <- human.result[,c(1:2,36,5:14)]
AID <- alpha$ID
alpha$Title_s[!is.na(alpha$sub.group)] <- paste(alpha[!is.na(alpha$sub.group),"Title_s"],alpha[!is.na(alpha$sub.group),"sub.group"],sep=":")
rownames(alpha) <- alpha$Title_s
str(alpha)
#study.d <- alpha[,3,drop = F]
alpha <- alpha[,-c(1:2)]
apply(alpha[,1:10],1,function(x){table(x)})
dim(alpha)

##seperate significant and enrichment
alpha.sig<-matrix("",nrow = nrow(alpha),ncol = ncol(alpha),dimnames=list(rownames(alpha),colnames(alpha)))
alpha.dir<-matrix("",nrow = nrow(alpha),ncol = ncol(alpha),dimnames=list(rownames(alpha),colnames(alpha)))
for( i in 1:nrow(alpha)){
  for(j in 1:ncol(alpha)){
    if(is.na(alpha[i,j])){
      #alpha.dir[i,j] <-0
    }
    else{
      if(alpha[i,j] == "Con" | alpha[i,j] == "MDD" | alpha[i,j] == "Treat"){
        alpha.sig[i,j] <- "*"
        alpha.dir[i,j] <- as.character(alpha[i,j])
      }
      if(grepl("NS",alpha[i,j])){
        alpha.sig[i,j]="NS"
        alpha.dir[i,j]<-""
      }
      if(grepl("NS; ",alpha[i,j])){
        alpha.sig[i,j]="NS"
        alpha.dir[i,j]<-gsub("NS; ","",alpha[i,j])
      }
    }
  }
}

apply(alpha,2,function(x){sum(grepl("NS",x))})
apply(alpha.sig,2,function(x){sum(grepl("NS",x))})

apply(alpha,2,function(x){sum(grepl("Con",x))})
apply(alpha.dir,2,function(x){sum(grepl("Con",x))})

apply(alpha,2,function(x){sum(grepl("MDD",x))})
apply(alpha.dir,2,function(x){sum(grepl("MDD",x))})


dist_letters = function(x, y) {
  x = strtoi(charToRaw(paste(x, collapse = "")), base = 16)
  y = strtoi(charToRaw(paste(y, collapse = "")), base = 16)
  sqrt(sum((x - y)^2))
}
library(ComplexHeatmap)
pdf("Alpha_diversity_heatmap.pdf",width = 6,height = 9)

##index
alpha.sig.index <- alpha.sig[,-2]
alpha.dir.indes <- alpha.dir[,-2]
h.alpha1 <- Heatmap(alpha.dir.indes, name = "letters", col = structure(c("gray",color[2],color[2],color[1],color[1],color[1]), names = c("NA","MDD","BPD","Con","BD_treat","Treat")),
                   clustering_distance_rows = dist_letters, clustering_distance_columns = dist_letters,
                   cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                     grid.text(alpha.sig.index[i, j], x, y)
                   })

h.alpha2 <- Heatmap(alpha.dir[,1], name = "letters", col = structure(c("gray",color[2],color[2],color[1],color[1],color[1]), names = c("NA","MDD","BPD","Con","BD_treat","Treat")),
                    clustering_distance_rows = dist_letters, clustering_distance_columns = dist_letters,
                    cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                      grid.text(alpha.sig[i, j], x, y)
                    })
h.alpha
dev.off()

###beta diversity
table(human.result$Overall.alpha.diversity)
table(human.result$Beta.diversity)
