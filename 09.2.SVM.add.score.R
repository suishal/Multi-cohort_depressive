rm(list=ls())
source("~/Project/Review/meta-analysis/taxname.reformat.R")
genus <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.count.genus.prof",sep="\t",check.names = F)
metadata <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.metadata.prof",sep="\t",check.names = F)

noise.removal <- function(data, percent=0.1, top=NULL){
  matrix <- data
  rmnoise.bac <- apply(matrix, 1, function(x){sum(x>0)/dim(matrix)[2]}) > percent 
  matrix_1 <- matrix[rmnoise.bac,]
  return(matrix_1)
}

dim(genus)
genus.denoinzed <- noise.removal(genus)
dim(genus.denoinzed)

occ.genus <- apply(genus.denoinzed, 1, function(x){sum(x>0)/dim(genus.denoinzed)[2]})
###remove unclassfied genus

rownames(genus.denoinzed)
genus.denoinzed <- genus.denoinzed[-1,]

##feature selection
library(pROC)
library(e1071)
library(plotROC)
library(mlbench)
library(caret)

genus.denoinzed <- t(genus.denoinzed)
correlationMatrix <- cor(genus.denoinzed)
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attributes
print(highlyCorrelated)
genus.denoinzed <- genus.denoinzed[,-highlyCorrelated]

pid <- pmatch(rownames(metadata),rownames(genus.denoinzed))
sum(is.na(pid))
genus.denoinzed <- as.data.frame(genus.denoinzed[pid,])
genus.denoinzed$Depression <-metadata$Depression
genus.denoinzed$Depression <- as.factor(genus.denoinzed$Depression)

# rank features by importance using the caret r packageR
# ensure results are repeatable
set.seed(1)
data_split <- createDataPartition(y=genus.denoinzed$Depression, p=0.8, list=F) #splits data
training <- genus.denoinzed[data_split[,1],] #call training data
testing <- genus.denoinzed[-data_split[,1],] #call testing or validation data

###metadata of train and test dataset
training.meta <- metadata[data_split[,1],]
testing.meta <- metadata[-data_split[,1],]

fisher.test(rbind(
table(training.meta$subjectID),table(testing.meta$subjectID)),simulate.p.value =  T)


# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=10, repeats=10)
# # train the model
# set.seed(1)
# model <- train(Depression~., data=training, method="lvq", preProcess="scale", trControl=control)
# 
# roc_imp <- filterVarImp(x = training[, -ncol(training)], y = training$Depression)
# head(roc_imp)
# # estimate variable importance
# importance <- varImp(model, scale=TRUE)
# # summarize importance
# print(importance,top=30)
# # plot importance
# plot(importance)


control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10)
# run the RFE algorithm
x <- scale(training[,-ncol(training)])
set.seed(1)
results_scale <- rfe(x, training[,"Depression"], sizes=c(1:50), rfeControl=control)
vimp <- varImp(results_scale,scale=FALSE)
vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
vimp[1:30,,drop=F]             

 #save(model,file="svm/RFE.feature.selection.model.rda")
 #save(vimp,file="svm/RFE.feature.selection.rda")
# summarize the results
print(results_scale)
# list the chosen features
feature.selected.scale <- predictors(results_scale)
##write.table(feature.selected.scale,file="svm/selected.top30.feature.txt",quote = F,sep="\t",row.names = F,col.names = F)
# plot the results
#pdf("svm/SVM_feature_selection.pdf",height = 6,width = 6)
plot(results_scale, type=c("g", "o"))
dev.off()

n = results_scale$bestSubset #n=30 have the highest accuracy

importance <- vimp[rownames(vimp) %in% feature.selected.scale,,drop=F]
importance$var <- extract.genus(rownames(importance))[,"genus"]
#pdf("svm/top30.importance.feature.pdf",width = 7,height = 5)
 ggdotchart(importance, x = "var", y = "Overall",
           color = "#17c3b2",
           sorting = "desc",                       
           add = "segments",                            
           add.params = list(color = "lightgray", size = 1), 
           dot.size = 3, xlab = "",ylab="Importance score"  ,             
           ggtheme = theme_classic2() 
) + font("xy.text", size = 10, vjust = 0.5,face = "bold") 

dev.off()

selected.train <- training[,c(feature.selected.scale,"Depression")]
selected.test <- testing[,c(feature.selected.scale,"Depression")]

set.seed(1)
train.m=svm(x=as.matrix(selected.train[,-ncol(selected.train)]),y=as.factor(selected.train$Depression), gamma = 0.3, probability = T,scale=T)
train.class=predict(train.m,selected.train[,-ncol(selected.train)],probability = T)
x=attributes(train.class)
roc.train=pROC::roc(selected.train$Depression~x$probabilities[,1])

test.class=predict(train.m,selected.test[,-ncol(selected.test)],probability = T)
x=attributes(test.class)
roc.test=pROC::roc(selected.test$Depression~x$probabilities[,1])

plot(roc.train,col="blue")
plot(roc.test,col="red",add=T)

### parameters training
library(e1071)

set.seed(1)
svm.tune.scale <- tune.svm(Depression~.,
                           data = selected.train,degree = c(2:10), cost=c(1:100),gamma=10^(-6:1),scale=T)
summary(svm.tune.scale)

best.para <- svm.tune.scale$best.parameters
set.seed(1)
svm.fit.train_scale = svm(Depression~., 
                    data = selected.train, 
                    degree=best.para[1],
                    gamma=best.para[2],
                    cost=best.para[3],
                    scale = T,probability = TRUE)

set.seed(1)
svm.fit.train_scale = svm(Depression~.,
  data = selected.train,
  degree=2,
  gamma=0.1,
  cost=1,
  scale = T,probability = TRUE)

##save(svm.tune.scale,file = "svm/svm.tune.rda")

summary(svm.fit.train_scale)
#save(svm.fit.train_scale,file = "svm/final.svm.rda")
x <- subset(selected.train, select = -Depression)
y <- selected.train$Depression

library(ROCR)
svm.fit.train.predict <- predict(svm.fit.train_scale,x, probability = TRUE)
n<-ifelse(svm.fit.train.predict==y,1,0)
sum(n)/dim(selected.train)[1]

probability.train <- attr(svm.fit.train.predict, "probabilities")
pr.train<-prediction(probability.train[,1], selected.train$Depression)
prf.train<- performance(pr.train, measure="tpr",x.measure="fpr")
plot(prf.train,col="blue")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr.train, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

#####test dataset
svm.fit.test.predict <- predict(svm.fit.train_scale,selected.test[,-31],probability = TRUE)
n<-ifelse(svm.fit.test.predict==selected.test$Depression,1,0)
sum(n)/dim(selected.test)

probability.test <- attr(svm.fit.test.predict, "probabilities")
pr<-prediction(probability.test[,1], selected.test$Depression)
prf<- performance(pr, measure="tpr",x.measure="fpr")
plot(prf,col="green")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

#write.table(selected.train,"svm/train.dataset.txt",sep="\t",quote = F)
#write.table(selected.test,"svm/test.dataset.txt",sep="\t",quote = F)
#write.table(importance,"svm/feature.importance.txt",sep="\t",quote = F)
#pdf("svm/final.model.roc.pdf",width = 10,height = 5)
par(mfrow=c(1:2))
roc.train <- pROC::roc(as.numeric(selected.train$Depression),probability.train[,1])
plot(roc.train, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     max.auc.polygon=TRUE,
     auc.polygon.col="lightblue",main="SVM model in training set")

roc.test <- pROC::roc(selected.test$Depression,probability.test[,1])

plot(roc.test,print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     max.auc.polygon=TRUE,
     auc.polygon.col="lightyellow",main="SVM model in test set")

dev.off()
###selceted genus abundance
library(clusterSim)
genus.nor <- t(data.Normalization(genus,type = "n10",normalization ="column"))

pid <- pmatch(rownames(metadata),rownames(genus.nor))
sum(is.na(pid))
genus.nor <- as.data.frame(genus.nor[pid,])
genus.nor$Depression <-metadata$Depression
genus.nor$Depression <- as.factor(genus.nor$Depression)
genus.nor.select <-  genus.nor[,c(rownames(importance),"Depression")]
colnames(genus.nor.select)[1:30] <- extract.genus(colnames(genus.nor.select)[1:30] )[,"genus"]

library(ggpubr)
par(mfcol = c(5, 6))
for (i in 1:30) {
  dat.tmp <- cbind.data.frame(genus.nor.select[, i],genus.nor.select$Depression)
  colnames(dat.tmp) <- c("taxo","Depression")
  p <- ggboxplot(dat.tmp, x = "Depression", y = "taxo",
                 fill = "Depression", palette = c("#ffb703","#219ebc")) +
    #stat_compare_means(method = "kruskal.test", label.y = max(dat.tmp$taxo)+0.03, label.x = 2)+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "wilcox.test") +
    ylab(colnames(genus.nor.select)[i])+xlab("Depression")
  txt <- paste("p",i," <- p",sep="")
  eval((parse(text =txt)))
}

##mean and median heatmap
mean.matrix <- apply(genus.nor.select[1:30], 2, function(x){by(x,genus.nor.select$Depression,mean)})
median.matrix <- apply(genus.nor.select[1:30], 2, function(x){by(x,genus.nor.select$Depression,median)})
pvalue.matrix <- apply(genus.nor.select[1:30], 2, function(x){wilcox.test(x~genus.nor.select$Depression)$p.value})

mean.dir <- as.matrix(as.numeric(mean.matrix[1,]<mean.matrix[2,]))
rownames(mean.dir) <- colnames(mean.matrix)

p.sig <- as.matrix(rep(" ",length(pvalue.matrix)))
p.sig[pvalue.matrix<0.05,1] <- "*"
rownames(p.sig) <- colnames(mean.matrix)
id <- which(mean.dir != median.dir)
median.matrix[,id]
mean.matrix[,id]
##plot the mean and p value
col_fun = colorRamp2(c(0,1), c("#ffb703", "#219ebc"))

ht1<-Heatmap(mean.dir,name="Mean",cluster_rows=F,cluster_columns=F,col = col_fun)

ht1<-Heatmap(mean.dir,name="Enrichment",cluster_rows=F,cluster_columns=F,col = col_fun,
             layer_fun = function(j, i, x, y, width, height, fill) {
               v = pindex(p.sig, i, j)
               grid.text(v, x, y, gp = gpar(fontsize = 10))
             })

col_fun = colorRamp2(c(0,0.01,0.12),c("#d9ed92","#1a759f","#184e77"))
ht2 <- Heatmap(t(mean.matrix),name="Mean",cluster_rows=F,cluster_columns=F,col = col_fun)
ht2
ht_list = ht1 + ht2
pdf("svm/top30.importance.feature.heatmap.pdf",height = 7,width = 4 )
ht_list
dev.off()

##mean and median heatmap (study)
genus.nor.study <- genus.nor.select
genus.nor.study$study <- metadata$subjectID
mean.study.matrix <- aggregate(genus.nor.study[1:30], by=list(genus.nor.study$Depression,genus.nor.study$study),mean)
for(i in study.id){
aggregate(genus.nor.study[1:30], by=list(genus.nor.study$Depression,genus.nor.study$study),function(x){wilcox.test$p.value})
}


#####add age.bmi and gender
sum(pmatch(rownames(metadata),rownames(genus.denoinzed)) != c(1:1827))
genus.denoinzed.meta <- genus.denoinzed
genus.denoinzed.meta$Age <-metadata$Age
genus.denoinzed.meta$BMI <- metadata$BMI
genus.denoinzed.meta$Sex <- metadata$Sex
genus.denoinzed.meta$Depression <- metadata$Depression
table(genus.denoinzed.meta$Sex)

genus.denoinzed.meta$Sex <- gsub("^male","Male",genus.denoinzed.meta$Sex)
genus.denoinzed.meta$Sex <- gsub("female","Female",genus.denoinzed.meta$Sex)

genus.denoinzed.meta$Depression <- as.factor(genus.denoinzed.meta$Depression)

dim(genus.denoinzed.meta)

training <- genus.denoinzed.meta[data_split[,1],] #call training data
testing <- genus.denoinzed.meta[-data_split[,1],] #

training.rmna <- training[apply(training, 1, function(x){sum(is.na(x))}) == 0,]
testing.rmna <- testing[apply(testing, 1, function(x){sum(is.na(x))}) == 0,]

selected.train <- training.rmna[,c(feature.selected.scale,c("Age","BMI","Sex","Depression"))]
selected.test <- testing.rmna[,c(feature.selected.scale,c("Age","BMI","Sex","Depression"))]

sum(apply(selected.train, 1, function(x){sum(is.na(x))}) != 0)

set.seed(1)
svm.fit.train_meta = svm(Depression~., 
                         data = selected.train, 
                         degree=2,
                         gamma=0.1,
                         cost=1,
                         scale = T,probability = TRUE)

summary(svm.fit.train_meta)
#save(svm.fit.train_meta,file = "svm/final.svm.meta.rda")
x <- subset(selected.train, select = -Depression)
y <- selected.train$Depression

library(ROCR)
svm.fit.train.predict <- predict(svm.fit.train_meta,x, probability = TRUE)
n<-ifelse(svm.fit.train.predict==y,1,0)
sum(n)/dim(selected.train)[1]

probability.train <- attr(svm.fit.train.predict, "probabilities")
pr.train<-prediction(probability.train[,2], selected.train$Depression)
prf.train<- performance(pr.train, measure="tpr",x.measure="fpr")
plot(prf.train,col="blue")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr.train, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

#####test dataset
svm.fit.test.predict <- predict(svm.fit.train_meta,selected.test[,-34],probability = TRUE)
n<-ifelse(svm.fit.test.predict==selected.test$Depression,1,0)
sum(n)/dim(selected.test)

probability.test <- attr(svm.fit.test.predict, "probabilities")
pr<-prediction(probability.test[,2], selected.test$Depression)
prf<- performance(pr, measure="tpr",x.measure="fpr")
plot(prf,col="green")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

#write.table(selected.train,"svm/train.dataset.meta.txt",sep="\t",quote = F)
#write.table(selected.test,"svm/test.dataset.meta.txt",sep="\t",quote = F)

#pdf("svm/final.model.roc.meta.pdf",width = 10,height = 5)
par(mfrow=c(1:2))
roc.train <- pROC::roc(selected.train$Depression,probability.train[,1])
plot(roc.train, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     max.auc.polygon=TRUE,
     auc.polygon.col="lightblue",main="SVM model in training set")

roc.test <- pROC::roc(selected.test$Depression,probability.test[,1])

plot(roc.test,print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     max.auc.polygon=TRUE,
     auc.polygon.col="lightyellow",main="SVM model in test set")

dev.off()


###only age bmi and gender
selected.train <- training.rmna[,c("Age","BMI","Sex","Depression")]
selected.test <- testing.rmna[,c("Age","BMI","Sex","Depression")]

dim(training.rmna)
dim(testing.rmna)
sum(apply(selected.train, 1, function(x){sum(is.na(x))}) != 0)

set.seed(1)
svm.fit.train_meta_only = svm(Depression~., 
                         data = selected.train, 
                         degree=2,
                         gamma=0.1,
                         cost=1,
                         scale = T,probability = TRUE)

summary(svm.fit.train_meta_only)
#save(svm.fit.train_meta_only,file = "svm/final.svm.meta.only.rda")
x <- subset(selected.train, select = -Depression)
y <- selected.train$Depression

library(ROCR)
svm.fit.train.predict <- predict(svm.fit.train_meta_only,x, probability = TRUE)
n<-ifelse(svm.fit.train.predict==y,1,0)
sum(n)/dim(selected.train)[1]

probability.train <- attr(svm.fit.train.predict, "probabilities")
pr.train<-prediction(probability.train[,1], selected.train$Depression)
prf.train<- performance(pr.train, measure="tpr",x.measure="fpr")
plot(prf.train,col="blue")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr.train, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

#####test dataset
svm.fit.test.predict <- predict(svm.fit.train_meta_only,selected.test[,-31],probability = TRUE)
n<-ifelse(svm.fit.test.predict==selected.test$Depression,1,0)
sum(n)/dim(selected.test)

probability.test <- attr(svm.fit.test.predict, "probabilities")
pr<-prediction(probability.test[,1], selected.test$Depression)
prf<- performance(pr, measure="tpr",x.measure="fpr")
plot(prf,col="green")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]


########selected feature in each study
sum(pmatch(rownames(metadata),rownames(genus.denoinzed)) != c(1:1827))
genus.denoinzed.study <- genus.denoinzed
genus.denoinzed.study$study <-metadata$subjectID

genus.denoinzed.study <-  genus.denoinzed.study[,c(feature.selected.scale,"study","Depression")]
study.id <- names(table(genus.denoinzed.study$study))

auc.tmp <- c()
for( i in study.id){
  svm.fit.study<- predict(svm.fit.train_scale,genus.denoinzed.study[genus.denoinzed.study$study == i,],probability = TRUE)
  probability.test <- attr(svm.fit.study, "probabilities")
  pr<-prediction(probability.test[,2], genus.denoinzed.study[genus.denoinzed.study$study == i,"Depression"])
  auc_ROCR <- performance(pr, measure = "auc")
  auc.study <- auc_ROCR@y.values[[1]]
  auc.tmp <- c(auc.tmp,auc.study)
}

#### predict power
selected.data <- genus.denoinzed[,c(feature.selected.scale,"Depression")]
svm.fit.all.predict <- predict(svm.fit.train_scale,selected.data[,-31], probability = TRUE)

n<-ifelse(svm.fit.all.predict==selected.data$Depression,1,0)
sum(n)/dim(selected.data)[1]


probability.data <- attr(svm.fit.all.predict, "probabilities")
pr.data<-prediction(probability.data[,1], selected.data$Depression)
all.prf.data<- performance(pr.data, measure="tpr",x.measure="fpr")
#pdf("svm/ROC_at_each_study.pdf")
plot(all.prf.data,col="black")
lines(x = c(0,1), y = c(0,1),col="gray")

auc_ROCR <- performance(pr.data, measure = "auc")
all.study.auc <- auc_ROCR@y.values[[1]]

sum(pmatch(rownames(metadata),rownames(selected.data)) != c(1:1827))
selected.data$study <- metadata$subjectID
study.id <- names(table(genus.denoinzed.study$study))

col <- rainbow(8)
roc.study <- c()
k <- 1
for(i in study.id){
  probability.data.study <- probability.data[selected.data$study == i,]
  pr.data<-prediction(probability.data.study[,1], selected.data[selected.data$study == i,"Depression"])
  auc_ROCR <- performance(pr.data, measure = "auc")
  roc.study <-c(roc.study, auc_ROCR@y.values[[1]])
  
  prf.data<- performance(pr.data, measure="tpr",x.measure="fpr")
  plot(prf.data,add = TRUE,col=col[k])
  k=k+1

}

id_title <- read.table("ID2Title_meta_analysis_study.txt",sep="\t",head=T)
id_title$id <- paste("Ref",id_title$id,sep="")
id_title$col <- col
id_title$auc <- roc.study
id_title$n <- table(genus.denoinzed.study$study)
id_title <- rbind(id_title,c("","All sample","#000000",all.study.auc,dim(genus.denoinzed.study)[1]))
id_title$auc <- as.numeric(id_title$auc)
legend("bottomright",
       legend = paste(paste(id_title$title," (",id_title$n,")",sep=""),round(id_title$auc,2),sep=": "),
       fill=id_title$col )
dev.off()


###predicted score with depression score
#Ref1 (HADS)
probability.data.study.ref1 <- probability.data[selected.data$study == "Ref1",]
ref1.score <- metadata[pmatch(rownames(probability.data.study.ref1),rownames(metadata)),"depression_score"]
cor.ref1 <- cor.test(probability.data.study.ref1[,1],ref1.score)

#8 and 10 suggests a mild disorder
#≥11 is considered a clinically significant disorder
#reference : Comparison of two self-rating scales to detect depression: HADS and PHQ-9
ref1.score.c <- rep(0,length(ref1.score))
ref1.score.c[7<ref1.score & ref1.score<11] <- 1
ref1.score.c[ ref1.score >= 11] <- 2

cor.test(probability.data.study.ref1[,1],ref1.score.c)

score.data.1 <- cbind.data.frame(probability.data.study.ref1[,1],ref1.score.c)
colnames(score.data.1) <- c("micro_score","Depression_level")
my_comparisons <- list( c("0", "1"),c("1","2"))
p1 <- ggboxplot(score.data.1, x = "Depression_level", y = "micro_score",
                fill = "Depression_level", palette = c("#caf0f8","#48cae4","#023e8a")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ylab("Micorbiota prediction")+xlab("Depression Level")


fit1t.aov <- aov(micro_score~Depression_level, data=score.data.1)
anova(fit1t.aov)

#ref40 (BDI-II)
probability.data.study.ref40 <- probability.data[selected.data$study == "Ref40",]
ref40.score <- metadata[pmatch(rownames(probability.data.study.ref40),rownames(metadata)),"BDI"]
cor.ref40 <- cor.test(probability.data.study.ref40[,1],ref40.score)


#reference : https://en.wikipedia.org/wiki/Beck_Depression_Inventory
#0–13: minimal depression
#14–19: mild depression
#20–28: moderate depression
#29–63: severe depression
ref40.score.c <- rep(20,length(ref40.score))
ref40.score.c[ref40.score<14] <- 0
ref40.score.c[ref40.score>=14 & ref40.score<=19] <- 1
ref40.score.c[ref40.score>=20 & ref40.score<=28] <- 2
ref40.score.c[ref40.score>=29 & ref40.score<=63] <- 3
cor.test(probability.data.study.ref40[,1],ref40.score.c)

score.data.40 <- cbind.data.frame(probability.data.study.ref40[,1],ref40.score.c)
colnames(score.data.40) <- c("micro_score","Depression_level")
my_comparisons <- list( c("0", "1"),c("1","2"),c("2","3"))
p40 <- ggboxplot(score.data.40, x = "Depression_level", y = "micro_score",
                 fill = "Depression_level", palette = c("#caf0f8","#48cae4","#023e8a","03045e")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ylab("Micorbiota prediction")+xlab("Depression Level")

#ref65 (HADS)
probability.data.study.ref65 <- probability.data[selected.data$study == "Ref65",]
ref65.score <- metadata[pmatch(rownames(probability.data.study.ref65),rownames(metadata)),"Depression_score"]
cor.ref65 <- cor.test(probability.data.study.ref65[,1],ref65.score)
plot(probability.data.study.ref65[,1],ref65.score)

#8 and 10 suggests a mild disorder
#≥11 is considered a clinically significant disorder
#reference : Comparison of two self-rating scales to detect depression: HADS and PHQ-9
ref65.score.c <- rep(0,length(ref65.score))
ref65.score.c[7<ref65.score & ref65.score<11] <- 1
ref65.score.c[ ref65.score >= 11] <- 2
cor.test(probability.data.study.ref65[,1],ref65.score.c)
plot(probability.data.study.ref65[,1],ref65.score.c)

score.data.65 <- cbind.data.frame(probability.data.study.ref65[,1],ref65.score.c)
colnames(score.data.65) <- c("micro_score","Depression_level")
my_comparisons <- list( c("0", "1"),c("1","2"))
p65 <- ggboxplot(score.data.65, x = "Depression_level", y = "micro_score",
                 fill = "Depression_level", palette = c("#caf0f8","#48cae4","#023e8a")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ylab("Micorbiota prediction")+xlab("Depression Level")




#Ref77
probability.data.study.ref77 <- probability.data[selected.data$study == "Ref77",]
ref77.score <- metadata[pmatch(rownames(probability.data.study.ref77),rownames(metadata)),"depression_score"]
cor.ref77 <- cor.test(probability.data.study.ref77[,1],ref77.score)

#reference : Validity, reliability and cut-offs of the Patient Health Questionnaire-9 as a screening tool for depression among patients living with epilepsy in Rwanda
#0–7	Absence of depression
#8–16	Mild
#17–23	Moderate
#≥24	Severe
ref77.score.c <- rep(20,length(ref77.score))
ref77.score.c[ref77.score<=7] <- 0
ref77.score.c[ref77.score>=8 & ref77.score<=16] <- 1
ref77.score.c[ref77.score>=17 & ref77.score<=23] <- 2
ref77.score.c[ref77.score>=24] <- 3
cor.test(probability.data.study.ref77[,1],ref77.score.c)
score.data.77 <- cbind.data.frame(probability.data.study.ref77[,1],ref77.score.c)
colnames(score.data.77) <- c("micro_score","Depression_level")
my_comparisons <- list( c("0", "1"),c("1","2"))
p77 <- ggboxplot(score.data.77, x = "Depression_level", y = "micro_score",
               fill = "Depression_level", palette = c("#caf0f8","#48cae4","#023e8a")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ylab("Micorbiota prediction")+xlab("Depression Level")


####all 4 study
ggarrange(p1,p40,p65,p77,ncol = 2,nrow = 2)

score.data.1 <- cbind.data.frame(probability.data.study.ref1[,1],ref1.score)
score.data.40 <- cbind.data.frame(probability.data.study.ref40[,1],ref40.score)
score.data.65 <- cbind.data.frame(probability.data.study.ref65[,1],ref65.score)
score.data.77<- cbind.data.frame(probability.data.study.ref77[,1],ref77.score)
colnames(score.data.1) <- colnames(score.data.40) <- colnames(score.data.65) <- colnames(score.data.77) <- c("micro_score","Depression")

p1 <- ggplot(score.data.1,aes(x=Depression,y=micro_score))+ geom_point() +theme_classic(base_size=15)+ 
  ylab("Microbiota score")+xlab(paste("Depression score",round(cor.ref1$estimate,4),round(cor.ref1$p.value,4),sep=":"))+ggtitle(id_title[1,2])+
  stat_smooth(method=lm, level=0.95)+theme_classic(base_size=15)

p40 <- ggplot(score.data.40,aes(x=Depression,y=micro_score))+ geom_point() +theme_classic(base_size=15)+ 
  ylab("Microbiota score")+xlab(paste("Depression score",round(cor.ref40$estimate,4),round(cor.ref40$p.value,4),sep=":"))+ggtitle(id_title[4,2])+
  stat_smooth(method=lm, level=0.95)+theme_classic(base_size=15)

p65 <- ggplot(score.data.65,aes(x=Depression,y=micro_score))+ geom_point() +theme_classic(base_size=15)+ 
  ylab("Microbiota score")+xlab(paste("Depression score",round(cor.ref65$estimate,4),round(cor.ref65$p.value,4),sep=":"))+ggtitle(id_title[6,2])+
  stat_smooth(method=lm, level=0.95)+theme_classic(base_size=15)

p77 <- ggplot(score.data.77,aes(x=Depression,y=micro_score))+ geom_point() +theme_classic(base_size=15)+ 
  ylab("Microbiota score")+xlab(paste("Depression score",round(cor.ref77$estimate,4),round(cor.ref77$p.value,4),sep=":"))+ggtitle(id_title[7,2])+
  stat_smooth(method=lm, level=0.95)+theme_classic(base_size=15)

pdf("svm/microbiota_score2depression_score.pdf",width = 10,height = 10)
ggarrange(p1,p40,p65,p77,ncol = 2,nrow = 2,labels="AUTO")
dev.off()
