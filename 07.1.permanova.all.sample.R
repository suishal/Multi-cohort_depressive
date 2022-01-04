###all sample permanova
dat.cout.cb <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.count.genus.prof",sep="\t",check.names = F)
dat.meta.cb <- read.table("~/Project/Review/meta-analysis/allstudy.1842.sample.metadata.prof",sep="\t",check.names = F)

dim(dat.cout.cb)
dim(dat.meta.cb)
pid <- pmatch(rownames(dat.meta.cb),colnames(dat.cout.cb)) 
dat.cout.cb <- dat.cout.cb[,pid]
dat.cout.cb.m <- data.Normalization(dat.cout.cb,type = "n10",normalization="column")
dat.cout.cb.m.s<-sqrt(dat.cout.cb.m)
colSums(dat.cout.cb.m.s)

bray.dist<-vegdist(t(dat.cout.cb.m.s), method='bray')
set.seed(1) 

bray.div.study<-adonis2(bray.dist~subjectID, data=dat.meta.cb, permutations = 999, method="bray")

dat.meta.cb.age <- dat.meta.cb[!is.na(dat.meta.cb$Age),]
bray.dist.age <-as.dist(as.matrix(bray.dist)[rownames(dat.meta.cb.age),rownames(dat.meta.cb.age)])
bray.div.age<-adonis2(bray.dist.age~Age, data=dat.meta.cb.age, permutations = 999, method="bray")


dat.meta.cb.Sex <- dat.meta.cb[!is.na(dat.meta.cb$Sex),]
bray.dist.Sex <-as.dist(as.matrix(bray.dist)[rownames(dat.meta.cb.Sex),rownames(dat.meta.cb.Sex)])
bray.div.Sex<-adonis2(bray.dist.Sex~Sex, data=dat.meta.cb.Sex, permutations = 999, method="bray")

dat.meta.cb.Depression <- dat.meta.cb[!is.na(dat.meta.cb$Depression),]
dat.meta.cb.Depression$Depression <- as.factor(dat.meta.cb.Depression$Depression)
bray.dist.Depression <-as.dist(as.matrix(bray.dist)[rownames(dat.meta.cb.Depression),rownames(dat.meta.cb.Depression)])
bray.div.Depression<-adonis2(bray.dist.Depression~Depression, data=dat.meta.cb.Depression, permutations = 999, method="bray")

dat.meta.cb.BMI <- dat.meta.cb[!is.na(dat.meta.cb$BMI),]
bray.dist.BMI <-as.dist(as.matrix(bray.dist)[rownames(dat.meta.cb.BMI),rownames(dat.meta.cb.BMI)])
bray.div.BMI<-adonis2(bray.dist.BMI~BMI, data=dat.meta.cb.BMI, permutations = 999, method="bray")

dat.meta.cb.comb <- dat.meta.cb[
  apply(dat.meta.cb[,c("Depression","Age","Sex","BMI","subjectID")],1,function(x){
  sum(is.na(x))==0}),]
bray.dist.comb <-as.dist(as.matrix(bray.dist)[rownames(dat.meta.cb.comb),rownames(dat.meta.cb.comb)])
bray.div.cb <- adonis2(bray.dist.comb~Depression+Age+Sex+BMI+subjectID, data=dat.meta.cb.comb, permutations = 999, method="bray")

permanova.result <- list(bray.div.age,
                         bray.div.Sex,
                         bray.div.BMI,
                         bray.div.study,
                         bray.div.Depression,
                         bray.div.cb
  
)
save(permanova.result,file = "pernomova/result.permanova.rda")
write.table(bray.div.cb,"pernomova/permanova.1634sample.txt",sep="\t",quote = F)
