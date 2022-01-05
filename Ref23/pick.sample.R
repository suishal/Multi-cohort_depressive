dat <- read.table("feces.full.record.list.pick",sep="\t",head=T)
na.sum <- apply(dat,1,function(x){sum(x=="Not provided")})
summary(na.sum)
dat.p <- dat[na.sum==0,]

dat.p <- dat.p[dat.p$env_package == "human-gut", ]
pick <- (dat.p$mental_illness_type_depression == "Yes"| dat.p$mental_illness_type_depression == "No") & (dat.p$subset_healthy == "True" | dat.p$subset_healthy == "TRUE" ) & (dat.p$sex == "female" | dat.p$sex == "male" )
dat.p <-dat.p[pick,]
write.table(dat.p,"feces.full.record.list.pick.pick.txt",sep="\t",quote=F,row.names=F)

dat <- read.table("feces.full.record.list.pick.pick.txt",sep="\t",head=T,stringsAsFactors=F)
dat$sample_name <-as.character( sprintf("%0.9f", dat$sample_name ))
rownames(dat) <- dat$sample_name
dat$mental_illness_type_depression_num <- 0
dat$mental_illness_type_depression_num[ dat$mental_illness_type_depression=="Yes"] <- 1

ebi.list <- read.table("study.txt",sep="\t",check.names=F,stringsAsFactors=F,colClasses=c("character","character"),head=T)
sp.str <- strsplit(as.character(ebi.list$description),":")
ebi.list$V3 <-  sapply(sp.str,function(x){x[2]})
ebi.list$V4 <-  sapply(sp.str,function(x){x[1]})

miss.id <- dat$sample_name[is.na(match(dat$sample_name,ebi.list$V3))]
dat <- dat[!dat$sample_name %in% miss.id,]


pe.sample <- read.table("pair_end.sample")
dat <- dat[!dat$sample_name %in% pe.sample$V1,]

read.count <- read.table("filereport_read_run_PRJEB11419_tsv.readcount.txt",sep="\t",head=T)
enough.read <- read.count[read.count$read_count > 10000 & !is.na(read.count$read_count ),]
enough.read.id.sep <-  strsplit(as.character(enough.read$sample_alias),":")
enough.read$sampleid <-  sapply(enough.read.id.sep,function(x){x[2]})

dat <- dat[dat$sample_name %in% enough.read$sampleid,]

m_ps <- glm(mental_illness_type_depression_num ~ age_years + bmi_corrected + sex,family = binomial(), data = dat)
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),depression = m_ps$model$mental_illness_type_depression_num)
head(prs_df)

library("MatchIt")
mod_match <- matchit(mental_illness_type_depression_num ~ age_years + bmi_corrected + sex,method = "nearest", data = dat)
dta_m <- match.data(mod_match)
#write.table(dta_m,"match_sample_metadata.txt",sep="\t",quote=F)

for(i in 1:dim(dta_m)[1]){
    sample <- dta_m$sample_name[i]
    ebi.path <- ebi.list[ebi.list$V3 == sample,]
    if(dim(ebi.path)[1]>1){
	x <- ebi.path[,1]
	dta_m$ebi[i] <- paste(x,collapse=";")
     }
     else{
	 dta_m$ebi[i] <- ebi.path[,1]
	 }
}

write.table(dta_m,"match_sample_metadata.txt",sep="\t",quote=F)
