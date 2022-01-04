
load("~/Project/Review/meta-analysis/pernomova/ref1.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref23.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref27.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref40.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref60.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref65.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref77.output.RData")
load("~/Project/Review/meta-analysis/pernomova/ref79.output.RData")

adonis.list.file <- ls(pattern="list.adonis.ref")
adonis.estimat.bray <- c()
adonis.estimat.jsd <- c()
adonis.estimat.unw <- c()
adonis.estimat.wuni <- c()
for(i in adonis.list.file){
  txt <- paste("estimat.list <-", i,sep="")
  eval((parse(text =txt)))
  
  adonis.estimat.bray <- rbind(adonis.estimat.bray,
                               c(estimat.list$adonis.output[1,3],
                                 estimat.list$adonis.R.CI.bootrap[1,1],
                                 estimat.list$adonis.R.CI.bootrap[1,2],
                                 estimat.list$adonis.R.CI.bootrap[1,3],
                                 estimat.list$adonis.output[1,5]))
  adonis.estimat.jsd <- rbind(adonis.estimat.jsd,
                              c(estimat.list$adonis.output[2,3],
                                estimat.list$adonis.R.CI.bootrap[2,1],
                                estimat.list$adonis.R.CI.bootrap[2,2],
                                estimat.list$adonis.R.CI.bootrap[2,3],
                                estimat.list$adonis.output[2,5]))
  adonis.estimat.unw <- rbind(adonis.estimat.unw,
                              c(estimat.list$adonis.output[3,3],
                                estimat.list$adonis.R.CI.bootrap[3,1],
                                estimat.list$adonis.R.CI.bootrap[3,2],
                                estimat.list$adonis.R.CI.bootrap[3,3],
                                estimat.list$adonis.output[3,5]))
  adonis.estimat.wuni <- rbind(adonis.estimat.wuni,
                               c(estimat.list$adonis.output[4,3],
                                 estimat.list$adonis.R.CI.bootrap[4,1],
                                 estimat.list$adonis.R.CI.bootrap[4,2],
                                 estimat.list$adonis.R.CI.bootrap[4,3],
                                 estimat.list$adonis.output[4,5]))
  
}
colnames(adonis.estimat.bray) <-colnames(adonis.estimat.jsd) <-colnames(adonis.estimat.unw) <-colnames(adonis.estimat.wuni) <-c("R2","boostrap_r2","boostrap_2.5%","boostrap_97.5%","pvalue")
id_title <- read.table("~/Project/Review/meta-analysis/ID2Title_meta_analysis_study.txt",sep="\t",head=T)

dat.R2 <- cbind(adonis.estimat.bray[,1],adonis.estimat.jsd[,1],
                adonis.estimat.unw[,1],adonis.estimat.wuni[,1] )
dat.p <- cbind(adonis.estimat.bray[,5],adonis.estimat.jsd[,5],
               adonis.estimat.unw[,5],adonis.estimat.wuni[,5] )
colnames(dat.R2) <- colnames(dat.p) <- c("Bray–Curtis dissimilarity","Jensen–Shannon divergence","Unweighted UniFrac distance","Weighted UniFrac distance")
rownames(dat.R2) <- rownames(dat.p) <- id_title$title
dat.p.sig <-dat.p
dat.p.sig[dat.p<0.05] <- "*"
dat.p.sig[dat.p>0.05] <- ""
pdf("permanova.heatmap.pdf")
pheatmap::pheatmap(mat = dat.R2,display_numbers = dat.p.sig)
dev.off()
#Beta-Dispersion heatmap

betadisper.estimat.bray <- c()
betadisper.estimat.jsd <- c()
betadisper.estimat.unw <- c()
betadisper.estimat.wuni <- c()
for(i in adonis.list.file){
  txt <- paste("estimat.list <-", i,sep="")
  eval((parse(text =txt)))
  
  betadisper.estimat.bray <- rbind(betadisper.estimat.bray,
                                   c(estimat.list$betadisper.HSD[1,1],
                                     estimat.list$betadisper.HSD[1,2],
                                     estimat.list$betadisper.HSD[1,3],
                                     estimat.list$betadisper.HSD[1,4]))
  betadisper.estimat.jsd <- rbind(betadisper.estimat.jsd,
                                  c(estimat.list$betadisper.HSD[2,1],
                                    estimat.list$betadisper.HSD[2,2],
                                    estimat.list$betadisper.HSD[2,3],
                                    estimat.list$betadisper.HSD[2,4]))
  betadisper.estimat.unw <- rbind(betadisper.estimat.unw,
                                  c(estimat.list$betadisper.HSD[3,1],
                                    estimat.list$betadisper.HSD[3,2],
                                    estimat.list$betadisper.HSD[3,3],
                                    estimat.list$betadisper.HSD[3,4]))
  betadisper.estimat.wuni <- rbind(betadisper.estimat.wuni,
                                   c(estimat.list$betadisper.HSD[4,1],
                                     estimat.list$betadisper.HSD[4,2],
                                     estimat.list$betadisper.HSD[4,3],
                                     estimat.list$betadisper.HSD[4,4]))
  
}
colnames(betadisper.estimat.bray) <-colnames(betadisper.estimat.jsd) <-colnames(betadisper.estimat.unw) <-colnames(betadisper.estimat.wuni) <-c("diff1-0","lwr","upr","pvalue")

##heatmap

dat.diff <- cbind(betadisper.estimat.bray[,1],betadisper.estimat.jsd[,1],
                  betadisper.estimat.unw[,1],betadisper.estimat.wuni[,1] )
dat.diff.p <- cbind(betadisper.estimat.bray[,4],betadisper.estimat.jsd[,4],
                    betadisper.estimat.unw[,4],betadisper.estimat.wuni[,4] )
colnames(dat.diff) <- colnames(dat.p) <- c("Bray–Curtis dissimilarity","Jensen–Shannon divergence","Unweighted UniFrac distance","Weighted UniFrac distance")
rownames(dat.diff) <- rownames(dat.p) <- id_title$title

dat.diff.p.sig <-dat.diff.p
dat.diff.p.sig[dat.diff.p<0.05] <- "*"
dat.diff.p.sig[dat.diff.p>0.05] <- ""
pdf("betadisper.heatmap.pdf")
pheatmap::pheatmap(mat = dat.diff,display_numbers = dat.diff.p.sig)
dev.off()


library(ComplexHeatmap)
dat.diff <- round(dat.diff,2)
col_fun = colorRamp2(c(-0.07, 0, 0.06), c( "#219ebc", "#023047","#fb8500"))
ht1<-Heatmap(dat.diff,name="betadisper",cluster_rows=T,cluster_columns=T,col = col_fun,
             layer_fun = function(j, i, x, y, width, height, fill) {
               v = pindex(dat.diff.p.sig, i, j)
               grid.text(v, x, y, gp = gpar(fontsize = 10))
             })
library(circlize)
col_fun = colorRamp2(c(0, 0.04, 0.1), c("#caf0f8", "#00b4d8", "#03045e"))
ht2<-Heatmap(dat.R2,name="R2",cluster_rows=T,cluster_columns=T,col = col_fun,
             layer_fun = function(j, i, x, y, width, height, fill) {
               v = pindex(dat.p.sig, i, j)
               grid.text(v, x, y, gp = gpar(fontsize = 10))
             })

ht_list =  ht2+ht1
pdf("~/Project/Review/meta-analysis/pernomova/permanova_R_betasper_diff.heatmap.pdf")
draw(ht_list)
dev.off()
