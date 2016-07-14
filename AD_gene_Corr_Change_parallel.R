rm(list=ls())
library(parallel)
library(cocor)
library(gtools)
options(stringsAsFactors = F)

# ==== Expression Background ====
load("/shared/hidelab2/user/md1wwxx/ADCorrChange/data/AMP_AD_ROSMAP_RNAseq_normalized1.RData")
colnames(expressionRank)
head(rownames(expressionRank))
iqr<-apply(expressionRank, 1, IQR)
qt<-quantile(iqr, probs=.995)
exprs_rnk<-expressionRank[iqr>qt,]
combination_of_genes<-combinations(nrow(exprs_rnk),2)
number_of_combinations<-nrow(combination_of_genes)

idxAD<-grep("AD", colnames(exprs_rnk))
t_exprs_rnk = exprs_rnk[,idxAD]
c_exprs_rnk = exprs_rnk[,-idxAD]
n.c<-ncol(c_exprs_rnk)
n.t<-ncol(t_exprs_rnk)
gene.names<-rownames(exprs_rnk)

ProcessElement <- function(ic){
  	A<-combination_of_genes[ic,1]
  	B<-combination_of_genes[ic,2]
    
  	c_A = c_exprs_rnk[A,]
  	c_B = c_exprs_rnk[B,]
  		
  	t_A = t_exprs_rnk[A,]
  	t_B = t_exprs_rnk[B,]
  	
  	# get correlation between the summaries for the unique genes
  	tmp = data.frame(Gene.A=gene.names[A],Gene.B=gene.names[B])
  	c_cortest<-cor.test(c_A, c_B, method="spearman")
  	t_cortest<-cor.test(t_A, t_B, method="spearman")
  	rc<-c_cortest$estimate
  	rt<-t_cortest$estimate
  	diffcor<-cocor.indep.groups(rc, rt, n.c, n.t)
  	tmp$r.c<-rc
  	tmp$p.c<-c_cortest$p.value
  	tmp$n.c<-n.c
  	tmp$r.t<-rt
  	tmp$p.t<-t_cortest$p.value
  	tmp$n.t<-n.t
  	tmp$p.cocor<-diffcor@fisher1925$p.value

  	setTxtProgressBar(pb,ic)
  	return(tmp)
}

nc = detectCores()
# loop through Gene combinations
input = 1:number_of_combinations
pb = txtProgressBar(min=0,max=number_of_combinations,style=3,initial=0)
cat("\n")
res = mclapply(input,ProcessElement,mc.cores=nc)
close(pb)

# save results
saveRDS(res,paste0("/shared/hidelab2/user/md1wwxx/ADCorrChange/output/hallmark/ROSMAP_AD_NCI__gene_corr_change.RDS"))
result<-data.frame()
for(i in 1:length(res)){
 	result<-rbind(result, res[[i]])
}
result$FDR<-p.adjust(result$p.cocor, method="fdr")
write.table(result, file="/shared/hidelab2/user/md1wwxx/ADCorrChange/output/hallmark/ROSMAP_AD_NCI_gene_corr_change.txt", sep="\t", row.names=FALSE, quote = FALSE)

rm(list=ls())
