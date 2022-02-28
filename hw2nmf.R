
## Load data -----

gene <- read.delim("GEUVADIS_normalized_expression_chr20")
geneexp <- gene[,c(1,5:ncol(gene))]
dim(geneexp)
genotype <- read.delim("GEUVADIS_chr20_processed.traw")
library(parallel)
library(doParallel)

## NMF -----
library(NMF)
geneexp_nmf <- t(geneexp[-1])
# transfer into nonnegative data
geneexp_nmf <- nneg(geneexp_nmf)

nmffit <- nmf(geneexp_nmf, rank=4, seed=123456) # higher rank -> better results
w <- nmffit@fit@W
h <- nmffit@fit@H
y_nmf <- w %*% h


## eQTL analysis -----

# regress each gene expression on each genotype
# resdf <- c()
# 
# eqtl <- function(gexp, gene=gene, genotype=genotype){
#   for ( i in 1:nrow(gene)){
#     gtype <- subset(genotype, POS >= (gene$start[i] - 1e6) & POS <= (gene$end[i] + 1e6))
#     gtype <- gtype[,-(1:6)]
#     y <- gexp[,i]
#     
#     resdf1 <- c()
#     for ( j in 1:nrow(gtype)){
#       x <- t(gtype[j,])
#       lmfit <- lm(y ~ x)
#       
#       if (nrow(summary(lmfit)$coefficients)>1) {
#         pval <- summary(lmfit)$coefficients[2,4]
#       } else {
#         pval <- NA
#       }
#       
#       out <- data.frame(gene=gene$gene_id[i],
#                         genotype = genotype$SNP[j],
#                         pval = pval)
#       resdf1 <- rbind(resdf1, out)
#     }
#     resdf1$pval <- p.adjust(resdf1$pval, method = "BH")
#     resdf <- rbind(resdf, resdf1)
#   }
#   
#   return(resdf)
# }
# 
# res_nmf <- eqtl(gexp = y_nmf, gene=gene, genotype=genotype)
# 
# saveRDS(res_nmf, "nmf.rds")

gexp <- y_nmf

cl <- parallel::makeForkCluster(8) 
doParallel::registerDoParallel(cl)

eqtl.analysis <- foreach(i = 1:nrow(gene), .combine = "rbind") %dopar% {
  gtype <- subset(genotype, POS >= (gene$start[i] - 500000) & POS <= (gene$end[i] + 500000))
  gtype <- gtype[,-(1:6)]
  y <- gexp[,i]
  
  all.i.eqtl <- foreach(j = 1:nrow(gtype), .combine = "rbind") %do% {
    x <- t(gtype[j,])
    lmfit <- lm(y ~ x)
    
    if (nrow(summary(lmfit)$coefficients)>1) {
      pval <- summary(lmfit)$coefficients[2,4]
    } else {
      pval <- NA
    }
    
    out <- data.frame(gene=gene$gene_id[i],
                      genotype = genotype$SNP[j],
                      pval = pval)
    return(out)
  }
  
  all.i.eqtl$pval <- p.adjust(all.i.eqtl$pval, method = "BH")
  return(all.i.eqtl)
}

# res_pca <- eqtl(gexp = y_pca, gene=gene, genotype=genotype)
saveRDS(eqtl.analysis, "nmf_dopar.rds")



