gene <- read.delim("GEUVADIS_normalized_expression_chr20")
geneexp <- gene[,c(1,5:ncol(gene))]
dim(geneexp)
genotype <- read.delim("GEUVADIS_chr20_processed.traw")

library(NMF)
geneexp_nmf <- t(geneexp[-1])
geneexp_nmf <- nneg(geneexp_nmf)

# nmffit <- nmf(geneexp_nmf, rank=50, "brunet", seed=123456) # higher rank -> better results
# w <- nmffit@fit@W
# h <- nmffit@fit@H
# x_nmf <- w %*% h

estim.r <- nmf(geneexp_nmf, 1:10, nrun=20, seed=123456)
saveRDS(estim.r, "estim.rds")