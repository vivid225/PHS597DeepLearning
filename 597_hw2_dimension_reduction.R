# In this HW, we will explore how these algorithms can preserve information for eQTL analysis. 
# eQTL analysis, in simple terms, is to regress gene expression level over genotypes. 
# For HW, please apply VQ, PCA and NMF to reduce the dimension of gene expression levels to a given 
# dimension (which will vary), and examine how many eQTLs can be recovered from compressed data as compared 
# to the analysis of uncompressed gene expression level. 

## Load data -----

gene <- read.delim("/Users/rrrrrita/Documents/GitHub/PHS597DeepLearning/gene_expression_sample/GEUVADIS_normalized_expression_chr20")
geneexp <- gene[,c(1,5:ncol(gene))]
dim(geneexp)
genotype <- read.delim("/Users/rrrrrita/Documents/GitHub/PHS597DeepLearning/gene_expression_sample/GEUVADIS_chr20_processed.traw")


## Vector quantization (VQ) -----
library(factoextra)
library(purrr)
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(t(geneexp[,-1]), k, nstart = 10 )$tot.withinss
}
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- seq(200,300, by =10)
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")

fviz_nbclust(t(geneexp[-1]), kmeans, method = "wss")
# 3 clusters are sufficient by elbow method

kres <- kmeans(t(geneexp[-1]), 3, nstart = 1)
fviz_cluster(kres, data = t(geneexp[-1]))
# the data itself is not very easy to be dicerned


y_vq <- kres$centers[kres$cluster,]
rownames(y_vq) <- colnames(geneexp[,-1])


## PCA -----
### pca on gene expression:
pca <- prcomp(t(geneexp[-1]), center=TRUE, scale = TRUE)
pc <- pca$x
rotation <- pca$rotation
plot(pca, npcs=10)
# first 10 pcs explain >90% of the variances 

npc <- 10
colm <- colMeans(t(geneexp[,-1]))

y_pca <- pc[,1:npc] %*% t(rotation[,1:npc]) + colm
# View(y_pca)


## NMF -----
library(NMF)
geneexp_nmf <- t(geneexp[-1])
# transfer into nonnegative data
geneexp_nmf <- nneg(geneexp_nmf)

nmffit <- nmf(geneexp_nmf, rank=10, "brunet", seed=123456) # higher rank -> better results
w <- nmffit@fit@W
h <- nmffit@fit@H
y_nmf <- w %*% h

# estim.r <- nmf(geneexp_nmf, seq(30,40, by=10), nrun=2, seed=123456)
# plot(estim.r)
# estim.r

## eQTL analysis -----

# regress each gene expression on each genotype
resdf <- c()

eqtl <- function(gexp){
  nrow(gene);nrow(gtype)
  for ( i in 1:100){
    gtype <- subset(genotype, POS >= (gene$start[i] - 1e6) & POS <= (gene$end[i] + 1e6))
    gtype <- gtype[,-(1:6)]
    y <- gexp[,i]
    
    resdf1 <- c()
    for ( j in 1:10){
      x <- t(gtype[j,])
      lmfit <- lm(y ~ x)
      
      out <- data.frame(gene=gene$gene_id[i],
                        genotype = genotype$SNP[j],
                        pval = summary(lmfit)$coefficients[2,4])
      resdf1 <- rbind(resdf1, out)
    }
    resdf1$pval <- p.adjust(resdf1$pval, method = "BH")
    resdf <- rbind(resdf, resdf1)
  }
  
  return(resdf)
}

res_pca <- eqtl(gexp = y_pca)





