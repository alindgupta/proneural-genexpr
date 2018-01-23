library(affy)
library(limma)
library(drosophila2.db)
library(topGO)
library(apcluster)
library(magrittr)
library(Rtsne)
library(dbscan)
library(RColorBrewer)

#########################################################################################
                                        # process celfiles
#########################################################################################
setwd("~/Documents/drosophila/data/jarman/")
affyDirs <- c(
    "atoGFP+/t1",
    "atoGFP-/t1",
    "mutGFP+", # dont care about mutant
    "mutGFP-", # dont care about mutant
    "atoGFP+/t2",
    "atoGFP-/t2",
    "atoGFP+/t3",
    "atoGFP-/t3")

affyBatch <- sapply(
    affyDirs,
    function(x) ReadAffy(celfile.path=paste("celfiles/", x, sep="")))
eset <- lapply(affyBatch, rma)
e <- lapply(eset, exprs)
rm(affyBatch, eset)


#########################################################################################
                                        # fit a linear model and perform empirical Bayes shrinkage
#########################################################################################
linearFit <- function(gfp, nogfp) {
    design <- model.matrix(~ 0
                           + c(rep(1, ncol(gfp)),
                               rep(0, ncol(nogfp)))
                           + c(rep(0, ncol(gfp)),
                               rep(1, ncol(nogfp))))
    colnames(design) <- c("gfp", "nogfp")
    glm <- lmFit(cbind(gfp, nogfp), design=design) %>% eBayes
    contrasts.matrix <- makeContrasts(dif=gfp-nogfp, levels=design)
    fit <- contrasts.fit(glm, contrasts.matrix) %>% eBayes
    return(list(glm$coefficients, fit$coefficients))
}

# execute linear fits on celfile subsets
# don't care about e[[3]], e[[4]] which is the ato mutant
out.a <- linearFit(e[[1]], e[[2]])
out.b <- linearFit(e[[5]], e[[6]])
out.c <- linearFit(e[[7]], e[[8]])
glm <- cbind(out.a[[1]], out.b[[1]], out.c[[1]])
fit <- cbind(out.a[[2]], out.b[[2]], out.c[[2]])
rm(out.a, out.b, out.c)

# replace probe id with gene symbols and remove probs that are NA
# map from id to gene symbol
tab <- AnnotationDbi::select(
                          drosophila2.db,
                          keys=keys(drosophila2.db),
                          columns=c("SYMBOL"))
row.names(glm) <- tab$SYMBOL[match(row.names(glm), tab$PROBEID)]
glm <- glm[!is.na(row.names(glm)),]
row.names(fit) <- tab$SYMBOL[match(row.names(fit), tab$PROBEID)]
fit <- fit[!is.na(row.names(fit)),]


#########################################################################################
                                        # reduce dataset size
#########################################################################################
threshold.abs <- 0.5 # threshold for absolute expression level
threshold.var <- 0.6 # threshold for variance over the three time points

temp <- cbind(glm[,1], glm[,3], glm[,5])
cutoff <- quantile(temp, threshold.abs)[[1]]
mask <- apply(temp, 1, function(x) any(x[1] > cutoff,
                                       x[2] > cutoff,
                                       x[3] > cutoff))
glm.trunc <- glm[mask,]
fit.trunc <- fit[mask,]

fit.variances <- apply(fit.trunc, 1, var)
cutoff <- quantile(fit.variances, threshold.var)[[1]]
mask <- fit.variances > cutoff
glm.trunc <- glm.trunc[mask,]
fit.trunc <- fit.trunc[mask,]

# check for presence of genes of interest
"sktl" %in% row.names(glm.trunc)
"INPP5E" %in% row.names(glm.trunc)
"unc" %in% row.names(glm.trunc)
"nompB" %in% row.names(glm.trunc)
"PI4KIIIalpha" %in% row.names(fit.trunc)
nrow(fit.trunc)

#########################################################################################
                                # run affinity propagation
#########################################################################################
if (nrow(fit.trunc) < 5000) {
    ap.out <- apcluster(negDistMat(r=2), fit.trunc)
} else {
    print("Warning: trying to run apcluster with >4000 datapoints")
}


#########################################################################################
                                # GO annotations and enrichments
#########################################################################################
clustNum <- 7 # number of cluster containing sktl
cluster <- ap.out[[clustNum]] %>% names

# the parameter intrs is a vector of gene symbols
analyzeGO <- function(intrs) {
    alls <- row.names(fit.trunc)
    geneList <- as.integer(alls %in% intrs) %>% factor
    names(geneList) <- alls
    GOdata <- new(
        "topGOdata",
        ontology="BP",
        allGenes=geneList,
        geneSel=function(x) x < 0.01,
        description="Test",
        annot=annFUN.org,
        mapping="org.Dm.eg.db",
        ID="SYMBOL")
    resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
    GenTable(GOdata, classicFisher=resultFisher, topNodes=10)
}
matplot(t(fit.trunc[row.names(fit.trunc) %in% names(ap.out[[clustNum]]),]), type="l")


#########################################################################################
                                        # run t-SNE
#########################################################################################
tsne.out <- Rtsne(fit.trunc, theta=0.0, verbose=TRUE) %>% (function(x) x$Y)
row.names(tsne.out) <- row.names(fit.trunc)
mask <- tsne.out[row.names(tsne.out) %in% names(ap.out[[clustNum]]),]

# run hdbscan
geneName <- "sktl"
minPts <- 5 # want this to be argmax<n>(num_clusters > arbitrary_num_clusters)
gene.at <- which(row.names(fit.trunc) == geneName)
hdbscan.out <- hdbscan(cbind(as.vector(tsne.out[,1]),
                             as.vector(tsne.out[,2])),
                       minPts=minPts)
clusterAssignment <- hdbscan.out$cluster[gene.at]
cluster.genes <- which(hdbscan.out$cluster == clusterAssignment)
analyzeGO(row.names(tsne.out)[cluster.genes])
plot(tsne.out, col=hdbscan.out$cluster)


#########################################################################################
                                        # plotting tsne embeddings
#########################################################################################
col1 <- brewer.pal(9, "Greys")
set1 <- c(brewer.pal(8, "Set1"))

par(bty="n",
    mfrow=c(1,1),
    mar=c(3, 3, 3, 3),   # plot margins b-l-t-r
    las=1,               # horizontal labels
    tcl=-.25,            # tick length
    font.main=1,         # plain font
    mgp=c(2.5, 0.5, 0),  # axis spacings
    cex.main=1.2)

plot(tsne.out,
     col=col1[4],
     pch=19,
     lwd=0.5, 
     xlab="", 
     ylab="", 
     xaxt="n", 
     yaxt="n")

# genes of interest
sktl.at = which(row.names(fit.trunc) == "sktl")
inpp5e.at = which(row.names(fit.trunc) == "INPP5E")

box(bty="L")
points(mask[,1], mask[,2], col=col1[7], pch=19)
points(tsne.out[sktl.at,1], tsne.out[sktl.at,2], col=set1[1], pch=19)
points(tsne.out[inpp5e.at,1], tsne.out[inpp5e.at,2], col=set1[3], pch=19)
title("t-SNE plot of gene expression profiles", font.main=1, line=1)
title(ylab="t-SNE dimension 2", line=1)
title(xlab="t-SNE dimension 1", line=1)

legend("topright", 
       legend=c("All genes",
                as.expression(bquote("AP cluster for"~italic("sktl"))),
                as.expression(bquote(italic("sktl"))),
                as.expression(bquote(italic("inpp5e")))),
       col=c(col1[4],
             col1[7],
             set1[1],
             set1[3]),
       border=NA,
       bty="n",
       pch=19)

