# Loading libraries
library(DESeq2)
library(vsn)
library(ggplot2)
library(limma)
library(ashr)


# Functions ####
plotPCAWithSampleNames = function(x, targets=targets, intgroup=colnames(targets)[1], ntop=500)
{
  library(RColorBrewer)
  library(genefilter)
  library(lattice)
  
  # pca
  #rv = rowVars(assay(x))
  rv = rowVars(x)
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  #pca = prcomp(t(assay(x)[select,]))
  pca = prcomp(t(x[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  names = colnames(x)
  #names = as.character(x$sample)
  
  # factor of groups
  #fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  fac = factor(apply(as.data.frame(targets[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  # colors
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors = brewer.pal(nlevels(fac), "Set1")
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # plot
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}

# LOADING FILES ####
countTable.g = read.csv("gene_count_matrix.csv",row.names=1)
countTable.t = read.csv("transcript_count_matrix.csv",row.names=1)
colTable = read.csv("design.csv",row.names=1)

#Sorting count tables by the col table rownames
countTable.g = countTable.g[,rownames(colTable)]
countTable.t = countTable.t[,rownames(colTable)]


# PROCESSING ####
# Creating dds dataset
dds.g = DESeqDataSetFromMatrix(countData=countTable.g, colData=colTable, design= ~ group + replicate)
dds.t = DESeqDataSetFromMatrix(countData=countTable.t, colData=colTable, design= ~ group + replicate)
#dds without batch (replicate) in the design for batch removal)
dds.g.group = DESeqDataSetFromMatrix(countData=countTable.g, colData=colTable, design= ~ group)

#Filtering genes/transcripts to have a minimum count of 10 in minimum 3 samples
minimum_samples = 3
filter.g = rowSums(counts(dds.g) >= 10) >= minimum_samples
dds.g = dds.g[filter.g,]

filter.g.group = rowSums(counts(dds.g.group) >= 10) >= minimum_samples
dds.g.group = dds.g.group[filter.g.group,]

filter.t = rowSums(counts(dds.t) >= 10) >= minimum_samples
dds.t = dds.t[filter.t,]

#Creating dds object
dds.g = DESeq(dds.g)
dds.g.group = DESeq(dds.g.group)
dds.t = DESeq(dds.t)

# MAKE PLOTS ####
#Dispersion plots
svg("Gene_dispersion.svg")
plotDispEsts(dds.g, main="Gene level")
dev.off()

svg("Transcript_dispersion.svg")
disp.t = plotDispEsts(dds.t, main="Transcript level")
dev.off()

#gene level rlog transformation
rld.g = rlog(dds.g,blind=TRUE)
rld_assay.g = assay(rld.g)

rld.g.group = rlog(dds.g.group,blind=TRUE)
rld_assay.g.group = assay(rld.g.group)

#PCA plot, several variants to assess which is more informative
svg("Gene_pca_group.svg")
plotPCAWithSampleNames(rld_assay.g, targets=colTable, intgroup=c("group"))
dev.off()
svg("Gene_pca_batch.svg")
plotPCAWithSampleNames(rld_assay.g, targets=colTable, intgroup=c("replicate"))
dev.off()
svg("Gene_pca_gb.svg")
plotPCAWithSampleNames(rld_assay.g, targets=colTable, intgroup=c("group", "replicate"))
dev.off()

#transcript level rlog transformation
rld.t = rlog(dds.t,blind=TRUE)
rld_assay.t = assay(rld.t)

#PCA plot
svg("Transcript_pca.svg")
plotPCAWithSampleNames(rld_assay.t, targets=colTable, intgroup=c("group"))
dev.off()


# GENE COUNT ONLY SECTION ####
#log2 transformation
lgc.norm = log2(counts(dds.g,normalized=TRUE)+0.0001)
#log2 SD vs mean plot
norm_msd=meanSdPlot(lgc.norm)
#Adding title to ggplot
norm_msd$gg + ggtitle("log2(normalised)")
ggsave("meanSDplot_log2.svg", width=5, height=4)
#rlog SD vs mean plot
rld_msd=meanSdPlot(rld_assay.g)
rld_msd$gg + ggtitle("rlog")
ggsave("meanSDplot_rlog.svg", width=5, height=4)

#Unshrunk results
res1 = results(dds.g, contrast=c("group","B","A"), lfcThreshold=0)
summary(res1,alpha=0.05)
#post hoc subset for comparison
res1.ph = subset(res1,abs(res1$log2FoldChange) > 1)
summary(res1.ph,alpha=0.05)

res2 = results(dds.g, contrast=c("group","B","A"), lfcThreshold=1)
summary(res2,alpha=0.05)

svg("MA_B_lfc0.svg")
DESeq2::plotMA(res1,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: AvsB, lfc=0")
dev.off()
svg("MA_B_lfc1.svg")
DESeq2::plotMA(res2,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: AvsB, lfc<1")
abline(h=c(-1,1), col="black")
dev.off()

res3 = results(dds.g, contrast=c("group","C","A"), lfcThreshold=0)
summary(res3,alpha=0.05)
res3.ph = subset(res3,abs(res3$log2FoldChange) >1)
summary(res3.ph,alpha=0.05)
res4 = results(dds.g, contrast=c("group","C","A"), lfcThreshold=1)
summary(res4,alpha=0.05)

svg("MA_C_lfc0.svg")
DESeq2::plotMA(res3,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: AvsC, lfc=0")
dev.off()
svg("MA_C_lfc1.svg")
DESeq2::plotMA(res4,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: AvsC, lfc<1")
abline(h=c(-1,1), col="black")
dev.off()

#Shrunken results
res1.2 = lfcShrink(dds.g, contrast=c("group","B","A"),type="ashr")
summary(res1.2,alpha=0.05)
res2.2 = lfcShrink(dds.g, contrast=c("group","B","A"),type="ashr", lfcThreshold=1)
summary(res2.2,alpha=0.05)

svg("MA_B_lfc0_shrink.svg")
DESeq2::plotMA(res1.2,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: shrunken AvsB, lfc=0")
dev.off()
svg("MA_B_lfc1_shrink.svg")
DESeq2::plotMA(res2.2,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: shrunken AvsB, lfc=1")
abline(h=c(-1,1), col="black")
dev.off()

res3.2 = lfcShrink(dds.g, contrast=c("group","C","A"),type="ashr")
summary(res3.2,alpha=0.05)
res4.2 = lfcShrink(dds.g, contrast=c("group","C","A"),type="ashr", lfcThreshold=1)
summary(res4.2,alpha=0.05)

svg("MA_C_lfc0_shrink.svg")
DESeq2::plotMA(res3.2,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: shrunken AvsC, lfc=0")
dev.off()
svg("MA_C_lfc1_shrink.svg")
DESeq2::plotMA(res4.2,alpha=0.05, cex.lab = 1.3, main="Gene level MvA: shrunken AvsC, lfc=1")
abline(h=c(-1,1), col="black")
dev.off()

#Filter, sort and write to file
res2.2.s = res2.2[order(res2.2$log2FoldChange),decreasing=FALSE]
res2.2_sig = subset(res2.2.s,res2.2.s$padj < 0.05)
write.csv(res2.2_sig,file="Significant.Genes.BvsA.csv",quote=FALSE)

res4.2.s = res4.2[order(res4.2$log2FoldChange),decreasing=FALSE]
res4.2_sig = subset(res4.2.s,res4.2.s$padj < 0.05)
write.csv(res4.2_sig,file="Significant.Genes.CvsA.csv",quote=FALSE)

# BATCH REMOVAL #####
mydesign = model.matrix(design(dds.g.group),colData(dds.g.group))
b.corrected = limma::removeBatchEffect(rld_assay.g.group, batch=colData(dds.g)$replicate, design=mydesign)
#Batch corrected PCA
svg("Gene_pca_batchcorrected.svg")
plotPCAWithSampleNames(b.corrected,targets=colTable, intgroup='group')
dev.off()
write.csv(b.corrected,file="BatchCorrected.Rlog.csv",quote=FALSE)
