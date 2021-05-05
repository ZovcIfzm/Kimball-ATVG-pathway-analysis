### Installing packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("statmod")

# Fix for org.Mm.eg.db package failing
remotes::install_version("RSQLite", version = "2.2.5")
options(connectionObserver = NULL)

# Packages to install
BiocManager::install("Rsubread")
BiocManager::install("edgeR")
BiocManager::install("RnaSeqGeneEdgeRQL")
BiocManager::install("org.Mm.eg.db")

### Setting up the groups variable
# library(RnaSeqGeneEdgeRQL)

#targetsFile <- system.file("extdata", "targets.txt", package="RnaSeqGeneEdgeRQL")
#targets <- read.delim(targetsFile, stringsAsFactors=FALSE)
targets <- read.delim("../results/groups.txt", stringsAsFactors=FALSE)
targets

group <- paste(targets$CellType, targets$Status, sep=".")
group <- group[c(1,3,2,4,5,7,6,8)]
group <- factor(group)
table(group)
group

### Preliminary analysis
'''
FileURL <- paste(
  "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
  "format=file",
  "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
  sep="&")
download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")
'''
#GenewiseCounts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz",
#                             row.names="EntrezGeneID")
GenewiseCounts <- read.delim("../results/read_counts.txt", row.names="Geneid")[,-1:-4]
GenewiseCounts <- GenewiseCounts[,c(1,2,4,3,5,6,8,7,9)]
#colnames(GenewiseCounts) <- colnames(GenewiseCounts)
dim(GenewiseCounts)
head(GenewiseCounts)

library(edgeR)
y <- DGEList(GenewiseCounts[,-1], group=group,
             genes=GenewiseCounts[,1,drop=FALSE])
options(digits=3)
y$samples

### Adding gene annotation
library(org.Mm.eg.db)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENTREZID", column="SYMBOL")
head(y$genes)

# Drop rows where there is no gene symbol
y <- y[!is.na(y$genes$Symbol), ]
dim(y)

### Filtering to remove low counts
keep <- rowSums(cpm(y) > 0.5) >= 2
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

### Normalization for composition bias
y <- calcNormFactors(y)
y$samples

### Exploring differences between libraries
pch <- c(0,1,15,16)
colors <- rep(c("darkgreen","blue"), 2)
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=2)
abline(h=0, col="red", lty=2, lwd=2)

plotMD(y, column=6)
abline(h=0, col="red", lty=2, lwd=2)

### Design Matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

### Dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

### CURRENT WORKING POINT

### Differential expression analysis
### Testing for differential expression
B.LvsP <- makeContrasts(DIO.Low-ND.Low, levels=design)
res <- glmQLFTest(fit, contrast=B.LvsP)
topTags(res)

is.de <- decideTestsDGE(res)
summary(is.de)

plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

### Differential expression above a fold-change threshold
tr <- glmTreat(fit, contrast=B.LvsP, lfc=log2(1.5))
topTags(tr)

is.de <- decideTestsDGE(tr)
summary(is.de)

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

### Heatmap clustering
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")

o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:30],]
logCPM <- t(scale(t(logCPM)))

library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

### Analysis of deviance
con <- makeContrasts(
  L.PvsL = L.pregnant - L.lactating,
  L.VvsL = L.virgin - L.lactating,
  L.VvsP = L.virgin - L.pregnant, levels=design)

res <- glmQLFTest(fit, contrast=con)
topTags(res)

### Complicated contrasts
con <- makeContrasts(
  (L.lactating-L.pregnant)-(B.lactating-B.pregnant), 
  levels=design)
res <- glmQLFTest(fit, contrast=con)
topTags(res)

### Pathway analysis
### Gene ontology analysis
go <- goana(tr, species="Mm")
topGO(go, n=15)

### KEGG pathway analysis
keg <- kegga(tr, species="Mm")
topKEGG(keg, n=15, truncate=34)

### FRY gene set tests
library(GO.db)
cyt.go <- c("GO:0032465", "GO:0000281", "GO:0000920")
term <- select(GO.db, keys=cyt.go, columns="TERM")
term

Rkeys(org.Mm.egGO2ALLEGS) <- cyt.go
cyt.go.genes <- as.list(org.Mm.egGO2ALLEGS)

B.VvsL <- makeContrasts(B.virgin-B.lactating, levels=design)
fry(y, index=cyt.go.genes, design=design, contrast=B.VvsL)

res <- glmQLFTest(fit, contrast=B.VvsL)
index <- rownames(fit) %in% cyt.go.genes[[1]]
barcodeplot(res$table$logFC, index=index, labels=c("B.virgin","B.lactating"), 
            main=cyt.go[1])
### Camera gene set enrichment analysis
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
idx <- ids2indices(Mm.c2,id=rownames(y))
BvsL.v <- makeContrasts(B.virgin - L.virgin, levels=design)
cam <- camera(y, idx, design, contrast=BvsL.v, inter.gene.cor=0.01)
options(digits=2)
head(cam,14)

res <- glmQLFTest(fit, contrast=BvsL.v)
barcodeplot(res$table$logFC,
            index=idx[["LIM_MAMMARY_STEM_CELL_UP"]],
            index2=idx[["LIM_MAMMARY_STEM_CELL_DN"]],
            labels=c("B.virgin","L.virgin"),
            main="LIM_MAMMARY_STEM_CELL",
            alpha=1)


