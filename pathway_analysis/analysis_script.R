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
targets <- read.delim("../results/dedup_read_counts.txt", stringsAsFactors=FALSE)
targets

group <- paste(targets$CellType, targets$Status, sep=".")
group <- group[c(1,3,2,4,5,7,6,8)]
group <- factor(group)
table(group)
group

### Preliminary analysis
GenewiseCounts <- read.delim("../results/read_counts.txt", row.names="Geneid")[,-1:-4]
GenewiseCounts <- GenewiseCounts[,c(1,2,4,3,5,6,8,7,9)]
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

### Differential expression analysis
### Testing for differential expression
con <- makeContrasts(
  (DIO.High-DIO.Low)-(ND.High-ND.Low), 
  levels=design)
tr <- glmQLFTest(fit, contrast=con)
#tr <- glmTreat(fit, contrast=con, lfc=log2(1.0))
tr_p <- tr[tr$table$PValue < 0.05,]
tr_p <- tr_p[abs(tr_p$table$logFC) > 1,]

tr_pos <- tr_p[tr_p$table$logFC > 0,]
# pos_tr <- tr
# pos_tr$table <- tr[tr$table$logFC > 0,]
# tr <- glmTreat(fit, contrast=con, lfc=log2(1.0))
pos_tr <- tr[tr$table$logFC > 0,]
pos_tr <- pos_tr[pos_tr$table$PValue < 0.05,]
neg_tr <- tr[tr$table$logFC < 0,]
neg_tr <- neg_tr[neg_tr$table$PValue < 0.05,]

topTags(pos_tr)
topTags(neg_tr)
pos_de <- decideTestsDGE(pos_tr)
neg_de <- decideTestsDGE(neg_tr)
summary(pos_de)
summary(neg_de)

plotMD(pos_tr, status=pos_de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
plotMD(neg_tr, status=neg_de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

### Heatmap clustering
logCPM_pos <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM_pos) <- y$genes$Symbol
colnames(logCPM_pos) <- paste(y$samples$group, 1:2, sep="-")

o <- order(pos_tr$table$PValue)
logCPM_pos <- logCPM_pos[o[1:30],]
logCPM_pos <- t(scale(t(logCPM_pos)))

logCPM_neg <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM_neg) <- y$genes$Symbol
colnames(logCPM_neg) <- paste(y$samples$group, 1:2, sep="-")

o <- order(pos_tr$table$PValue)
logCPM_neg <- logCPM_neg[o[1:30],]
logCPM_neg <- t(scale(t(logCPM_neg)))

library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM_pos, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM_neg, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))


### Pathway analysis
### Gene ontology analysis
go <- goana(tr_p, species="Mm", FDR=1)
topGO(go, n=15)

### KEGG pathway analysis
keg <- kegga(tr_p, species="Mm", FDR=1)
topKEGG(keg, n=15, truncate=34)

### CURRENT WORKING POINT
'''
### FRY gene set tests
library(GO.db)
cyt.go <- c("GO:0032465", "GO:0000281")
term <- select(GO.db, keys=cyt.go, columns="TERM")
term

Rkeys(org.Mm.egGO2ALLEGS) <- cyt.go
cyt.go.genes <- as.list(org.Mm.egGO2ALLEGS)

B.VvsL <- makeContrasts(DIO.High-DIO.Low, levels=design)
fry(y, index=cyt.go.genes, design=design, contrast=B.VvsL)

res <- glmQLFTest(fit, contrast=B.VvsL)
index <- rownames(fit) %in% cyt.go.genes[[1]]
barcodeplot(res$table$logFC, index=index, labels=c("DIO.High","DIO.Low"), 
            main=cyt.go[1])
### Camera gene set enrichment analysis
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
idx <- ids2indices(Mm.c2,id=rownames(y))
High.DIOvsNDCam <- makeContrasts(DIO.High - ND.High, levels=design)
cam <- camera(y, idx, design, contrast=High.DIOvsNDCam, inter.gene.cor=0.01)
options(digits=2)
head(cam,14)

res <- glmQLFTest(fit, contrast=High.DIOvsNDCam)
barcodeplot(res$table$logFC,
            index=idx[["LIM_MAMMARY_STEM_CELL_UP"]],
            index2=idx[["LIM_MAMMARY_STEM_CELL_DN"]],
            labels=c("B.virgin","L.virgin"),
            main="LIM_MAMMARY_STEM_CELL",
            alpha=1)

'''
