########## Differential Expression Analysis: CTRL-SOR vs NR4ADN-SOR (i.e., the effect of NR4ADN in the context of SOR).

library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0
library(RUVSeq) # RUVSeq_1.20.0

########## load data.
### read count data.
ct = as.matrix(read.csv("data/counts_all.csv", row.names=1, stringsAsFactors=FALSE))
### read sample info.
sample.info = read.csv("data/samples_information.csv", row.names=1, stringsAsFactors=FALSE)
### read gene info.
gene.info = read.csv("data/gene_annotation.csv", row.names=1, stringsAsFactors=FALSE)

identical(rownames(ct), rownames(gene.info))
identical(colnames(ct), rownames(sample.info))

########## create unfiltered SeqExpressionSet object.
# create seqExpressionSet.
set.all = newSeqExpressionSet(counts=ct, phenoData=sample.info, featureData=gene.info)

# make a copy.
set = set.all

### select samples from experiment2.
set = set[ , pData(set)$experiment %in% c("experiment1")]
pData(set) = droplevels(pData(set))

# ordering levels for interpreting sign of fold change.
# positive fold change indicates elevated expression in DN+.
pData(set)$DN = factor(pData(set)$DN, levels=c("DN-", "DN+"))

########## gene filtering.
### filtering with 'filterByExpr'.
y <- DGEList(counts(set), samples=pData(set), genes=fData(set))
keep = filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
set = set[keep,]; rm(keep)

########## EDASeq.
### EDASeq normalization for GC content, and then sequencing depth.
dataWithin <- withinLaneNormalization(set, "gc", which="full", offset=TRUE)
dataNorm <- betweenLaneNormalization(dataWithin, which="upper", offset=TRUE)

### create list to store EDASeq data.
exp1.edaseq = list(set=set, dataWithin=dataWithin, dataNorm=dataNorm)

########## edgeR (QL).
### construct the DGEList and add offsets from EDASeq.
y = DGEList(counts=counts(exp1.edaseq$dataNorm), samples=pData(exp1.edaseq$dataNorm), genes=fData(exp1.edaseq$dataNorm))
y$offset = -offst(exp1.edaseq$dataNorm)

### create design matrix.
design = model.matrix(~ DN, data=pData(exp1.edaseq$dataNorm))

### estimate dispersion and fit model.
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)

### conduct statistical test. 
test = glmQLFTest(fit)

### get differential expression test results.
top = topTags(test, Inf)$table
top$effect.size = sign(top$logFC) * (abs(top$logFC) / (sqrt(1/(top$logCPM + y[rownames(top),]$trended.dispersion + abs(min(top$logCPM)))))   )

### add first pass to 'out' variable.
exp1.results = list(W0=list(set=exp1.edaseq$dataNorm, y=y, fit=fit, test=test, top=top))

########## RUVSeq.
### get residuals from first pass fit (RUVr).
residuals = residuals(exp1.results$W0$fit, type="deviance")

### RUVr
set.ruv = RUVr(exp1.edaseq$dataNorm, rownames(exp1.edaseq$dataNorm), k=6, residuals)
design = model.matrix(~ W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + DN, data=pData(set.ruv))

### repeat edgeR workflow with RUV component(s) included as covariates.
y = DGEList(counts=counts(exp1.edaseq$dataNorm), samples=pData(exp1.edaseq$dataNorm), genes=fData(exp1.edaseq$dataNorm))
y$offset = -offst(exp1.edaseq$dataNorm)
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)
test = glmQLFTest(fit)
top = topTags(test, Inf)$table
top$effect.size = sign(top$logFC) * (abs(top$logFC) / (sqrt(1/(top$logCPM + y[rownames(top),]$trended.dispersion + abs(min(top$logCPM)))))   )

### add first pass to 'out' variable.
exp1.results$W6 = list(set=set.ruv, y=y, fit=fit, test=test, top=top)

### save differential expression results for experiment 1 (i.e., the effect of NR4ADN in the context of SOR).
# write.csv(exp1.results$W6$top, file="results/experiment1_diffexp_results.csv", quote=FALSE)