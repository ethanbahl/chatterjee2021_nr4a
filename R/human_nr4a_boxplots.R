
#################################################################################### Cerad/Braak boxplots.

########################################################### Read in data from the Aging, Dementia and TBI Study.
# URL: http://aging.brain-map.org

### read expression (FPKM).
fpkm.norm = as.matrix(read.csv("gene_expression_matrix_2016-03-03/fpkm_table_normalized.csv", row.names=1))
fpkm.norm = log2(fpkm.norm+1)

### gene info.
genes = read.csv("gene_expression_matrix_2016-03-03/rows-genes.csv", row.names=1)
nr4a = subset(genes, gene_symbol %in% c("NR4A1", "NR4A2", "NR4A3"))
nr4a = nr4a[order(nr4a$gene_symbol),]

### sample info.
samples = read.csv("gene_expression_matrix_2016-03-03/columns-samples.csv", row.names=1)
rownames(samples) = make.names(rownames(samples))

### donor info.
donors = read.csv("DonorInformation.csv", row.names=1, stringsAsFactors = FALSE)
donors$cerad = factor(donors$cerad)
donors$braak = factor(donors$braak)

samples$cerad = donors[match(samples$donor_name, donors$name), "cerad"]
samples$braak = donors[match(samples$donor_name, donors$name), "braak"]


########################################################### Fit linear models.
### index for hippocampus samples.
hip.idx = which(samples$structure_acronym == "HIP")

############## CERAD.
### fit models.
cerad.fit.list = lapply(rownames(nr4a), function(i) {
   out = glm(fpkm.norm[i,hip.idx] ~ as.numeric(as.character(samples$cerad[hip.idx]))  )
})
names(cerad.fit.list) = nr4a$gene_symbol

### extract model stats.
cerad.fit.stats = do.call('rbind', lapply(cerad.fit.list, function(x) summary(x)$coefficients[2,]))

############## BRAAK.
### fit linear models.
braak.fit.list = lapply(rownames(nr4a), function(i) {
   out = glm(fpkm.norm[i,hip.idx] ~ as.numeric(as.character(samples$braak[hip.idx]))  )
})
names(braak.fit.list) = nr4a$gene_symbol

### extract model stats.
braak.fit.stats = do.call('rbind', lapply(braak.fit.list, function(x) summary(x)$coefficients[2,]))

########################################################### Fit linear models.

############## CERAD.
m = matrix(c(1,1,1,2,3,4), nrow = 2, ncol = 3,byrow = TRUE)
layout(mat = m, heights = c(0.15, 0.85))
par(mar=c(0,0,0,0))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("NR4A family gene expression across Cerad scores"), 
     cex = 2.5, font=2, col = "black")

par(mar=c(5,5,3,2))
for(i in c(1:3)){
    boxplot(fpkm.norm[rownames(nr4a)[i], hip.idx] ~ samples$cerad[hip.idx],
        col=c("grey", RColorBrewer::brewer.pal(9, "YlOrRd")[c(3,5,7)]), las=1, xlab="Cerad score", ylab="log2(FPKM+1)",
        main=paste0(nr4a$gene_symbol[i], ":  P = ", round(fit[i,4], 4) )
    )
    abline(cerad.fit.list[[i]], col=alpha("dodgerblue", 0.7), lwd=3)
}

############## BRAAK.
m = matrix(c(1,1,1,2,3,4), nrow = 2, ncol = 3,byrow = TRUE)
layout(mat = m, heights = c(0.15, 0.85))
par(mar=c(0,0,0,0))

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("NR4A family gene expression across Braak stages"), 
     cex = 2.5, font=2, col = "black")

par(mar=c(5,5,3,2))
for(i in c(1:3)){
    boxplot(fpkm.norm[rownames(nr4a)[i], hip.idx] ~ samples$braak[hip.idx],
        col=c("grey", RColorBrewer::brewer.pal(9, "YlOrRd")[c(2:7)]), las=1, xlab="Braak stage", ylab="log2(FPKM+1)",
        main=paste0(nr4a$gene_symbol[i], ":  P = ", round(fit[i,4], 4) )
    )
    abline(braak.fit.list[[i]], col=alpha("dodgerblue", 0.7), lwd=3)
}
