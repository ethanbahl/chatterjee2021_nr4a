########## NORMALIZATION, EXPERIMENT 2

library(EDASeq) # EDASeq_2.20.0
library(scatterplot3d) # scatterplot3d_0.3-41
library(RColorBrewer) # RColorBrewer_1.1-2
library(svglite) # svglite_1.2.3

### grab the 3 sets.
set = exp2.edaseq$set
dataWithin = exp2.edaseq$dataWithin
dataNorm = exp2.edaseq$dataNorm

### set up plotting objects.
colors = brewer.pal(8, "Dark2")[c(3,1)]
legend.cols = colors
group_color = factor(pData(set)$exposure, levels=c("homecage", "SOR"))
plot.cols = colors[group_color]
use_3d_pca_plots = FALSE

##### reproducing SVD data generated by EDASeq::plotPCA().
# set.
Y <- apply(log(counts(set)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s1 <- svd(Y)
percent <- s1$d^2/sum(s1$d^2)*100
labs1 <- sapply(seq_along(percent), function(i) {
    paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
})
# dataWithin.
Y <- apply(log(normCounts(dataWithin)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s2 <- svd(Y)
percent <- s2$d^2/sum(s2$d^2)*100
labs2 <- sapply(seq_along(percent), function(i) {
    paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
})
# dataNorm.
Y <- apply(log(normCounts(dataNorm)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s3 <- svd(Y)
percent <- s3$d^2/sum(s3$d^2)*100
labs3 <- sapply(seq_along(percent), function(i) {
    paste("PC ", i, " (", round(percent[i], 2), "%)", sep="")
})

################################################# PLOTTING.
#svglite("figures/edaseq_exp2.svg", height=16, width=18, system_fonts = list(sans = "Arial"))
cex.lab = 1.75
cex.main = 2
cex.axis=1.3

m = matrix(
    c(
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
    ),
    byrow=TRUE,
    nrow=4
)
layout(mat = m, heights = c(0.04, rep(0.24, 3)))
par(mar=c(0,0,0,0))#, family="Arial")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("raw"), 
     cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC normalization"), 
     cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC / depth normalization"), 
     cex = 3, font=2, col = "black")  
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

par(
  mar=c(5,5,4,2)
)

# GC content parameters.
ylab = "log(gene counts)"
xlab = "GC content"
ylim=c(2,8)
xlim=c(min(fData(set)$gc),max(fData(set)$gc))
gc.col = "#333333"
c1 = alpha(gc.col, 0.8)
c2 = alpha(gc.col, 0.2)
c3 = alpha(gc.col, 0.8)

########## GC bias plots.
biasPlot(set, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(set)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

biasPlot(dataWithin, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(group_color), fill=legend.cols, bty="n", cex=1.1)
boxplot(fData(dataWithin)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(dataWithin)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

biasPlot(dataNorm, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(dataNorm)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(dataNorm)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

########## RLE
ylim=c(-2, 2)
plotRLE(set, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

plotRLE(dataWithin, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

plotRLE(dataNorm, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")


########## PCA
xlim=ylim=zlim = c(-1,1)

if(use_3d_pca_plots == FALSE) {
    plotPCA(set, col=plot.cols, pch=16, xlim=xlim, ylim=ylim, main="PCA plot", labels=FALSE)
    legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

    plotPCA(dataWithin, col=plot.cols, pch=16, xlim=xlim, ylim=ylim, main="PCA plot", labels=FALSE)
    legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

    plotPCA(dataNorm, col=plot.cols, pch=16, xlim=xlim, ylim=ylim, main="PCA plot", labels=FALSE)
    legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")
} else {

    s = s1
    #s$u = s$u[,-2]
    scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC3", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
    text(x = 3.5, y = -0.85, labels="PC2", srt = 55, cex=cex.lab-0.05)
    #mtext(text="", side=1, line=3, at = c(0), cex=0.8)

    s = s2
    #s$u = s$u[,-2]
    scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC3", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
    text(x = 3.5, y = -0.85, labels="PC2", srt = 55, cex=cex.lab-0.05)
    #mtext(text="", side=1, line=3, at = c(0), cex=0.8)

    s = s3
    #s$u = s$u[,-1]
    scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC3", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
    text(x = 3.5, y = -0.85, labels="PC2", srt = 55, cex=cex.lab-0.05)
    #mtext(text="", side=1, line=3, at = c(0), cex=0.8)

}

#dev.off()

###################################################################################################
############################# RUVSEQ
########## prepare data for plotting.
# get differential expression results, with and without RUV normalization.
tt = exp2.results$W0$top
tt.ruv = exp2.results$W1$top
set = exp2.results$W0$set
set.ruv = exp2.results$W1$set

### set up plotting objects.
colors = brewer.pal(8, "Dark2")[c(3,1)]
legend.cols = colors
group_color = factor(pData(set)$exposure, levels=c("homecage", "SOR"))
plot.cols = colors[group_color]
use_3d_pca_plots = TRUE

# get PCA plot data using adapted source code from EDASeq's plotPCA().
object = normCounts(set)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S <- svd(Y)
object = normCounts(set.ruv)
Y <- apply(log(object+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
S.ruv <- svd(Y)

########## plot set up.
# graphical parameters.
cex.lab = 1.75
cex.main = 2
cex.axis=1.3

# layout matrix.
m = matrix(
    c(
        1,2,
        3,4,
        5,6,
        7,8,
        9,10
    ),
    byrow=TRUE,
    nrow=5
)

########## save plot as SVG.
#svglite("figures/ruvseq_exp2.svg", height=20, width=12, system_fonts = list(sans = "Arial"))
layout(mat = m, heights = c(0.04, rep(0.24, 4)))
par(mar=c(0,0,0,0))
# column 1 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth normalization"), cex = 3, font=2, col = "black")
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# column 2 heading.
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, paste("GC / depth / RUV normalization"), cex = 3, font=2, col = "black")  
segments(x0=0.05, x1=0.95, y0=0.2, y1=0.2, lwd=2.5, col="black")
# update parameters for the main plots in the figure.
par(mar=c(5,5,3,2))

########## PValue distribution plots.
# PValue distribution of differential expression results (without RUV normalization).
hist(tt$PValue, xlab="P-Value", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")
# PValue distribution of differential expression results (with RUV normalization).
hist(tt.ruv$PValue, xlab="P-Value", ylab="number of genes", main="distribution of P-values", cex.main=cex.main, ylim=c(0,3000), cex.lab=cex.lab, cex.axis=cex.axis, col="#d3d3d3")

########## mean-difference plots.
# mean-difference plot (without RUV normalization).
sig = tt[,"FDR"] <= 0.05
plot(tt[-which(sig),"logCPM"], tt[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt$logCPM), max(tt$logCPM)),
    ylim=c(min(tt$logFC), max(tt$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt[sig, "logCPM"],tt[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, "red", "dodgerblue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "dodgerblue", "#333333"), bty="n", cex=1.2)
# mean-difference plot (with RUV normalization).
sig = tt.ruv[,"FDR"] <= 0.05
plot(tt.ruv[-which(sig),"logCPM"], tt.ruv[-which(sig),"logFC"],
    pch=16,cex=0.3,
    ylab="logFC (relative to control)",
    xlab="average log(CPM)",
    xlim=c(min(tt.ruv$logCPM), max(tt.ruv$logCPM)),
    ylim=c(min(tt.ruv$logFC), max(tt.ruv$logFC)),
    main="mean-difference plot",
    col=scales::alpha("#333333", 0.6),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)
points(tt.ruv[sig, "logCPM"],tt.ruv[sig, "logFC"],
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt.ruv[sig,"logFC"] > 0, "red", "dodgerblue"), 0.8)
)
abline(h=0, lty=2, col="#AFAFAF")
legend("topright", legend=c("upregulated", "downregulated", "not significant"), fill=c("red", "dodgerblue", "#333333"), bty="n", cex=1.2)

########## RLE plots.
# RLE plot (without RUV normalization).
ylim=c(-2, 2)
plotRLE(set, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

# RLE plot (with RUV normalization).
plotRLE(set.ruv, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topright", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

########## 3D PCA plots.
xlim=ylim=zlim = c(-1,1)
# 3D PCA plot (without RUV normalization).
s = S
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC3", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.5, y = -0.85, labels="PC2", srt = 55, cex=cex.lab-0.05)

# 3D PCA plot (with RUV normalization).
s = S.ruv
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC3", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.5, y = -0.85, labels="PC2", srt = 55, cex=cex.lab-0.05)

# save plot.
#dev.off()