########################################################################################################################################
############################################################################################################# METANALYSIS: SCATTER PLOT.
library(grid)
library(ggplot2) # ggplot2_3.3.0
library(ggrepel) # ggrepel_0.8.2      
library(RColorBrewer) # RColorBrewer_1.1-2
library(svglite) # svglite_1.2.3

# read in results data.
res.exp1 = read.csv("results/experiment1_diffexp_results.csv", stringsAsFactors=FALSE, row.names=1)
res.exp2 = read.csv("results/experiment2_diffexp_results.csv", stringsAsFactors=FALSE, row.names=1)
rownames(res.exp2) = res.exp2$ensembl_gene_id
rownames(res.exp1) = res.exp1$ensembl_gene_id

### filter rows to common tested gene list, then merge.
int.genes = intersect(rownames(res.exp2), rownames(res.exp1))
res = cbind(res.exp2[int.genes,], res.exp1[int.genes, c("logFC", "logCPM", "F", "PValue", "FDR", "effect.size")])
colnames(res)[16:21] = paste0("exp2.", colnames(res)[16:21])
colnames(res)[22:27] = paste0("exp1.", colnames(res)[22:27])
# switch logFC sign for experiment1 so CTRL+SOR samples have higher expression of positive fold change genes in both lists.
res$exp1.logFC = -1 * res$exp1.logFC
# create signed log10(FDR) for visualization.
res$exp2.slfdr = -log10(res$exp2.FDR)*sign(res$exp2.logFC)
res$exp1.slfdr = -log10(res$exp1.FDR)*sign(res$exp1.logFC)
# create max FDR between both analyses for mapping to color/size.
res$max.fdr = pmax(res$exp2.FDR, res$exp1.FDR)

### set up color palette.
palette.rdylbu = (RColorBrewer::brewer.pal(11, "RdYlBu"))
max.fdr = 1 - res$max.fdr

# if signs are concordant, use RdYlBu palette. Otherwise, use gray.
sign.concord = sign(res$exp1.slfdr) * (sign(res$exp2.slfdr) )
sign.direction = sign(res$exp1.slfdr)
color.data = max.fdr * sign.direction

### function to map numeric to color vectors. Found here: https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r/18749392, but unsure of original source.
map2color = function (x, pal, limits = NULL) 
{
    if (is.null(limits)) 
        limits = range(x)
    pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 
        1), all.inside = TRUE)]
}

max.fdr.colors = map2color(x=color.data, pal=rev(colorRampPalette(palette.rdylbu)(101)), limits=c(-1, 1))
# change non-concordant genes to grey.
idx.discordant = sign.concord == -1
max.fdr.colors[idx.discordant] = scales::alpha(rep("#AAAAAA", sum(idx.discordant)), max.fdr[idx.discordant])

p = ggplot(res, aes(exp1.slfdr, exp2.slfdr)) +
  geom_point(color = ifelse(res$exp2.FDR <= 0.05 & res$exp1.FDR <= 0.05, "red", "black"))

p1 = p + geom_label_repel(
    data = subset(res, exp1.slfdr >= -log10(0.05) & exp2.slfdr >= -log10(0.05)),
    aes(label=gene_name, size=0.1),
    #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
    segment.size  = 0.1,
    segment.color = "grey50",
    force=.1,
    max.iter=10000,
   # direction     = "x",
   seed=777,
   xlim=c(4.65, 5.45),
   ylim=c(0.4, 2.95),
   box.padding=0.18,
   label.padding=0.05,
   point.padding=0.09,
   label.size=0.1,
   direction="y"
  )
  
  p2 = p + geom_label_repel(
    data = subset(res, exp1.slfdr <= log10(0.05) & exp2.slfdr <= log10(0.05)),
    aes(label=gene_name, size=0.2),
    #nudge_y       = 36 - subset(dat, mpg > 30)$mpg,
    segment.size  = 0.1,
    segment.color = "grey50",
    force=0.4,
    max.iter=10000,
   # direction     = "x",
   seed=777,
   xlim=c(-2.45, log10(0.05)),
   ylim=c(-2.95, log10(0.05)),
   box.padding=0.7,
   label.padding=0.3,
   point.padding=0.1,
   label.size=0.5
  ) 
  
# Function: get the x and y positions of a single ggrepel label
get.xy.pos.labs <- function(n) {
  grb <- grid.get(n)
  data.frame(
  x = xrg[1]+diff(xrg)*convertX(grb$x, "native", valueOnly = TRUE),
  y = yrg[1]+diff(yrg)*convertY(grb$y, "native", valueOnly = TRUE)
  )
}

# Get x and y plot ranges 
xrg <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
yrg <- ggplot_build(p1)$layout$panel_params[[1]]$y.range

# manual extraction of coordinates.
p1
grid.force()
kids1 <- childNames(grid.get("labelrepeltree", grep = TRUE))
bulk.lab1 <- do.call(rbind, lapply(kids1, get.xy.pos.labs))
bulk.lab1 = cbind(bulk.lab1, subset(res, exp1.slfdr >= -log10(0.05) & exp2.slfdr >= -log10(0.05))[c(10,28,29)])
p2
grid.force()
kids2 <- childNames(grid.get("labelrepeltree", grep = TRUE))
bulk.lab2 <- do.call(rbind, lapply(kids2, get.xy.pos.labs))
bulk.lab2 = cbind(bulk.lab2, subset(res, exp1.slfdr <= log10(0.05) & exp2.slfdr <= log10(0.05))[c(10,28,29)])

# getting labels, formatting coordinates.
if(exists("bulk.lab2")) {
    bulk.lab = rbind(bulk.lab1, bulk.lab2)
} else {
    bulk.lab = bulk.lab1
}
colnames(bulk.lab) = c("x1", "y1", "label", "y0", "x0")
bulk.lab = bulk.lab[, c(3,4,1,5,2)]
bulk.lab$distance = sqrt( (bulk.lab$x1-bulk.lab$x0)^2 + (bulk.lab$y1-bulk.lab$y0)^2 )
d = bulk.lab$distance

# prep color data.
channels.min = 0.1
channel.coef = 1.2
palette.rdylbu = (RColorBrewer::brewer.pal(11, "RdBu"))[-6]
color.data = max.fdr * sign.direction
max.fdr.colors = map2color(x=color.data, pal=rev(colorRampPalette(palette.rdylbu)(101)), limits=c(log10(0.05), -log10(0.05)))

# change non-concordant genes to grey.
idx.discordant = sign.concord == -1
max.fdr.colors[idx.discordant] = scales::alpha(rep("#AAAAAA" ,  sum(idx.discordant)), pmax(max.fdr[idx.discordant], channels.min))

## shifting the Fos label up slightly.
bulk.lab[13,5] = bulk.lab[12,5] + ((bulk.lab[11,5] - bulk.lab[12,5])/ 2)

# svglite("figures/bulk_bulk_scatterplot.svg", height=9.5, width=10.75, system_fonts = list(sans = "Arial"))
par(mar=c(6,6,6,2), xpd=FALSE)
plot( 
    x=res$exp1.slfdr,
    y=res$exp2.slfdr,
    pch=16, cex=0, cex.lab=1.2, cex.main=1.8,
    col=scales::alpha(max.fdr.colors, pmax(max.fdr, channels.min)*channel.coef),
    xlim=c(-2.85,5.1), ylim=c(-2.35,3),
    ylab=expression("signed -log"[10]*"(FDR) :  CTRL-SOR relative to CTRL-homecage"),
    xlab=expression("signed -log"[10]*"(FDR) :  CTRL-SOR relative to NR4ADN-SOR"),
    main="Control SOR Vs Control Homecage: RNA seq\nNR4ADN SOR Vs Control SOR: RNA seq"
)

# add gridlines.
abline(h=0, v=0, lty=1, lwd=1.5, col=scales::alpha("#333333", 0.8))
abline(h=c(log10(0.05), -log10(0.05)), v=c(log10(0.05), -log10(0.05)), col=scales::alpha("#333333", 0.8), lty=3, lwd=1.2)

# add segments for text boxes.
seg.mindist = 0.27
for(i in c(1:nrow(bulk.lab))) {
    if(bulk.lab$distance[i] >= seg.mindist)  segments(x0=bulk.lab$x0[i], y0=bulk.lab$y0[i], x1=4.3, y1=bulk.lab$y1[i], lwd=0.4, col=scales::alpha("#444444", 0.3))
}

# add points.
points(
    x=res$exp1.slfdr,
    y=res$exp2.slfdr,
    pch=16, cex=pmax(max.fdr, channels.min)*channel.coef,
    col=scales::alpha(max.fdr.colors, pmax(max.fdr, channels.min)*channel.coef)
)

# add labels.
text(x=4.4, y=bulk.lab$y1, labels=bulk.lab$label, cex=1.05, col="#444444", font=2, adj=0)

# add legend.
legend(3.25, -1.27, expression(bold("FDR = 0.05")), lwd=1.3, lty=3, cex=0.9, bty="n")
legend(2.8, -1.7, c("discordant between experiments", "concordant: higher in WT-SOR", "concordant: lower in WT-SOR"), title=expression(bold("direction of effect")), fill=c("#AAAAAA", palette.rdylbu[c(3,8)]), cex=0.9)

#dev.off()