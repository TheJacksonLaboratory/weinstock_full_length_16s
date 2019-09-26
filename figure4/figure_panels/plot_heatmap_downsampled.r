
library('ape')
library('gplots')
library('scales')
library('DECIPHER')
library('dendextend')

df <- read.table('downsampled_isolate_error_profiles.tsv', row.names=1, sep='\t', header=TRUE)
df <- as.data.frame(t(df))
head(df)

# remove the OTU prefix from the data frame so that the sequence IDs are compatible with those in the phylogeny
row.names(df) <- gsub('OTU_[0-9]+_', '', row.names(df))
head(df)


# Generate a panel for x-axis showing the intervals of each variable region, based on the position of primer sites in
# B. vulgatus 16S gene sequence
## I1 = 1-9
I1 = rep('white', 9)
# V1 = 9-103
V1 = rep('dark grey', 103-9)
# V2 = 103-249
V2 = rep('dark grey', 249-103)
## I2 = 249-341
I2 = rep('white', 341-249)
# V3 = 341-514
V3 = rep('dark grey', 514-341)
# V4 = 514-685
V4 = rep('dark grey', 685-514)
## I3 = 685-778
I3 = rep('white', 778-685)
# V5 = 778-909
V5 = rep('dark grey', 909-778)
## I4 = 909-968
I4 = rep('white', 968-909)
# V6 = 968-1061
V6 = rep('dark grey', 1061-968)
## I5 = 1061-1098
I5 = rep('white', 1098-1061)
# V7 = 1098-1177
V7 = rep('dark grey', 1177-1098)
## I6 = 1177-1222
I6 = rep('white', 1222-1177)
# V8 = 1222-1370
V8 = rep('dark grey', 1370-1222)
## I7 = 1370-1377
I7 = rep('white', 1377-1370)
# V9 = 1377-1495
V9 = rep('dark grey', 1495-1377)
## I8 = 1495-1500
I8 = rep('white', 1500-1495)


variable.regions <- c(I1, V1, V2, I2, V3, V4, I3, V5, I4, V6, I5, V7, I6, V8, I7, V9, I8)

# Convert the dataframe to be binary
df.1 <- df
df.1[df.1 > 0] <- 1

# Plot the image, but making all the SNPs binary
heatmap.2(as.matrix(df.1), 
          Rowv=FALSE, Colv=FALSE, dendrogram="none", 
          trace="none", 
          key=FALSE, 
          col=colorRampPalette(c('ivory', 'midnightblue'))(n=20),
          labCol=FALSE, rowCol=FALSE,
          margins = c(1,1),
          ColSideColors=alpha(variable.regions, 0.2))

isolate.tree <- read.tree('isolate.tre')
class(isolate.tree)

length(isolate.tree$tip.label)

length(row.names(df.1))

row.names(df.1) %in% isolate.tree$tip.label

isolate.tree.sub <- keep.tip(isolate.tree, row.names(df.1))

length(isolate.tree.sub$tip.label)

identical(isolate.tree.sub$tip.label, row.names(df.1))

plot(isolate.tree.sub, cex=0.4)

write.tree(isolate.tree.sub, file='isolate_sub.tre')

# Using decipher
isolate.tree.dend <- ReadDendrogram('isolate_sub.tre', internalLabels=F)
isolate.tree.dend <- hang.dendrogram(isolate.tree.dend, hang=1.0)
class(isolate.tree.dend)

# Check that the resulting dendrogram has the same labels in the same order as the tree
# For some reason decipher decides to remove underscores and replace them with spaces... 
labels(isolate.tree.dend) <- gsub(' ', '_', labels(isolate.tree.dend))
identical(labels(isolate.tree.dend), isolate.tree.sub$tip.label)

# Check if the row.names of the dataframe and the tip labels are identical
identical(labels(isolate.tree.dend), row.names(df.1) )

identical(sort(labels(isolate.tree.dend)), sort(row.names(df.1)))

df.col <- read.table('../sample_colours.tsv', row.names=1, header=TRUE, sep='\t', comment.char='', 
                     stringsAsFactors = FALSE)
head(df.col)

# Order the colours by the samples in dendrogram
leaf.cols <- df.col[labels(isolate.tree.dend), 'FCol']

isolate.tree.dend.col <- color_branches(isolate.tree.dend, 
                                        k=length(labels(isolate.tree.dend)),
                                        , col=alpha(leaf.cols, 0.8))

plot(isolate.tree.dend.col)

# Manually re-order the rows in the snp matrix, this seems to be necessary (!?)
df.2 <- df.1[labels(isolate.tree.dend),]

identical(row.names(df.1), row.names(df.2))

identical(labels(isolate.tree.dend), row.names(df.2))

heatmap.2(as.matrix(df.2), 
          Rowv=isolate.tree.dend.col, 
          dendrogram='row',
          Colv=FALSE, 
          trace="none", 
          key=FALSE, 
          col=colorRampPalette(c('ivory', 'midnightblue'))(n=20),
          labCol=FALSE,
          labRow=FALSE,
          margins = c(1,1),
          ColSideColors=alpha(variable.regions, 0.2))

pdf('04_Figure4a.pdf', height=4, width=10)
heatmap.2(as.matrix(df.2), 
          Rowv=isolate.tree.dend.col, 
          dendrogram='row',
          Colv=FALSE, 
          trace="none", 
          key=FALSE, 
          col=colorRampPalette(c('ivory', 'midnightblue'))(n=20),
          labCol=FALSE,
          labRow=FALSE,
          margins = c(1,1),
          ColSideColors=alpha(variable.regions, 0.2))
dev.off()

f = c('Erysipelotrichaceae',
'Lactobacillaceae',
'Enterobacteriaceae',
'Ruminococcaceae',
'Bacteroidaceae',
'Clostridiaceae',
'Coriobacteriaceae',
'Lachnospiraceae',
'Bifidobacteriaceae',
'Other')

cl = c("#E41A1C", 
"#984EA3",
"#F781BF",
"#377EB8",
"#FF7F00",
"#FFFF33",
"#A65628",
"#4DAF4A",
"#999999",
"#000000")

pdf('04_Figure4a_legend.pdf')
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", 
       legend =f,
       pch=15, pt.cex=4, cex=1.5, bty='n',
                col = cl)
mtext("Family", at=0.2, cex=2)
dev.off()
