##################################################################################################
#
# Plot the entropy across 16S gene
#
##################################################################################################
library('ggplot2')

df = read.table('../00_data/V1_V9_smoothed_entropy.tsv', header=T, sep='\t')
head(df)

# The basic entropy plot. 
plt <- ggplot(df, aes(x=Base_Position, y=Entropy)) + geom_line()
plt


##################################################################################################
# The same plot, but with panels. 
# parameters
# region title offset
n = 0.01
# region title alpha
ra=0.8
# region title colour
rc='dark grey'
# panel alpha
pa=0.2
#panel colour
pc='grey'
#amplicon colour
ac = 'dark red'


df$value_n <- df$Entropy / max(df$Entropy, na.rm=T)
plt <- ggplot(df)
plt <- plt + geom_line(aes(x=Base_Position, y=value_n))
plt <- plt + theme_bw() + ylim(-0.1, 1.02)
# add colour panels
plt <- plt + annotate('rect', xmin=9, xmax=104, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=103, xmax=249, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=341, xmax=515, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=514, xmax=685, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=778, xmax=909, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=968, xmax=1061, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=1098, xmax=1177, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=1222, xmax=1370, ymin=0, ymax=1, fill=pc, alpha=pa)
plt <- plt + annotate('rect', xmin=1377, xmax=1492, ymin=0, ymax=1, fill=pc, alpha=pa)

# add individual bars
#V1
plt <- plt + annotate('segment', x=9, xend=104, y=1, yend=1, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=9 + (104-9)/2, y=1+n, label='V1', cex=3)
#V2
plt <- plt + annotate('segment', x=103, xend=249, y=1+n, yend=1+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=103 + (249-103)/2, y=1+n+n, label='V2', cex=3)
#V3
plt <- plt + annotate('segment', x=341, xend=515, y=1, yend=1, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=341 + (515-341)/2, y=1+n, label='V3', cex=3)
#V4
plt <- plt + annotate('segment', x=514, xend=685, y=1+n, yend=1+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=514 + (685-514)/2, y=1+n+n, label='V4', cex=3)
#V5
plt <- plt + annotate('segment', x=778, xend=909, y=1, yend=1, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=778 + (909-778)/2, y=1+n, label='V5', cex=3)
#V6
plt <- plt + annotate('segment', x=968, xend=1061, y=1+n, yend=1+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=968 + (1061-968)/2, y=1+n+n, label='V6', cex=3)
#V7
plt <- plt + annotate('segment', x=1098, xend=1177, y=1, yend=1, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1098 + (1177-1098)/2, y=1+n, label='V7', cex=3)
#V8
plt <- plt + annotate('segment', x=1222, xend=1370, y=1+n, yend=1+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1222 + (1370-1222)/2, y=1+n+n, label='V8', cex=3)
#V9
plt <- plt + annotate('segment', x=1377, xend=1495, y=1, yend=1, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1377 + (1492-1377)/2, y=1+n, label='V9', cex=3)

# Add amplicon regions
#V1-V2
plt <- plt + annotate('segment', x=28, xend=341, y=-0.01, yend=-0.01, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=28 + (341-28)/2, y=-0.01-n, label='V1-V2', colour=ac, cex=3)
#V1-V3
plt <- plt + annotate('segment', x=28, xend=517, y=-0.02-n, yend=-0.02-n, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=28 + (517-28)/2, y=-0.02-2*n, label='V1-V3', colour=ac, cex=3)
#V3-V5
plt <- plt + annotate('segment', x=358, xend=908, y=-0.03-2*n, yend=-0.03-2*n, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=358 + (908-358)/2, y=-0.03-3*n, label='V3-V5', colour=ac, cex=3)
#V4
plt <- plt + annotate('segment', x=534, xend=786, y=-0.04-3*n, yend=-0.04-3*n, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=534 + (786-534)/2, y=-0.04-4*n, label='V4', colour=ac, cex=3)
#V6-V9
plt <- plt + annotate('segment', x=985, xend=1491, y=-0.02-n, yend=-0.02-n, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=985 + (1491-985)/2, y=-0.02-2*n, label='V6-V9', colour=ac, cex=3)
#V6-V9
plt <- plt + annotate('segment', x=28, xend=1491, y=-0.05-4*n, yend=-0.05-4*n, colour=ac, alpha=ra)
plt <- plt + annotate('text', x=28 + (1491-28)/2, y=-0.05-5*n, label='V1-V9', colour=ac, cex=3)

plt <- plt + theme(panel.grid = element_blank())

plot(plt)

##################################################################################################
# Save the plot... 
pdf('../01_entropy_plot.pdf', height=4, width=8)
plot(plt)
dev.off()


