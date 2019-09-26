###################################################################################################
#
# Plot the first panel of the E. coli substitution figure. 
#
###################################################################################################
library('ggplot2')

## Load the substitution data matrices... 
df <- read.table('01_pb_il_substitution_matrix.tsv', header=T, sep='\t')

# Forgot to reset the index when concatenating PB and IL data
df <- df[,-1]
head(df)

# Also need to fix the apparent offset in the gene position coordinates
df$Position[df$Panel=="lower"] <- df$Position[df$Panel=="lower"] - 2

###################################################################################################
## Plot the panel A figure
df.a <- df[df$Panel=='upper',]
plt <- ggplot(df.a, aes(x=Position, y=Substitutions)) + geom_line()
plt <- plt + theme_bw()
plt <- plt + theme(panel.grid = element_blank())
plt <- plt + ylim(c(-0, 1.2))
plt

plt <- plt + geom_hline(yintercept=seq(0,1, by=1/7),
                        linetype='dotted', alpha=0.5)
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

plt <- plt + xlab('Position along 16S gene') + ylab('Percent Substitutions')

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
plt

pdf('02A_PBIL_ECOLI.pdf', height=3, width=3.5)
plot(plt)
dev.off()

###################################################################################################
## Plot the panel B figure
df.b <- df[df$Panel!='upper',]
plt <- ggplot(df.b, aes(x=Position, y=Substitutions)) + geom_line()
plt <- plt + theme_bw()
plt <- plt + theme(panel.grid = element_blank())
plt <- plt + ylim(c(-0, 1.2))
plt

plt <- plt + geom_hline(yintercept=seq(0,1, by=1/7),
                        linetype='dotted', alpha=0.5)
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

plt <- plt + xlab('Position along 16S gene') + ylab('Percent Substitutions')

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
plt


pl <- plt
pl <- pl + geom_vline(xintercept=c(70, 100))
pl <- pl + geom_vline(xintercept=c(200, 230))
pl <- pl + geom_vline(xintercept=c(245, 275))
pl <- pl + geom_vline(xintercept=c(1000, 1050))
pl

pdf('02B_PBIL_ECOLI.pdf', height=3, width=3.5)
plot(pl)
dev.off()
