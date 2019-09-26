###########################################################################################
#
# Plot the alignment profiles of bacteroides vulgatus mock vs. in vivo
#
###########################################################################################
library('ggplot2')

df.a <- read.table('cp000139_snp_profile.tsv', sep='\t', stringsAsFactors=F, header=T,
                   row.names=1)
head(df.a)

df.b <- read.table('cp013020_snp_profile.tsv', sep='\t', stringsAsFactors=F, header=T,
                   row.names=1)
head(df.b)

########################################################################################
# Note that the primer coordinates were taken from the original bacteroides vulgatus
# entropy plot that formed part of previous figure 3 panels
# parameters
# region title offset
n = 1
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
########################################################################################
# Check the length of the x-axis...
head(df.a)
tail(df.a)


head(df.b)
tail(df.b)

df.a <- df.a[1:1500,]
df.b <- df.b[1:1500,]
########################################################################################

#  plot

# Note that the original pile-ups determined by Shana included the reference sequence
# in the calculation. The SNP frequencies are therefore out of 8, rather than 7. 
head(df.a)
df.a$b_vulgatus_adj = (df.a$bacteroides_vulgatus*8)/7
head(df.a)

plt <- ggplot(df.a, aes(x=BasePosition, y=b_vulgatus_adj)) + geom_line()
plt <- plt + theme_bw()
plt <- plt + theme(panel.grid = element_blank())
plt <- plt + ggtitle('B. vulgatus ATCC 8482') + xlab('') + ylab('')
plt <- plt + ylim(c(0, 102))
plt

plt <- plt + geom_hline(yintercept=seq(-100,100, by=100/7),
                        linetype='dotted', alpha=0.5)
plt

# add individual bars
#V1
plt <- plt + annotate('rect', xmin=31, xmax=110, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=31, xend=110, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=9 + (110-31)/2, y=100+n, label='V1', cex=3)
#V2
plt <- plt + annotate('rect', xmin=126, xmax=258, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=126, xend=258, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=126 + (258-126)/2, y=100+n+n, label='V2', cex=3)
#V3
plt <- plt + annotate('rect', xmin=365, xmax=519, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=365, xend=519, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=365 + (519-365)/2, y=100+n, label='V3', cex=3)
#V4
plt <- plt + annotate('rect', xmin=537, xmax=689, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=537, xend=689, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=537 + (689-537)/2, y=100+n+n, label='V4', cex=3)
#V5
plt <- plt + annotate('rect', xmin=801, xmax=910, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=801, xend=910, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=801 + (910-801)/2, y=100+n, label='V5', cex=3)
#V6
plt <- plt + annotate('rect', xmin=986, xmax=1059, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=986, xend=1059, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=986 + (1059-986)/2, y=100+n+n, label='V6', cex=3)
#V7
plt <- plt + annotate('rect', xmin=1114, xmax=1176, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1114, xend=1176, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1114 + (1176-1114)/2, y=100+n, label='V7', cex=3)
#V8
plt <- plt + annotate('rect', xmin=1240, xmax=1369, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1240, xend=1369, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1240 + (1369-1240)/2, y=100+n+n, label='V8', cex=3)
#V9
plt <- plt + annotate('rect', xmin=1394, xmax=1482, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1394, xend=1482, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1394 + (1482-1394)/2, y=100+n, label='V9', cex=3)
plt

pdf('03_bvulgatus_atcc8482_only_panel_corrected.pdf', height=3, width=5.5)
plot(plt)
dev.off()



########################################################################################
# plot
# Note that the original pile-ups determined by Shana included the reference sequence
# in the calculation. The SNP frequencies are therefore out of 8, rather than 7. 
head(df.b)
df.b$b_vulgatus_adj = (df.b$bacteroides_vulgatus*-8)/7
head(df.b)

# Note also that the first SNV at position 2 only has six sequences covering it, 
# I'm therefore masking it in this plot...
df.b[df.b$BasePosition==2, 'b_vulgatus_adj'] <- 0


plt <- ggplot(df.b, aes(x=BasePosition, y=b_vulgatus_adj)) + geom_line()
plt <- plt + theme_bw()
plt <- plt + theme(panel.grid = element_blank())
plt <- plt + ggtitle('B. vulgatus mpk') + xlab('') + ylab('')
plt <- plt + ylim(c(0, 102))
plt

plt <- plt + geom_hline(yintercept=seq(-100,100, by=100/7),
                        linetype='dotted', alpha=0.5)
plt

# add individual bars
#V1
plt <- plt + annotate('rect', xmin=31, xmax=110, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=31, xend=110, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=9 + (110-31)/2, y=100+n, label='V1', cex=3)
#V2
plt <- plt + annotate('rect', xmin=126, xmax=258, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=126, xend=258, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=126 + (258-126)/2, y=100+n+n, label='V2', cex=3)
#V3
plt <- plt + annotate('rect', xmin=365, xmax=519, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=365, xend=519, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=365 + (519-365)/2, y=100+n, label='V3', cex=3)
#V4
plt <- plt + annotate('rect', xmin=537, xmax=689, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=537, xend=689, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=537 + (689-537)/2, y=100+n+n, label='V4', cex=3)
#V5
plt <- plt + annotate('rect', xmin=801, xmax=910, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=801, xend=910, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=801 + (910-801)/2, y=100+n, label='V5', cex=3)
#V6
plt <- plt + annotate('rect', xmin=986, xmax=1059, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=986, xend=1059, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=986 + (1059-986)/2, y=100+n+n, label='V6', cex=3)
#V7
plt <- plt + annotate('rect', xmin=1114, xmax=1176, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1114, xend=1176, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1114 + (1176-1114)/2, y=100+n, label='V7', cex=3)
#V8
plt <- plt + annotate('rect', xmin=1240, xmax=1369, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1240, xend=1369, y=100+n, yend=100+n, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1240 + (1369-1240)/2, y=100+n+n, label='V8', cex=3)
#V9
plt <- plt + annotate('rect', xmin=1394, xmax=1482, ymin=0, ymax=100, fill=pc, alpha=pa)
plt <- plt + annotate('segment', x=1394, xend=1482, y=100, yend=100, colour=rc, alpha=ra)
plt <- plt + annotate('text', x=1394 + (1482-1394)/2, y=100+n, label='V9', cex=3)
plt

pdf('03_bvulgatus_mpk_only_panel_corrected.pdf', height=3, width=5.5)
plot(plt)
dev.off()




