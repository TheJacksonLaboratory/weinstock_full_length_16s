###################################################################################################
#
# Plot the number sequences that were not able to be correctly identified to species level
#
###################################################################################################

# Fetch a dataframe of a breakdown of taxonomic assignment success at 80% confidence
# This is data for all sequences (clustered at 99% identity) that have genus level classification
library("RSQLite")
con = dbConnect(drv=RSQLite::SQLite(), dbname='../00_data/genera_csvdb')
statement = "SELECT * FROM greengenes_insilico_tax_assignment_summary_table WHERE RDP_Threshold = 80"
df = dbGetQuery(con, statement)
dbDisconnect(con)
df

# calculate the % sequences that were not assigned to species level
rownames(df) <- df$bin
df$bin <- NULL
df$RDP_Threshold <- NULL
df.sum <- sweep(df, 2, apply(df, 2, sum), '/')
df.sum <- df.sum[rownames(df.sum) != 'Species', ]
pptn <- apply(df.sum, 2, sum)*100
pptn
class(pptn)

# Convert the proportion of species-level failures to a dataframe
df.pptn <- as.data.frame(pptn)


#################################################################################################
# Fetch a dataframe of entropy values... based on a single sequence for each known species.
con = dbConnect(drv=RSQLite::SQLite(), dbname='../00_data/csvdb')
statement = "SELECT * FROM greengenes_insilico_entropy"
df.ent = dbGetQuery(con, statement)
dbDisconnect(con)
row.names(df.ent) <- df.ent$Amplicon
df.ent

# and merge the proprotions with the entropy calculations
df.plt <- merge(df.pptn, df.ent, by='row.names')
row.names(df.plt) <- df.plt$Row.names
df.plt <- df.plt[,-1]
df.plt


################################################################################################
## Make a nice looking plot

# Convert the median entropy to a proportion of the maximum encountered
df.plt$perBaseEntropy <- df.plt$MedianEntropy/(max(df.plt$MedianEntropy))
df.plt

pallete <- colorRampPalette(c('azure', 'azure4'))(100)

par(mar=c(3,3,1,1), las=1)
# make a nice plot 
barplot(df.plt$pptn, 
        names.arg=df.plt$Amplicon,
        ylab = '% Sequences',
        xlab = 'Amplicon',
        col=pallete[cut(df.plt$perBaseEntropy, 100)],
        ylim=c(0, 60),
        las=1,
        cex.axis=0.5, 
        cex.names=0.5,
        cex.lab=0.5)
box()

pdf('../01b_pptn_sequences_unsucceffully_assigned.pdf',
    width=4, height=2.5)
barplot(df.plt$pptn, 
        names.arg=df.plt$Amplicon,
        ylab = '% Sequences',
        xlab = 'Amplicon',
        col=pallete[cut(df.plt$perBaseEntropy, 100)],
        ylim=c(0, 60),
        las=1,
        cex.axis=0.5, 
        cex.names=0.5,
        cex.lab=0.5)
box()
dev.off()
