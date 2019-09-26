#########################################################################################
#
# Plot the number of OTUs that were generated for different subregions.
#
#########################################################################################
# Fetch a dataframe of the number of OTUs generated 
library("RSQLite")
con = dbConnect(drv=RSQLite::SQLite(), dbname='../00_data/genera_csvdb')
statement = "SELECT * FROM greengenes_otu_summary WHERE Distance=='dist01'"
df = dbGetQuery(con, statement)
row.names(df) <- df$Amplicon
dbDisconnect(con)
df

########################################################################################
# Fetch a dataframe of entropy values... based on a single sequence for each known species.
con = dbConnect(drv=RSQLite::SQLite(), dbname='../00_data/csvdb')
statement = "SELECT * FROM greengenes_insilico_entropy"
df.ent = dbGetQuery(con, statement)
dbDisconnect(con)
row.names(df.ent) <- df.ent$Amplicon
df.ent

# and merge the proprotions with the entropy calculations
df.plt <- merge(df, df.ent, by='row.names')
row.names(df.plt) <- df.plt$Row.names
df.plt <- df.plt[,-1]
df.plt

# Set the per base entropy as a proportion
df.plt$PerBaseEntropy <- df.plt$MedianEntropy/max(df.plt$MedianEntropy)
df.plt


######################################################################################
# make a nice plot 
library('plotrix')
# set up a scaled colour palette
palette <- colorRampPalette(c('azure', 'azure4'))(100)
pdf('../01c_n_otus_generated_per_amplicon.pdf',
    width=4, height=2.5)
par(mar=c(3,3,1,1), las=1)
barplot(df.plt$OTU_Count, 
        names.arg=df$Amplicon,
        ylab = '# Sequences',
        xlab = 'Amplicon',
        col=palette[cut(df.plt$PerBaseEntropy, 100)],
        ylim=c(0,15000), 
        las=1,
        cex.axis=0.5, 
        cex.names=0.5,
        cex.lab=0.5)
box()
# The number of sequences in the greengenes dataset. 
abline(h=13587, lty=5)
dev.off()


# create the legend... based on the values in the dataframe
z = matrix(1:21, nrow=1)
x = 1
y = seq(0.8, 1, 0.01)

pdf('../01c_n_otus_generated_per_amplicon_key.pdf', width=2, height=6)
par(mar=c(0,3,0,0), las=1)
image(x, y, z, col=palette[cut(seq(0.8, 1, 0.01), 100)], axes=FALSE, xlab="", ylab="")
axis(2, las=1)
box()
dev.off()
