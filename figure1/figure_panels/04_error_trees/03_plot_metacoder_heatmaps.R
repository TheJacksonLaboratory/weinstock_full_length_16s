##################################################################################################
#
# Plot the full taxonomic tree with the proportion of failed amplicons as colour
# (assumes working in the same r environment as 00_plot_error_trees.R)
##################################################################################################
library('metacoder')

# Read in the table containing the heatmap info
df.h <- read.table('taxonomic_assignment_percent_failure.tsv', row.names=2, header=T, sep='\t')
head(df.h)

# sort the columns based on the entries in data
str(data)
df.h = df.h[data$taxa,]
head(df.h)

# create a dummy table for plotting v1v9
df.d <- df.h
df.d[1656, 'TaxonCount_v1v9'] <- 100

##################################################################################################
# plot metacoder figures with the colour set as percent error and size as taxon level counts
# Master
pdf('heat_tree_master.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = n_obs,
          node_label = name,
          node_color_range = c("#FFFFFF", "darkorange3", "#4e567d", "gold"),
          edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()


# V1V2
pdf('heat_tree_v1v2.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.h$TaxonCount_v1v2)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()

# V1V3
pdf('heat_tree_v1v3.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.h$TaxonCount_v1v3)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()

# V1V9
pdf('heat_tree_v1v9.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.d$TaxonCount_v1v9)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()

# V3V5
pdf('heat_tree_v3v5.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.h$TaxonCount_v3v5)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()

# V4
pdf('heat_tree_v4.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.h$TaxonCount_v4)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()

# V6V9
pdf('heat_tree_v6v9.pdf', height=5, width=5)
heat_tree(data, node_size = n_obs,
          node_color = (df.h$TaxonCount_v6v9)^2,
          node_label = name,
          node_color_range = c("white", "red"),
          #edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
dev.off()
