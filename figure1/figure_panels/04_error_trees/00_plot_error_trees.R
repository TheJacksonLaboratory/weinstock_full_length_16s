##################################################################################################
#
# Plot the full taxonomic tree with the proportion of failed amplicons as colour
#
##################################################################################################
library('ape')
library('metacoder')

# Read in the metacoder tax objects for all, and for each subregion
fasta_path <- file.path('greengenes_metacoder.fasta')
seqs <- ape::read.FASTA(fasta_path)
data <- extract_taxonomy(seqs, regex = "^(.*)\\t(.*)",
                         key = c(id = "obs_info", "class"),
                         class_sep = ";")
#V1V2
fasta_path_v1v2 <- file.path('V1_V2_amplicons_annotated.fasta')
seqs_v1v2 <- ape::read.FASTA(fasta_path_v1v2 )
data.v1v2 <- extract_taxonomy(seqs_v1v2, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")
#V1V3
fasta_path_v1v3 <- file.path('V1_V3_amplicons_annotated.fasta')
seqs_v1v3 <- ape::read.FASTA(fasta_path_v1v3 )
data.v1v3 <- extract_taxonomy(seqs_v1v3, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")

#V1V9
fasta_path_v1v9 <- file.path('V1_V9_amplicons_annotated.fasta')
seqs_v1v9 <- ape::read.FASTA(fasta_path_v1v9 )
data.v1v9 <- extract_taxonomy(seqs_v1v9, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")
#V3V5
fasta_path_v3v5 <- file.path('V3_V5_amplicons_annotated.fasta')
seqs_v3v5 <- ape::read.FASTA(fasta_path_v3v5 )
data.v3v5 <- extract_taxonomy(seqs_v3v5, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")

#V4
fasta_path_v4 <- file.path('V4_amplicons_annotated.fasta')
seqs_v4 <- ape::read.FASTA(fasta_path_v4 )
data.v4 <- extract_taxonomy(seqs_v4, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")
#V6V9
fasta_path_v6v9 <- file.path('V6_V9_amplicons_annotated.fasta')
seqs_v6v9 <- ape::read.FASTA(fasta_path_v6v9 )
data.v6v9 <- extract_taxonomy(seqs_v6v9, regex = "^(.*)\\t(.*)",
                              key = c(id = "obs_info", "class"),
                              class_sep = ";")

##################################################################################################
# Try plotting the master tree... 
heat_tree(data, node_size = n_obs,
          node_color = n_obs,
          node_label = name,
          node_color_range = c("#FFFFFF", "darkorange3", "#4e567d", "gold"),
          edge_color='grey',
          node_label_max = 100,
          node_label_size=20)
heat_tree(data.v1v2, node_size=n_obs)
##################################################################################################
### Create a table of the taxon counts for each subregion.

# Have a look at the taxmap datastructure
str(data)
data$taxon_data$taxon_ids
data$taxon_data$name
n_obs(data)

hierarchies(data)

identical(names(n_obs(data)), data$taxon_data$taxon_ids)

# Fetch the taxonomic levels and their counts...
df.master <- data.frame('TaxonName'=data$taxon_data$name,
                        'TaxonCount'=n_obs(data),
                        'TaxonHierarchy'=hierarchies(data),
                        row.names=data$taxon_ids,
                        stringsAsFactors=FALSE)
write.table(df.master, file='master_hierachy.tsv', sep='\t', quote=F)
head(df.master)
nrow(df.master)

##Fetch for their sub-regions
#V1V2
df.v1v2 <- data.frame('TaxonName_v1v2'=data.v1v2$taxon_data$name, 
                      'TaxonCount_v1v2'=n_obs(data.v1v2),
                      'TaxonHierarchy_v1v2'=hierarchies(data.v1v2),
                      row.names=data.v1v2$taxon_ids,
                      stringsAsFactors=FALSE)
write.table(df.v1v2, file='v1v2_hierachy.tsv', sep='\t', quote=F)
head(df.v1v2)
nrow(df.v1v2)

#V1V3
df.v1v3 <- data.frame('TaxonName_v1v3'=data.v1v3$taxon_data$name, 
                      'TaxonCount_v1v3'=n_obs(data.v1v3),
                      'TaxonHierarchy_v1v3'=hierarchies(data.v1v3),
                      row.names=data.v1v3$taxon_ids,
                      stringsAsFactors=FALSE)
write.table(df.v1v3, file='v1v3_hierachy.tsv', sep='\t', quote=F)
head(df.v1v3)
nrow(df.v1v3)

#V1V9
df.v1v9 <- data.frame('TaxonName_v1v9'=data.v1v9$taxon_data$name, 
                      'TaxonCount_v1v9'=n_obs(data.v1v9),
                      'TaxonHierarchy_v1v9'=hierarchies(data.v1v9),
                      row.names=data.v1v9$taxon_ids,
                      stringsAsFactors=FALSE)
write.table(df.v1v9, file='v1v9_hierachy.tsv', sep='\t', quote=F)
head(df.v1v9)
nrow(df.v1v9)

#V3V5
df.v3v5 <- data.frame('TaxonName_v3v5'=data.v3v5$taxon_data$name, 
                      'TaxonCount_v3v5'=n_obs(data.v3v5),
                      'TaxonHierarchy_v3v5'=hierarchies(data.v3v5),
                      row.names=data.v3v5$taxon_ids,
                      stringsAsFactors=FALSE)
write.table(df.v3v5, file='v3v5_hierachy.tsv', sep='\t', quote=F)
head(df.v3v5)
nrow(df.v3v5)

#V4
df.v4 <- data.frame('TaxonName_v4'=data.v4$taxon_data$name, 
                    'TaxonCount_v4'=n_obs(data.v4),
                    'TaxonHierarchy_v4'=hierarchies(data.v4),
                    row.names=data.v4$taxon_ids,
                    stringsAsFactors=FALSE)
write.table(df.v4, file='v4_hierachy.tsv', sep='\t', quote=F)
head(df.v4)
nrow(df.v4)

#V6V9
df.v6v9 <- data.frame('TaxonName_v6v9'=data.v6v9$taxon_data$name, 
                      'TaxonCount_v6v9'=n_obs(data.v6v9),
                      'TaxonHierarchy_v6v9'=hierarchies(data.v6v9),
                      row.names=data.v6v9$taxon_ids,
                      stringsAsFactors=FALSE)
write.table(df.v6v9, file='v6v9_hierachy.tsv', sep='\t', quote=F)
head(df.v6v9)
nrow(df.v6v9)

### Merge the output...
df <- merge(x=df.master, y=df.v1v2, by.x='TaxonName', by.y='TaxonName_v1v2', all.x=TRUE, sort=FALSE)
head(df)
nrow(df) - nrow(df.master)

# checking
df[is.na(df$TaxonName),]
df[is.na(df$TaxonCount),]
df[is.na(df$TaxonCount_v1v2),]

df[duplicated(df$TaxonName),]

df.master[duplicated(df.master$TaxonName),]
df.v1v2[duplicated(df.v1v2$TaxonName_v1v2),]
df.master[duplicated(df.master$TaxonName),]
df.master[duplicated(df.master$TaxonName),]
df.master[duplicated(df.master$TaxonName),]
df.master[duplicated(df.master$TaxonName),]
