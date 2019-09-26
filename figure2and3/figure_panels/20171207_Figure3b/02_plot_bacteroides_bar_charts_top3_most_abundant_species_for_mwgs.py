
# coding: utf-8

# In[12]:


# there's a bug in the code that generates the abundance tables, the scientific notation e-XX is converted to e_XX
# this needs to be fixed for numbers to be recognized. %%!


# In[13]:


get_ipython().run_cell_magic('bash', '', "\ncat 00_data/distance_1.0_bacteroides_species_abundance.tsv | sed 's/e_/e-/g' | \\\nsed 's/Horse-/Horse_/g' > 02_bacteroides_species_abundance_modified.tsv")


# In[14]:


import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from collections import Counter
import re


# In[15]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[16]:


# Fetch the corrected table and drop the non-numeric columns
df = pd.read_table('02_bacteroides_species_abundance_modified.tsv', index_col=0)
df.head()


# In[17]:


# reset index to contain unique identifiers for the OTUs that were missing from mWGS data and drop columns
new_index = []
for row in df.iterrows():
    if row[0] == 'Missing':
        if row[1]['V1V3'] != 'Missing':
            new_index.append('V1V3_' + row[1]['V1V3'])
        else:
            new_index.append('V1V9_' + row[1]['V1V9'])
    else:
        new_index.append(row[0])
df.index = new_index
df = df.drop(['V1V3', 'V1V9'], axis=1)


# In[18]:


# How to normalize row values in a pandas dataframe 
df = df.div(df.sum(axis=0), axis=1)


# In[19]:


# Select the top three most abundant taxa in each sample for mWGS data
species = []

for name in ['Breezer', 'Commencal', 'IronHorse', 'Scott']:
    mwgs = df[name +'_MGWS'].order(ascending=False)
    m = list(mwgs[:3].index)
    species.extend(m)
    
for k, v in Counter(species).items():
    print k, v


# In[20]:


# Subset only those taxa that are in the top three most abundant for one or more sample as quantified by mWGS
dfm = df[df.index.isin(species)]
dfm.index.name = 'MWGS'
dfm.reset_index(inplace=True)


# In[21]:


dfm.head()


# In[22]:


dfm = pd.melt(dfm, id_vars=['MWGS'])


# In[23]:


dfm.head()


# In[24]:


# clean up table
sample = [x.split('_')[0] for x in dfm['variable']]
platform = [x.split('_')[1] for x in dfm['variable']]
platform = ['MWGS' if x == 'MGWS' else x for x in platform ]
species = [re.sub('Bacteroides', 'B.', x) for x in dfm['MWGS']]
platform
dfm['platform'] = platform
dfm['sample'] = sample
dfm['species'] = species
dfm = dfm.drop(['variable', 'MWGS'], axis=1)
dfm.head()


# In[25]:


get_ipython().run_line_magic('Rpush', 'dfm')


# In[26]:


get_ipython().run_cell_magic('R', '', "library('ggplot2')")


# In[27]:


get_ipython().run_cell_magic('R', '', '\npl <- ggplot(dfm, aes(x=species, y=value, fill=platform)) + geom_bar(stat=\'identity\', position=\'dodge\')\npl <- pl + facet_wrap(~sample, ncol=1)\npl <- pl + theme_bw()\npl <- pl + theme(axis.text.x = element_text(angle = 90, hjust = 1))\npl <- pl + ylab(\'Proportion of Bacteroides reads\') + xlab("")\nplot(pl)')


# In[28]:


get_ipython().run_cell_magic('R', '', "pdf('02_BarChart_proportion_of_bacteroides_reads.pdf', width=5, height=6)\nplot(pl)\ndev.off()")

