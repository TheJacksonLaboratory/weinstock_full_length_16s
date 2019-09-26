
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


# Read in the table... NA values occur for gener that were absent from V1V3 or V1V9 data, but not both. 
df = pd.read_table('00_data/ampds_distance_1.0_genus_correlations.txt', index_col=0)
# convert the NA values to 0
df = df.fillna(value=0)
df.head()


# In[3]:


# Normalize the data
df = df.div(df.sum(axis=0)).multiply(100)
df.head()


# In[4]:


# Fetch the V1-V3 reads
v1v3 = df.filter(like='V1_V3').loc['Bacteroides',]
v1v3.index = [x.split('_')[0] for x in v1v3.index]
v1v3.name = 'V1_V3'
v1v3


# In[5]:


v1v9 = df.filter(like='V1_V9').loc['Bacteroides',]
v1v9.index = [x.split('_')[0] for x in v1v9.index]
v1v9.name = 'V1_V9'
v1v9


# In[6]:


df_b = pd.DataFrame([v1v9, v1v3]).T
df_b.index = ['Scott', 'IronHorse', 'Commencal', 'Breezer']
df_b


# In[7]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[8]:


from rpy2.robjects import pandas2ri
pandas2ri.activate()


# In[9]:


get_ipython().run_line_magic('Rpush', 'df_b')


# In[10]:


get_ipython().run_cell_magic('R', '', "library('ggplot2')")


# In[12]:


get_ipython().run_cell_magic('R', '', "pl <- ggplot(df_b, aes(x=V1_V9, y=V1_V3, label=rownames(df_b))) + geom_point() + geom_abline(intercept=0, lty=2)\npl <- pl + geom_text(hjust=0, vjust=0)\npl <- pl + xlim(0, 100) + ylim(0, 100)\npl <- pl + xlab('V1-V9 Bacteroides Abundance (%)')\npl <- pl + ylab('V1-V3 Bacteroides Abundance (%)')\npl <- pl + theme_bw()\npl")

