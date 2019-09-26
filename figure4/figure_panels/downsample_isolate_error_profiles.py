
# coding: utf-8

# # Randomly downsample the distinct error profiles

# In[9]:


import random
import pandas as pd


# In[3]:


df = pd.read_table('isolate_error_profiles.tsv', index_col=0, sep='\t')
df = df.transpose()
df.head()


# In[5]:


df['OTU'] = ['_'.join(x.split('_')[:2]) for x in df.index]
df.head()


# In[20]:


to_keep = []
for otu in set(df['OTU']):
    samples = df[df['OTU'] == otu].index.tolist()
    
    if len(samples) < 5:
        to_keep.extend(samples)
    else:
        to_keep.extend(random.sample(samples, 5))

to_keep


# In[25]:


df_out = df[df.index.isin(to_keep)].transpose()
df_out.drop('OTU', axis=0, inplace=True)
df_out.head()


# In[26]:


df_out.to_csv('downsampled_isolate_error_profiles.tsv', sep='\t')

