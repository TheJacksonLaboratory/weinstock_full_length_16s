
# coding: utf-8

# In[1]:


# Parse the pacbio substitution matrix to get the per-base substitution percentages
# Load the per-base substitution percentages from illumina data
# create a combined dataframe for plotting in R


# In[2]:


import pandas as pd


# In[5]:


# Fetch the illumina substitution data
df_il = pd.read_table('01_illumina_substitutions.tsv', sep='\t')
df_il.head()


# In[9]:


# Fetch the pacbio substitution matrix
df_pb = pd.read_table('00_data/Substitution_matrix.csv', sep=',', index_col=0, header=None)
df_pb.head()


# In[10]:


# Parse the substitution matrix into the desired format
df_pb = pd.DataFrame(df_pb.sum(axis=1))
df_pb.columns = ['Substitutions',]
df_pb.reset_index(inplace=True)

df_pb.columns = ['Position', 'Substitutions']
df_pb['Panel'] = ['upper',]*len(df_pb.index)

# Converting the substitution rate to percentages
df_pb['Substitutions'] = df_pb['Substitutions']/25676.0

df_pb.head()


# In[11]:


df = pd.concat([df_pb, df_il])
print len(df_pb.index)
print len(df_il.index)
print len(df.index)


# In[13]:


df.to_csv('02_pb_il_substitution_matrix.tsv', sep='\t')

