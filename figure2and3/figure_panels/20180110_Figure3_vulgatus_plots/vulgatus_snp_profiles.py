
# coding: utf-8

# In[1]:


import pandas as pd
from rpy2.robjects import pandas2ri
pandas2ri.activate()


# In[2]:


get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[3]:


# read the (slightly modified) data into python...
df = pd.read_excel('Bvulgatis_rRNAgenes_Variants_edited.xlsx', sheetname='Sheet2')
df.head()


# In[4]:


df_1 = df[df['Strain'] == 'A']

dict1 = {}

for row in df_1.iterrows():
    p = str(row[1]['Reference Position'])
    pptn = float(row[1]['Pptn'])
    dict1[p] = pptn

for i in range(1, 1573):
    i = str(i)
    if i not in dict1.keys():
        dict1[i] = 0.0


# In[5]:


# Parse the data into the form used for the pb and il data... 

df_2 = df[df['Strain'] == 'B']



dict2 = {}
for row in df_2.iterrows():
    p = str(row[1]['Reference Position'])
    pptn = float(row[1]['Pptn'])
    dict2[p] = pptn
    
for i in range(1, 1573):
    i = str(i)
    
    if i not in dict2.keys():
        dict2[i] = 0.0


# In[6]:


# Create dataframes
mg = pd.DataFrame(pd.Series(dict1))
mg.index = [int(x) for x in mg.index]
mg = mg.sort_index()
mg.reset_index(inplace=True)
mg.columns = ['Position', 'Substitutions']
mg['Panel'] = 'CP000139.1'
mg
sum(mg['Substitutions'])


# In[8]:


mg.columns = ['BasePosition', 'bacteroides_vulgatus', 'panel']
mg.head()


# In[9]:


mg.to_csv('cp000139_snp_profile.tsv', sep='\t')


# In[10]:


# Create dataframes
sk = pd.DataFrame(pd.Series(dict2))
sk.index = [int(x) for x in sk.index]
sk = sk.sort_index()
sk.reset_index(inplace=True)
sk.columns = ['Position', 'Substitutions']
sk['Panel'] = 'CP013020.1'
sk['Substitutions'] = -1*sk['Substitutions']
sk.head()


# In[11]:


sk.columns = ['BasePosition', 'bacteroides_vulgatus', 'panel']
sk.head()


# In[12]:


sk.to_csv('cp013020_snp_profile.tsv', sep='\t')


# In[55]:


df = pd.concat([mg, sk])
print len(mg.index)
print len(sk.index)
print len(df.index)


# In[56]:


get_ipython().run_line_magic('Rpush', 'df')


# In[57]:


get_ipython().run_cell_magic('R', '', "library('ggplot2')")


# In[64]:


get_ipython().run_cell_magic('R', '', "df$Substitutions = df$Substitutions/100\n\n###################################################################################################\n## Plot the panel figure\nplt <- ggplot(df, aes(x=Position, y=Substitutions)) + geom_line()\nplt <- plt + theme_bw()\nplt <- plt + theme(panel.grid = element_blank())\nplt <- plt + ggtitle('Bacteroides vulgatus') + xlab('') + ylab('')\nplt")


# In[62]:


df.to_csv('04_mg1655_vs_sakai_counts.tsv', sep='\t')


# In[65]:


get_ipython().run_cell_magic('R', '', "plt <- plt + ylim(c(-1.2, 1.2))\n\n# parameters\n# region title offset\nn = 0.01\n# region title alpha\nra=0.8\n# region title colour\nrc='dark grey'\n# panel alpha\npa=0.2\n#panel colour\npc='grey'\n#amplicon colour\nac = 'dark red'\n\n\n\n# add colour panels\nplt <- plt + annotate('rect', xmin=9, xmax=104, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=103, xmax=249, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=341, xmax=515, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=514, xmax=685, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=778, xmax=909, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=968, xmax=1061, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=1098, xmax=1177, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=1222, xmax=1370, ymin=-1, ymax=1, fill=pc, alpha=pa)\nplt <- plt + annotate('rect', xmin=1377, xmax=1492, ymin=-1, ymax=1, fill=pc, alpha=pa)\n\nplt <- plt + xlab('Position along 16S gene') + ylab('Percent Substitutions')\n\n# add individual bars\n#V1\nplt <- plt + annotate('segment', x=9, xend=104, y=1, yend=1, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=9 + (104-9)/2, y=1+n, label='V1', cex=3)\n#V2\nplt <- plt + annotate('segment', x=103, xend=249, y=1+n, yend=1+n, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=103 + (249-103)/2, y=1+n+n, label='V2', cex=3)\n#V3\nplt <- plt + annotate('segment', x=341, xend=515, y=1, yend=1, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=341 + (515-341)/2, y=1+n, label='V3', cex=3)\n#V4\nplt <- plt + annotate('segment', x=514, xend=685, y=1+n, yend=1+n, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=514 + (685-514)/2, y=1+n+n, label='V4', cex=3)\n#V5\nplt <- plt + annotate('segment', x=778, xend=909, y=1, yend=1, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=778 + (909-778)/2, y=1+n, label='V5', cex=3)\n#V6\nplt <- plt + annotate('segment', x=968, xend=1061, y=1+n, yend=1+n, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=968 + (1061-968)/2, y=1+n+n, label='V6', cex=3)\n#V7\nplt <- plt + annotate('segment', x=1098, xend=1177, y=1, yend=1, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=1098 + (1177-1098)/2, y=1+n, label='V7', cex=3)\n#V8\nplt <- plt + annotate('segment', x=1222, xend=1370, y=1+n, yend=1+n, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=1222 + (1370-1222)/2, y=1+n+n, label='V8', cex=3)\n#V9\nplt <- plt + annotate('segment', x=1377, xend=1495, y=1, yend=1, colour=rc, alpha=ra)\nplt <- plt + annotate('text', x=1377 + (1492-1377)/2, y=1+n, label='V9', cex=3)\nplt\n")


# In[66]:


get_ipython().run_cell_magic('R', '', "pdf('mg1655_vs_sakai_substitutions.pdf', height=3, width=3.5)\nplot(plt)\ndev.off()")

