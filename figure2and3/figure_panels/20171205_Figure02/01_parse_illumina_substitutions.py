
# coding: utf-8

# In[1]:


import re


# In[94]:


infile = "00_data/diaz.hiseq.sort.pileup.proc"


# In[62]:


# Don't try to get into screening out indels... it's complicated...
# ... in indels greater than 1, how to account for the exact number?
d = ".$,$.,.....,,,.....,..,.,..,..,..,,.,,.,,,,,,,.,,.,,,.,-1a.,,,,,,,.,.,,..,,...,,,,,..,,.,,,..,,..,..,,.,,...,.,,,.,,^].^]."
s = ".TTT......T..TT..TT.TTTTTT.TTT,,TT...T.T.TTTTT....T.TT.t,T....T.TTT.^].-1a."
i = ".$,$.,.....,,,.....,..,.,..,..,..,,.,,.,,,,,,,.,,.,,,.,+1a.,,,,,,,.,.,,..,,...,,,,,..,,.,,,..,,..,..,,.,,...,.,,,.,,^].^]."

re.search('(?<![\-\+][1-9])[ACTGNactgn]', i)


# In[5]:


# Instead make the assumption that $5 is number of substitutions
n = 0
for line in open(pileup):
    if re.search('[ACTGNactgn]', line):
        n += 1
        print line
        
    if n >= 10:
        break


# In[93]:


print len(re.findall('G', '.GGG...G.G.G.GGGGGGG.GGGGGG..GGG,gGG.GGGG.GGGGGG.G.GGGG.g,G..G.G.GGG.'.upper()))
print len(re.findall('A', '.AAA.......A..AA..A.AAAAAA.$.AAA,,AA...A.A.AAAAA....A.AA.a,A....A.AAA.'.upper()))
print len(re.findall('T', '.TTT......T..TT..TT.TTTTTT.TTT,,TT...T.T.TTTTT....T.TT.t,T....T.TTT..'))


# In[96]:


# $5 is 0 when there are deletions but no substitutions
for line in open(infile):
    if re.search('-[0-9]+[ACGTNacgtn]+', line.split('\t')[-1]): # , line.split('\t')[-1]
        print line
        break


# In[97]:


# There are no insertions
for line in open(infile):
    if re.search('\+[0-9]+[ACGTNacgtn]+', line.split('\t')[-1]): # , line.split('\t')[-1]
        print line
        break


# In[100]:


outfile = '01_illumina_substitutions.tsv'

with open(outfile, 'w') as outf:
    outf.write('Position\tSubstitutions\tPanel\n')
    for line in open(infile):
        line = line.split('\t')
        position = line[0]
        coverage = line[1]
        substitutions = line[4]

        percent_substitutions = float(substitutions)/float(coverage)
        
        outf.write('\t'.join([position, str(percent_substitutions), 'lower']) + '\n')
        

