###############################################################################
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.  
###############################################################################
"""
===========================================================================
pipeline_blast_16s.py

A pipeline to determine the taxonomic identity of full length 16S sequences
===========================================================================

Jethro Johnson

"""
# load modules
from ruffus import *

import sys
import os
import re
import glob
import math
import pickle
import gzip
import sqlite3
import collections
import shutil
import pandas as pd
import numpy as np
import random
import logging as L

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
pandas2ri.activate()

import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools

import PipelineTools as P
import PipelineBlast16S as PB

# Pipeline configuration
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS


def connect():
    """
    Connect to default sqlite database, returning a database handle
    """
    #dbh = sqlite3.connect(PARAMS['database'])
    dbh = sqlite3.connect('csvdb')
    return dbh


###############################################################################
# Pipeline
###############################################################################
# @transform('raindance_otus.xlsx', suffix('.xlsx'), '.load')
# def loadOriginalOTUClassifications(infile, outfile):

#     df = pd.read_excel(infile, index_col=0)
#     df.to_sql(name='original_classifications',
#               con=connect(),
#               index=True,
#               if_exists='replace')
#     open(outfile, 'w').close()

@follows(mkdir('blast_matches.dir'))
@transform('input.dir/*fast*',
           regex('.*/(.+).fast.*'),
           r'blast_matches.dir/\1_blast.nt')
def blastToFullNucleotideDatabase(infile, outfile):
    '''Blast fasta sequences against full nucleotide database'''

    if infile.endswith('.gz'):
         x = 'zcat'
    else:
         x = 'cat'
    statement = ("blastn"
                 "  -db {database_blast_db}"
                 "  -query <({x} {infile})"
                 "  -outfmt \"6 qseqid sseqid sscinames pident length evalue"
                 "            bitscore qstart qend qlen sstart send\""
                 " {blastn_options}"
                 " -out {outfile}"
                 " &> {outfile}.log")
    cluster_options='-l walltime=24:00:00,mem=10Gb'
    P.run()


@jobs_limit(4)
@transform(blastToFullNucleotideDatabase,
           regex('(.+)_blast.nt'),
           r'\1_gi_numbers.tsv')
def extractGINumbers(infile, outfile):
    '''Extract GI numbers from blast table,
    taking only the hit with lowest evalue for each GI'''
    tmpf = P.getTempFilename('.')
    statement = ("cat {infile} |"
                 " sort -k6,6 -k7nr,7nr" # sort on evalue, then bitscore
                 " > {tmpf}")
    to_cluster = False
    P.run()
    seq_ids = set()
    with IOTools.openFile(outfile, 'w') as outf:
        for line in IOTools.openFile(tmpf):
            seq_id, match_id = line.split()[:2]
            # only take the best blast match
            if seq_id in seq_ids:
                continue
            else:
                gi = match_id.split('|')[1]
                outf.write(gi + '\t' + seq_id + '\n')
                seq_ids.add(seq_id)
    os.unlink(tmpf)
    

@transform(extractGINumbers,
           suffix('_numbers.tsv'),
           add_inputs(PARAMS['blastn_gi2taxid']),
           '2taxid.tsv')
def fetchGITaxIDs(infiles, outfile):
    '''Fetch the taxid for each gi number using the gi_taxid_nucl.dmp.
    Output gi number, tax id. This assumes the gi number is present in
    the gi2taxid mapping file, and present only once. 
    '''
    # infile is gi #, seq_id, tax_map_file is gi, taxid
    infile, tax_map_file = infiles
    outfile = os.path.abspath(outfile)
    
    if os.path.exists(outfile):
        os.unlink(outfile)

    statement = (" cat {infile} |"
                 "  cut -f1 |"
                 "  while read X;"
                 "   do grep -m 1 \"^$X\\s\" {tax_map_file}"
                 "  >> {outfile};"
                 " done;")
    cluster_options='-l walltime=04:00:00,mem=10Gb'
    P.run()


@follows(mkdir('read_taxonomy.dir'))
@jobs_limit(4)
@transform(fetchGITaxIDs,
           regex('.+/(.+)_gi2taxid.tsv'),
           add_inputs(PARAMS['blastn_names_dmp'], PARAMS['blastn_nodes_dmp']),
           r'read_taxonomy.dir/\1_taxonomy.tsv')
def fetchGITaxonomy(infiles, outfile):
    '''Parse the ncbi taxonomy tree to retrieve the full taxonomy for
    each unmapped read.'''
    infile, names_dmp, nodes_dmp = infiles
    
    PB.fetchTaxonomy(infile, names_dmp, nodes_dmp, outfile)
    

@merge(fetchGITaxonomy, 'read_taxonomy.tsv.gz')
def combineTaxonomy(infiles, outfile):
    '''Combine the taxonomic information for each sample'''

    infiles = ' '.join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --log={outfile}.log"
                 " --cat=SampleID"
                 " --regex-filename='.+/(.+)_taxonomy.tsv'"
                 " {infiles} |"
                 " gzip > {outfile}")
    to_cluster=False
    P.run()
    
@transform(fetchGITaxonomy, suffix('.tsv'), '.load')
def loadTaxonomy(infile, outfile):
    P.load(infile, outfile)
           

@jobs_limit(1)
@files(fetchGITaxonomy, 'species_summary.pdf')
def plotUnmappedReadSummary(infiles, outfile):
    '''Plot a stacked bar chart of unmapped read classification'''

    R('''rm(list=ls())''')
    
    def _collapseDataFrame(df, limit=5, outf=None):
        '''Convert to %, collapse anything less than 5%,
        return the melted dataframe'''
        df.fillna(0, inplace=True)
        # convert to %
        df = df.apply(lambda x: x/float(x.sum())*100, axis=0)
        if outf:
            df.to_csv(IOTools.openFile(outf, 'w'), sep='\t')
        
        # collapse anything less than limit
        new_index = []
        for row in df.iterrows():
            if max(row[1]) < limit:
                new_index.append('Other')
            else:
                new_index.append(row[0])
        df.index = new_index
        df = df.groupby(level=0).sum()
        # melt
        df.index.name = 'Taxon'
        df = df.reset_index()
        df = pd.melt(df,
                     id_vars=['Taxon',],
                     var_name='SampleID',
                     value_name='Count')

        return df

    first = True
    for infile in infiles:
        sampleID = os.path.basename(infile).split('_')[0]
        df = pd.read_table(infile)
        df_s = pd.DataFrame(df['S'].value_counts(), columns=[sampleID,])
        df_f = pd.DataFrame(df['F'].value_counts(), columns=[sampleID,])
        if first:
            first = False
            df_s_out = df_s
            df_f_out = df_f
            continue
        else:
            df_s_out = df_s_out.join(df_s, how='outer')
            df_f_out = df_f_out.join(df_f, how='outer')

    # collapse dataframes to %...
    # ...for species...
    sp_tab = P.snip(outfile, '.pdf') + '.tsv.gz'
    df_s = _collapseDataFrame(df_s_out, outf=sp_tab)
    # ...and family
    f_tab = P.snip(outfile, '_species_summary.pdf') + '_family_summary.tsv.gz'
    df_f = _collapseDataFrame(df_f_out, outf=f_tab)

    # plot
    f_out = P.snip(f_tab, '.tsv.gz') + '.pdf'
    R.assign('df_s', df_s)
    R.assign('df_f', df_f)

    R('''
    require('ggplot2')
    s <- ggplot(df_s, aes(x=SampleID, y=Count, fill=Taxon)) + geom_bar(stat='identity')
    s <- s + theme_bw() + ylab('Percent of OTUs')
    s <- s + scale_fill_brewer(type='div', palette='Spectral')

    f <- ggplot(df_f, aes(x=SampleID, y=Count, fill=Taxon)) + geom_bar(stat='identity')
    f <- f + theme_bw() + ylab('Percent of OTUs')
    f <- f + scale_fill_brewer(type='div', palette='Spectral')

    pdf('{outfile}')
    plot(s)
    dev.off()

    pdf('{f_out}')
    plot(f)
    dev.off()
    rm(list=ls())'''.format(**locals()))

###############################################################################
# Pipeline options
###############################################################################
@follows(loadTaxonomy)
def full():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
