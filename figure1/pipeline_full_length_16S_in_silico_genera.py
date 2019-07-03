
"""

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
import pysam

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects import r as R

import PipelineTools as P

import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator
import CGAT.BamTools as BamTools

import PipelineFullLength16S as P16S

# Pipeline configuration
PARAMS = P.getParameters(["%s/pipeline.ini" % os.path.splitext(__file__)[0],
                          "pipeline.ini"],)

###############################################################################
# Utility functions
###############################################################################
def connect():
    """
    Connect to default sqlite database, returning a database handle
    """
    #dbh = sqlite3.connect(PARAMS['database'])
    dbh = sqlite3.connect('csvdb')
    return dbh

# to be exported
def multiSplitFasta(infile, outfile_stub, n_outfiles, random=False):
    '''split fasta into n smaller files'''
    
    outfile_handles = []
    outfile_names = []
    
    # create list of upper bounds for intervals
    intervals = []
    lower = 0
    for i in range(n_outfiles):
        upper = lower + 1.0/n_outfiles
        intervals.append(upper)
        # add an outfile handle to list of outfile handles
        outf = outfile_stub + '_' + str(i).zfill(2) + '.fasta.gz'
        outfile_names.append(outf)
        outfile_handles.append(gzip.GzipFile(outf, 'w'))
        lower = upper

    if not random:
        n_seq = len([n.title for n in FastaIterator.FastaIterator(
            IOTools.openFile(infile))])
        intervals = [int(math.ceil(x*n_seq)) for x in intervals]

    def _get_number(n):
        if random:
            n = random.random()
        return n

    # iterate through reads in file and write them to the outfiles
    seq_num = 0
    for seq in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        seq_num += 1
        num = _get_number(seq_num)

        for i in range(len(intervals)):
            if num <= intervals[i]:
                outfile_handles[i].write('\n'.join(['>' + seq.title,
                                                        seq.sequence]) + '\n')
                break
            else:
                continue


###############################################################################
# Metatask: QC OF REFERENCE DATABASE
###############################################################################
## Task: QC of the mock community database. 
###############################################################################
### Subtask: Summarize sequence lengths in mock community fasta
###############################################################################
@follows(mkdir('mock_ref_db_qc.dir'))
@transform(PARAMS['data_mock_ref_db'],
           regex('.+/(.+).fa'),
           r'mock_ref_db_qc.dir/\1_summary.tsv')
def summarizeRefDB(infile, outfile):
    """
    Run fasta2table to get summary of sequence read lengths and composition of
    mock community reference database. 
    """

    statement = ("python {scriptsdir}/fasta2table.py"
                 " --section=length,na,gaps"
                 " --split-fasta-identifier"
                 " --log={outfile}.log"
                 " < {infile} > {outfile}")
    P.run()


@transform(summarizeRefDB,
           regex('.+/(.+).tsv'),
           r'\1.load')
def loadRefDBSummary(infile, outfile):

    table_name = P.snip(infile, '.tsv',  strip_path=True)
    table_name = re.sub('-', '_', table_name)
    to_cluster=False
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 " > {outfile}")
    P.run()
    

@jobs_limit(1, 'RGlobal')
@transform(summarizeRefDB, suffix('.tsv'), '.png')
def plotRefDBSummary(infile, outfile):
    """plot histograms of mock community reference database sequence
    characteristics"""

    df = pd.read_table(infile, sep='\t', index_col='id')
    df = pd.melt(df)
    R.assign('df', df)

    R('''
      require('ggplot2')
      hist <- ggplot(df, aes(x=value))
      hist <- hist  + geom_histogram(fill=c('darkblue'))
      hist <- hist  + facet_wrap(~variable, scales='free')
      hist <- hist + theme_bw() + theme(panel.grid=element_blank())
      png('{outfile}')
      plot(hist)
      dev.off()
      rm(list=ls())'''.format(**locals()))


@transform(PARAMS['data_mock_ref_db_taxonomy'],
           regex('.+/(.+).tsv'),
           r'\1.load')
def loadRefDBTaxSummary(infile, outfile):
    '''load summary of reference database taxonomy'''

    table_name = P.snip(infile, '.tsv',  strip_path=True)

    df = pd.read_table(infile, sep='\t', index_col='StrainID')

    df.to_sql(name=re.sub('-', '_', table_name),
              con=connect(),
              index=True,
              if_exists='replace')

    open(outfile, 'w').close()


@transform(loadRefDBTaxSummary,
           regex('(.+).load'),
           add_inputs(loadRefDBSummary),
           r'mock_ref_db_qc.dir/\1.xls')
def fetchRefDBLengthTable(infiles, outfile):
    '''Get a taxonomic table containing ref sequence lengths'''
    tab_sum, tab_tax = infiles
    tab_sum = re.sub('-', '_', P.snip(tab_sum, '.load', strip_path=True))
    tab_tax = re.sub('-', '_', P.snip(tab_tax, '.load', strip_path=True))
    
    statement = ("SELECT b.StrainID,b.StrainName,b.phylum,b.class,"
                 "b.`order`,b.family,b.genus,a.length"
                 " FROM {tab_sum} AS b"
                 " INNER JOIN {tab_tax} AS a"
                 " ON a.id == b.StrainID ORDER BY a.length".format(**locals()))
    df = pd.read_sql_query(sql=statement, con=connect(), index_col='StrainID')
    df.to_excel(outfile)


###############################################################################
### Subtask: Count the homopolymer runs in reference database sequences
###############################################################################
@follows(mkdir('mock_ref_db_qc.dir'))
@subdivide(PARAMS['data_mock_ref_db'],
           regex('.+.fa'),
           r'mock_ref_db_qc.dir/*_homopolymers.txt')
def findReferenceHomopolymerRuns(infile, outfiles):
    '''Output a table of homopolymer runs for each of the sequences in
    the reference database. Position is zero-based.
    '''

    # write a homopolymer iterator
    def homopolymer_iterator(seq):
        last = seq[0]
        loc = 0
        homopolymer = []
        for i, n in enumerate(seq.upper()):
            if n != last:
                yield (loc, homopolymer)
                homopolymer = []
                last = n
                loc = i
            homopolymer.append(n)

    out_dir = 'mock_ref_db_qc.dir'
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        strain = fasta.title.split()[0]
        outfile = os.path.join(out_dir, strain + '_homopolymers.txt')
        # write summary table
        with open(outfile, 'w') as outf:
            outf.write('Position\tNucleotide\tHomopolymerLength\n')
            for homopolymer in homopolymer_iterator(fasta.sequence):
                if len(homopolymer[1]) >= int(PARAMS['cross_match_homopolymer_length']):
                    outf.write('\t'.join(map(str, [homopolymer[0], 
                                                   homopolymer[1][0], 
                                                   len(homopolymer[1])])) + '\n')


@collate(findReferenceHomopolymerRuns,
         regex('(.+)/.+_homopolymers.txt'),
         r'\1/reference_homopolymers.txt')
def collateHomopolymerRuns(infiles, outfile):
    '''Create a stacked table of homopolymer run information '''
    to_cluster = False
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 "  --cat=Strain"
                 "  --regex-filename='.+/(.+)_homopolymers.txt'"
                 "  --log={outfile}.log"
                 " {infiles}"
                 " > {outfile}")
    P.run()


@transform(collateHomopolymerRuns,
           regex('.+/(.+).txt'),
           r'\1.load')
def loadHomopolymerRuns(infile, outfile):
    
    table_name = P.snip(os.path.basename(infile), '.txt')

    to_cluster=False
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 " > {outfile}")
    P.run()


###############################################################################
# Metatask: IN SILICO AMPLICON COMPARISON
###############################################################################
# Amplicons with which to perform in silico analysis
AMPLICON_TARGETS = {}
header = True
for line in IOTools.openFile(PARAMS['in_silico_amplicons']).readlines():
    if header:
        assert line == 'Region,Forward,ForwardName,Reverse,ReverseName,ExpectedSize\n', \
            'Unexpected file format: {}'.format(line)
        header = False
        continue
    AMPLICON_TARGETS[line.split(',')[0]] = line.split(',')[1:]


###############################################################################
## Task: Create database specific fasta and taxonomy files
###############################################################################
### Subtask: GreenGenes
###############################################################################
@follows(mkdir('greengenes_insilico_16S.dir'))
@files(PARAMS['green_genes_gg_taxonomy'],
       'greengenes_insilico_16S.dir/greengenes_reference_taxonomy.txt.gz')
def parseGreenGenesTaxonomy(infile, outfile):
    '''Convert greengenes database into a format suitable for in silico
    analysis. Outfiles must have suffix *reference* and prefix to match
    parent directory.
    '''
    
    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')
    
    statement = ("cat {infile} |"
                 " python {tax_parser} -t '; ' -g -s -m 6 -x ID -r"
                 " --log={outfile}.log |"
                 " gzip > {outfile}")
    P.run()


@transform(PARAMS['green_genes_gg_fasta'],
           regex('.+/.+.fasta'),
           add_inputs(parseGreenGenesTaxonomy),
           'greengenes_insilico_16S.dir/greengenes_reference.fasta.gz')
def parseGreenGenesFasta(infiles, outfile):
    '''Parse the greengenes fasta file so that it contains only sequences
    that also appear in the filtered taxonomy file. An 'ID' prefix is
    also added to the sequence IDs. 
    '''

    fasta_file, tax_file = infiles

    # get a list of the sequence IDs to be retained
    IDs = []
    for line in IOTools.openFile(tax_file):
        IDs.append(line.split()[0])

    # parse the fasta file
    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
            if 'ID' + fasta.title in IDs:
                outf.write('>ID' + fasta.title + '\n' + fasta.sequence + '\n')
            else:
                continue


@transform(parseGreenGenesFasta,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/subsampled.fasta.gz')
def subSampleGreenGenes(infile, outfile):
    '''Downsample green genes database to specified number of sequences'''

    infile = os.path.abspath(infile)

    # number of genes in greengenes database
    n = 0
    for i in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        n += 1
    L.info('There are %i sequences in %s' % (n, infile))

    if PARAMS['green_genes_n_seq']:
        pptn = str(float(PARAMS['green_genes_n_seq'])/n)

        statement = ("zcat {infile} |"
                     " python {scriptsdir}/fasta2fasta.py"
                     " --method=sample"
                     " --sample-proportion={pptn}"
                     " --random-seed={green_genes_seed}"
                     " --log={outfile}.log |"
                     " gzip > {outfile}")
    else:
        to_cluster = False
        statement = ("ln -s {infile} {outfile}")

    P.run()


@follows(subSampleGreenGenes)
@split(None, [os.path.join('greengenes_insilico_16S.dir', x, x + '_primers.txt')
              for x in AMPLICON_TARGETS.keys()])
def makeGreenGenesDirs(infile, outfiles):
    '''Create directories for generating in silico amplicons'''
    P16S.makeInSilicoDirs(outfiles, AMPLICON_TARGETS)


###############################################################################
### Subtask: mothur's greengenes database
###############################################################################
@follows(mkdir('mothur_greengenes_insilico_16S.dir'))
@files([PARAMS['mothur_green_genes_fasta'], PARAMS['mothur_green_genes_taxonomy']],
       ['mothur_greengenes_insilico_16S.dir/mothur_greengenes_reference.fasta.gz',
        'mothur_greengenes_insilico_16S.dir/mothur_greengenes_reference_taxonomy.txt.gz'])
def parseMothurGreenGenesDatabase(infiles, outfiles):
    '''Convert greengenes database into a format suitable for in silico
    analysis. Outfiles must have suffix *reference* and prefix to match
    parent directory.
    '''
    fasta_in, tax_in = infiles
    fasta_out, tax_out = outfiles
    
    P16S.FormatMothurGreenGenes(fasta_in, tax_in, fasta_out, tax_out).run()
    

@transform(parseMothurGreenGenesDatabase,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/subsampled.fasta.gz')
def subSampleMothurGreenGenes(infiles, outfile):
    '''Downsample green genes database to specified number of sequences'''

    infile, tax_file = infiles
    infile = os.path.abspath(infile)

    # number of genes in greengenes database
    n = 0
    for i in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        n += 1
    L.info('There are %i sequences in %s' % (n, infile))

    if PARAMS['mothur_green_genes_n_seq']:
        pptn = str(float(PARAMS['mothur_green_genes_n_seq'])/n)

        statement = ("zcat {infile} |"
                     " python {scriptsdir}/fasta2fasta.py"
                     " --method=sample"
                     " --sample-proportion={pptn}"
                     " --random-seed={mothur_green_genes_seed}"
                     " --log={outfile}.log |"
                     " gzip > {outfile}")
    else:
        statement = ("ln -s {infile} {outfile}")

    P.run()


@follows(subSampleMothurGreenGenes)
@split(None, [os.path.join('mothur_greengenes_insilico_16S.dir', x, x + '_primers.txt')
              for x in AMPLICON_TARGETS.keys()])
def makeMothurGreenGenesDirs(infile, outfiles):
    '''Create directories for generating in silico amplicons'''
    P16S.makeInSilicoDirs(outfiles, AMPLICON_TARGETS)


###############################################################################
### Subtask: mothur's RDP database
###############################################################################
@follows(mkdir('mothur_rdp_insilico_16S.dir'))
@files([PARAMS['mothur_rdp_fasta'], PARAMS['mothur_rdp_taxonomy']],
       ['mothur_rdp_insilico_16S.dir/mothur_rdp_reference.fasta.gz',
        'mothur_rdp_insilico_16S.dir/mothur_rdp_reference_taxonomy.txt.gz'])
def parseMothurRDPDatabase(infiles, outfiles):
    '''Convert rdp database into a format suitable for in silico
    analysis. Outfiles must have suffix *reference* and prefix to match
    parent directory.
    '''
    fasta_in, tax_in = infiles
    fasta_out, tax_out = outfiles
    
    P16S.FormatRDP(fasta_in, tax_in, fasta_out, tax_out).run()


@transform(parseMothurRDPDatabase,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/subsampled.fasta.gz')
def subSampleMothurRDP(infiles, outfile):
    '''Downsample RDP database to specified number of sequences'''

    infile, tax_file = infiles
    infile = os.path.abspath(infile)

    # number of genes in greengenes database
    n = 0
    for i in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        n += 1
    L.info('There are %i sequences in %s' % (n, infile))

    if PARAMS['mothur_rdp_n_seq']:
        pptn = str(float(PARAMS['mothur_rdp_n_seq'])/n)

        statement = ("zcat {infile} |"
                     " python {scriptsdir}/fasta2fasta.py"
                     " --method=sample"
                     " --sample-proportion={pptn}"
                     " --random-seed={mothur_rdp_seed}"
                     " --log={outfile}.log |"
                     " gzip > {outfile}")
    else:
        statement = ("ln -s {infile} {outfile}")

    P.run()


@follows(subSampleMothurRDP)
@split(None, [os.path.join('mothur_rdp_insilico_16S.dir', x, x + '_primers.txt')
              for x in AMPLICON_TARGETS.keys()])
def makeMothurRDPDirs(infile, outfiles):
    '''Create directories for generating in silico amplicons'''
    P16S.makeInSilicoDirs(outfiles, AMPLICON_TARGETS)


###############################################################################
### Subtask: HOMD database
###############################################################################
@follows(mkdir('homd_insilico_16S.dir'))
@files([PARAMS['homd_fasta'], PARAMS['homd_taxonomy']],
       ['homd_insilico_16S.dir/homd_reference.fasta.gz',
        'homd_insilico_16S.dir/homd_reference_taxonomy.txt.gz'])
def parseHOMDDatabase(infiles, outfiles):
    '''Convert HOMD database into a format suitable for in silico
    analysis. Outfiles must have suffix *reference*.
    '''
    fasta_in, tax_in = infiles
    fasta_out, tax_out = outfiles
    
    # this is the same format at the E. coli db (i.e. nothing needs doing to it)
    P16S.FormatEColi(fasta_in, tax_in, fasta_out, tax_out).run()


@transform(parseHOMDDatabase,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/subsampled.fasta.gz')
def subSampleHOMD(infiles, outfile):
    '''No downsampling is done for the HOMD database'''

    infile, tax_file = infiles
    infile = os.path.abspath(infile)

    # number of genes in HOMD database
    n = 0
    for i in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        n += 1
    L.info('There are %i sequences in %s' % (n, infile))

    statement = ("ln -s {infile} {outfile}")
    P.run()


@follows(subSampleHOMD)
@split(None, [os.path.join('homd_insilico_16S.dir', x, x + '_primers.txt')
              for x in AMPLICON_TARGETS.keys()])
def makeHOMDDirs(infile, outfiles):
    '''Create directories for generating in silico amplicons'''
    P16S.makeInSilicoDirs(outfiles, AMPLICON_TARGETS)

###############################################################################
# 16S databases on which to perform in silco analysis
db_targets = {'greengenes': makeGreenGenesDirs,
              'mothur_greengenes': makeMothurGreenGenesDirs,
              'mothur_rdp': makeMothurRDPDirs,
              'homd': makeHOMDDirs
              }
DB_TARGETS=[]

for x in PARAMS['in_silico_databases'].split(','):
    DB_TARGETS.append(db_targets[x])


@follows(makeGreenGenesDirs, makeHOMDDirs)
def fetchDatabases():
    pass
    
###############################################################################
## Task: Extract amplicons
###############################################################################
@transform(DB_TARGETS,
           regex('(.+)/(.+)/(.+)_primers.txt'),
           add_inputs(r'\1/subsampled.fasta.gz'),
           r'\1/\2/\3_amplicons_s.fasta.gz')
def generateAmplicons(infiles, outfile):
    '''Use cutadapt to trim amplicon regions. At most cutadapt finds two primer
    matches. However, it outputs sequences for which either one, or two primer
    matches have been found. 
    '''
    
    primer_file, fasta_file = infiles

    # other outfiles
    out_untrimmed = P.snip(outfile, '.fasta.gz') + '_untrimmed.fasta.gz'
    out_details = P.snip(outfile, '.fasta.gz') + '_info.txt'

    # suffix for fasta headers
    suffix = os.path.basename(os.path.dirname(outfile))

    # get primers out of text file: forward\nreverse
    forward, reverse = [x.strip() for x in IOTools.openFile(primer_file).readlines()]
    reverse = P16S.reverseComplement(reverse)
    min_len = min(len(forward), len(reverse))
    error_rate = str(float(PARAMS['in_silico_errors'])/min_len)
    overlap = str(min_len - int(PARAMS['in_silico_overhang']))

    # run cutadapt
    statement = ("cutadapt"
                 " -g forward={forward}" # 5' adapter
                 " -a reverse={reverse}" # 3' adapter
                 " -n 2" # number of adapters to trim
                 " --error-rate={error_rate} "
                 " --overlap={overlap}"
                 " --match-read-wildcards"
                 " --info-file={out_details}"
                 " --untrimmed-output={out_untrimmed}"
                 " --suffix=_{suffix}"
                 " --output={outfile}"
                 " {fasta_file}")
    to_cluster = False
    P.run()


@transform(generateAmplicons,
           regex('(.+)/(.+)/(.+)_amplicons_s.fasta.gz'),
           add_inputs(r'\1/\2/\2_primers.txt'),
           r'\1/\2/\3_amplicons_as.fasta.gz')
def generateAntisenseAmplicons(infiles, outfile):
    '''For those intervals that were not trimmed, generate antisense
    sequence and re-search for primer matches, this is a sanity check'''
    fasta_file, primer_file = infiles

    in_untrimmed = P.snip(fasta_file, '_s.fasta.gz') + '_s_untrimmed.fasta.gz'

    out_untrimmed = P.snip(outfile, '.fasta.gz') + '_untrimmed.fasta.gz'
    out_details = P.snip(outfile, '.fasta.gz') + '_info.txt'

    # reverse complement the reads in the untrimmed fasta file 

    tmpf = os.path.join(os.path.dirname(outfile), 'temporary.fasta.gz')

    statement = ("zcat {in_untrimmed} |"
                 " python {scriptsdir}/fasta2fasta.py"
                 "  --method=reverse-complement"
                 "  --log={outfile}.log |"
                 " gzip > {tmpf}") 
    P.run()

    # suffix for fasta headers
    suffix = os.path.basename(os.path.dirname(outfile))

    # get primers out of text file: forward\nreverse
    forward, reverse = [x.strip() for x in IOTools.openFile(primer_file).readlines()]
    forward = P16S.reverseComplement(forward)
    min_len = min(len(forward), len(reverse))
    error_rate = str(float(PARAMS['in_silico_errors'])/min_len)
    overlap = str(min_len - int(PARAMS['in_silico_overhang']))

    # run cutadapt
    statement = ("cutadapt"
                 " -g forward={forward}" # 5' adapter
                 " -a reverse={reverse}" # 3' adapter
                 " -n 2" # number of adapters to trim
                 " --error-rate={error_rate} "
                 " --overlap={overlap}"
                 " --match-read-wildcards"
                 " --info-file={out_details}"
                 " --untrimmed-output={out_untrimmed}"
                 " --suffix=_{suffix}"
                 " --output={outfile}"
                 " {tmpf}")
    P.run()

    os.unlink(tmpf)


@collate([generateAmplicons, generateAntisenseAmplicons],
         regex('(.+)/(.+)(_s|_as).fasta.gz'),
         r'\1/\2_unfiltered.fasta.gz')
def collateAmplicons(infiles, outfile):
    '''Combine those amplicons on the forward strand with those
    found in the reverse complement
    '''
    to_cluster = False
    infiles = " ".join(infiles)
    statement = ("cat {infiles} > {outfile}")
    P.run()


@transform(collateAmplicons,
           suffix('_unfiltered.fasta.gz'),
           '.fasta.gz')
def filterAmplicons(amplicon_fasta, outfile):
    '''Parse cutadapt output to find sequences with both forward 
    and reverse primer match. Filter the fasta files to retain only
    these sequences.
    '''
    s_matches = P.snip(amplicon_fasta, '_unfiltered.fasta.gz') + '_s_info.txt'
    as_matches = P.snip(amplicon_fasta, '_unfiltered.fasta.gz') + '_as_info.txt'

    # find sense matches for both forward and reverse primer
    s_forward = set()
    s_reverse = set()

    for line in IOTools.openFile(s_matches):
        line = line.split()
        # sequences for which there was no primer match
        if int(line[1]) == -1:
            continue
        # cutadapt output can have varying number of fields...
        # ...search for primer 'name' instead (see above)
        elif 'forward' in line[2:]:
            s_forward.add(line[0])
        elif 'reverse' in line[2:]:
            s_reverse.add(line[0])
        else:
            raise ValueError('Unexpected line format: {}'.format(" ".join(line)))

    # get the sequences with both forward and reverse primer match
    s_sequence_matches = set.intersection(s_forward, s_reverse)

    # find antisense matches for both forward and reverse primer
    as_forward = set()
    as_reverse = set()

    for line in IOTools.openFile(as_matches):
        line = line.split()
        if int(line[1]) == -1:
            continue
        elif 'forward' in line[2:]:
            as_forward.add(line[0])
        elif 'reverse' in line[2:]:
            as_reverse.add(line[0])
        else:
            raise ValueError('Unexpected line format: {}'.format(" ".join(line)))

    as_sequence_matches = set.intersection(as_forward, as_reverse)


    # sanity check... there should be no sequences in both sense and antisense sets
    assert len(set.intersection(s_sequence_matches, as_sequence_matches)) == 0, \
        'There are sequences for which both sense and antisense primer matches occur'

    sequences_to_keep = tuple(s_sequence_matches.union(as_sequence_matches))

    # iterate over the cutadapt fasta and only output those sequences that have
    # both primers matching. 
    outf = IOTools.openFile(outfile, 'w')
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(amplicon_fasta)):
        header = re.sub('_V[1-9]*', '', fasta.title)
        if header in sequences_to_keep:
            outf.write('>' + fasta.title + "\n" + fasta.sequence + "\n")
        else:
            continue


@transform(filterAmplicons, suffix('.fasta.gz'), '_summary.txt.gz')
def summarizeAmplicons(infile, outfile):
    '''Get a summary of amplicon lengths'''

    statement = ("zcat {infile} |"
                 " python {scriptsdir}/fasta2table.py"
                 "  --section=length"
                 "  --split-fasta-identifier"
                 "  --log={outfile}.log |"
                 " gzip > {outfile}")
    P.run()


@collate(summarizeAmplicons,
         regex('(.+)/(.+)/(.+)_amplicons_summary.txt.gz'),
         r'\1/amplicon_length_summary.txt.gz')
def collateAmpliconLengthSummary(infiles, outfile):
    '''Collate summaries of amplicon lengths for different primers'''

    to_cluster = False
    infiles = " ".join(infiles)
    statement = ("sleep 2m;"
                 " python {scriptsdir}/combine_tables.py"
                 "  --cat=Region"
                 "  --regex-filename='.+/(.+)_amplicons_summary.txt.gz'"
                 "  --log={outfile}.log"
                 " {infiles} |"
                 " gzip > {outfile}")
    P.run()


@jobs_limit(1, 'RGlobal')
@transform(collateAmpliconLengthSummary, suffix('.txt.gz'), '.png')
def plotAmpliconLengthSummary(infile, outfile):
    '''Plot histogram of distribution of amplicon lengths'''

    R('''rm(list=ls())
         require('ggplot2')
         df = read.table('{infile}', header=TRUE, sep='\t',
                         stringsAsFactors=FALSE)
         hist <- ggplot(df, aes(x=length))
         hist <- hist  + geom_histogram(fill=c('darkblue'))
         hist <- hist  + facet_wrap(~Region, scales='free')
         hist <- hist + theme_bw() + theme(panel.grid=element_blank())
         png('{outfile}')
         plot(hist)
         dev.off()
         rm(list=ls())'''.format(**locals()))

    
@jobs_limit(1)
@follows(collateAmpliconLengthSummary)
@transform(summarizeAmplicons, suffix('_summary.txt.gz'), '_outliers.txt')
def findAmpliconLengthOutliers(infile, outfile):
    '''Create a list of IDs for those amplicons that have a length
    greater than 2 standard devations from the mean'''

    df = pd.read_table(infile, compression='gzip')

    # get those entries that are too long
    df_long = df[df['length'] > np.mean(df['length']) + 2*np.std(df['length'])]
    # and those that are too short
    df_short = df[df['length'] < np.mean(df['length']) - 2*np.std(df['length'])]

    df_out = pd.concat([df_long, df_short])

    L.info('There are %i outliers in file %s' % (len(df_out.index), outfile))
    df_out.to_csv(outfile, sep='\t', index=False)


@follows(findAmpliconLengthOutliers)
@transform(filterAmplicons,
           regex('(.+)/(.+).fasta.gz'),
           add_inputs(r'\1/\2_outliers.txt'),
           r'\1/\2_filtered01.fasta.gz')
def filterAmpliconsOnLength(infiles, outfile):
    '''Filter amplicon fasta files to remove those of spurious length '''
    fasta_file, length_file = infiles

    # a list of the sequence IDs to be removed from fasta file
    to_remove = []
    header = True
    for a in IOTools.openFile(length_file):
        if header:
            header = False
            continue
        to_remove.append(a.split()[0])

    # open outfiles
    out_filtered = IOTools.openFile(outfile, 'w')
    out_removed = P.snip(outfile, '_filtered01.fasta.gz') + '_removed01.fasta.gz'
    out_removed = IOTools.openFile(out_removed, 'w')

    # write out fasta sequences
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        if fasta.title in to_remove:
            out_removed.write( "\n".join(['>' + fasta.title, fasta.sequence]) + '\n')
        else:
            out_filtered.write( "\n".join(['>' + fasta.title, fasta.sequence]) + '\n')
    
    out_filtered.close()
    out_removed.close()


@transform(filterAmpliconsOnLength,
           suffix('_filtered01.fasta.gz'),
           '_filtered02.fasta.gz')
def filterAmpliconsForNs(infile, outfile):
    '''Remove those amplicons that contain Ns in sequence'''
    outf_discarded = P.snip(outfile, '_filtered02.fasta.gz') + '_removed02.fasta.gz'

    outf_f = IOTools.openFile(outfile, 'w')
    outf_d = IOTools.openFile(outf_discarded, 'w')

    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        if re.search('N+', fasta.sequence):
            outf_d.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
        else:
            outf_f.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
    
    outf_f.close()
    outf_d.close()


@collate(filterAmpliconsForNs,
         regex('(.+)/(.+)/(.+)_filtered02.fasta.gz'),
         r'\1/amplicon_sequence_intersections.txt')
def findAmpliconIntersection(infiles, outfile):
    '''Find the those sequence ids that are present in all the
    subsampled amplicon fasta files. Also output a file containing
    the count of the number of sequences in each amplicon set. 
    '''
    count_table = P.snip(outfile, '_intersections.txt') + '_counts.txt'
    count_table = IOTools.openFile(count_table, 'w')
    count_table.write('Region\tSequenceNumber\n')

    first = True
    for infile in infiles:
        region = os.path.basename(os.path.dirname(infile))
        fasta_headers = [x.title for x in 
                         FastaIterator.FastaIterator(IOTools.openFile(infile))]
        fasta_headers = [re.sub('_V[1-9]*', '', x) for x in fasta_headers]
        count_table.write(region + '\t' + str(len(fasta_headers)) + '\n')
        if first:
            first = False
            headers_to_keep = set(fasta_headers)
        else:
            headers_to_keep &= set(fasta_headers)

    count_table.close()

    with IOTools.openFile(outfile, 'w') as outf:
        for header in headers_to_keep:
            outf.write(header + '\n')


@merge(findAmpliconIntersection, 'insilico_db_sequence_counts.tsv')
def countAmpliconIntersections(infiles, outfile):
    '''For each database, count the number of sequences that are present
    in each of the subsampled amplicon fasta files'''

    outf = IOTools.openFile(outfile, 'w')
    outf.write('Database\tNSeq\n')
    for infile in infiles:
        n_seq = len(IOTools.openFile(infile).readlines())
        dataset = os.path.basename(os.path.dirname(infile))
        dataset = P.snip(dataset, '_insilico_16S.dir')
        outf.write(dataset + '\t' + str(n_seq) + '\n')
    outf.close()


@follows(findAmpliconIntersection)
@transform(filterAmpliconsForNs,
           regex('(.+)/(.+)/(.+)_filtered02.fasta.gz'),
           add_inputs(r'\1/amplicon_sequence_intersections.txt'),
           r'\1/\2/\3_filtered03.fasta.gz')
def filterAmpliconIntersection(infiles, outfile):
    '''Filter fasta files containing amplicons of subregions so they only
    contain sequences that are in the intersection between all subregion
    fasta files.
    '''
    fasta_file, to_keep = infiles
    
    to_keep = [x.strip() for x in IOTools.openFile(to_keep).readlines()]

    outf = IOTools.openFile(outfile, 'w')
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        title = re.sub('_V[1-9]*', '', fasta.title)
        if title in to_keep:
            outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
        else:
            continue
    outf.close()


@collate(filterAmpliconIntersection,
         regex('(.+)/.+/.+_filtered03.fasta.gz'),
         r'\1/amplicons.sentinel')
def checkAmpliconFastas(infiles, outfile):
    '''A sanity check to ensure the same sequences are present in each
    subregion fasta'''
    
    infiles = list(infiles)
    first_inf = infiles.pop(0)
    ref_headers = [re.sub('_V[1-9]*', '', x.title) for x in 
                   FastaIterator.FastaIterator(IOTools.openFile(first_inf))]
    ref_set = set(ref_headers)

    # check for duplicates in first fasta
    assert len(ref_headers) == len(ref_set), 'Duplicate headers in fasta file {}'.format(first_inf)

    for infile in infiles:
        headers = [re.sub('_V[1-9]*', '', x.title) for x in 
                   FastaIterator.FastaIterator(IOTools.openFile(infile))]
        header_set = set(headers)
        # check for duplicates in fasta
        assert len(headers) == len(header_set), 'Duplicate IDs in fasta: {}'.format(infile)
        # check that fasta contains same headers as reference
        assert set(headers) == ref_set, ('Different sequence IDs in fasta: {}'
                                         ' compared to {}'.format(infile, first_inf))
    
    open(outfile, 'w').close()    


###################
# Parse the filtered fasta files for their headers. filterAmpliconIntersection
# Find the full taxonomy for each header parse
# regex(homd_insilico_16S.dir/homd_reference_taxonomy.txt.gz)
# For each full taxonomy, select one header at random
# Output a list of the selected fasta headers with taxonomy. 

@transform(filterAmpliconIntersection,
           regex('(.+)_insilico_16S.dir/V1_V9/V1_V9_amplicons_filtered03.fasta.gz'),
           add_inputs(r'\1_insilico_16S.dir/\1_reference_taxonomy.txt.gz'),
           r'\1_insilico_16S.dir/amplicons_to_keep.tsv')
def selectSingleAmpliconForEachSpecies(infiles, outfile):
    '''Parse the filtered amplicons for the V1-V9 data. Fetch their
    taxonomy. Select only a single amplicon for each unique species
    '''

    fasta_file, tax_file = infiles

    assert os.path.exists(tax_file)
    
    # parse the fasta file
    amplicon_ids = []
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        ID = P.snip(fasta.title, '_V1_V9')
        amplicon_ids.append(ID)
    
    # parse the taxonomy file 
    taxonomy_dict = collections.defaultdict(list)
    for line in IOTools.openFile(tax_file):
        ID, taxonomy = line.split()
        if ID in amplicon_ids:
            taxonomy_dict[taxonomy].append(ID)
        else:
            continue

    # iterate over the taxonomy dictionary, select a single ID for each taxon
    # output ID, taxonomy, number of IDs
    with IOTools.openFile(outfile, 'w') as outf:
        for taxon, IDs in taxonomy_dict.iteritems():
            n = len(IDs) 
            if n == 1:
                outf.write('\t'.join([IDs[0], str(n), taxon]) + '\n')
            else:
                random.shuffle(IDs)
                outf.write('\t'.join([IDs.pop(), str(n), taxon]) + '\n')

@follows(selectSingleAmpliconForEachSpecies)
@transform(filterAmpliconIntersection,
           regex('(.+)/(.+)/(.+)_filtered03.fasta.gz'),
           add_inputs(r'\1/amplicons_to_keep.tsv'),
           r'\1/\2/\3_filtered04.fasta.gz')
def filterUniqueSpecies(infiles, outfile):
    '''Select only a single amplicon sequence for each species'''

    fasta_file, to_keep = infiles
    to_keep = [x.split()[0] for x in IOTools.openFile(to_keep)]

    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
            # There are sometimes underscores in the fasta headers
            title = re.sub('_V1_V2|_V1_V3|_V1_V9|_V3_V5|_V4|_V6_V9', '', fasta.title)
            if title in to_keep:
                outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
            else:
                continue

           
            
@follows(plotAmpliconLengthSummary,
         countAmpliconIntersections,
         checkAmpliconFastas)
def createInSilicoAmplicons():
    pass


###############################################################################
### Subtask: Summarize In silico primers lost due to filtering
###############################################################################
@follows(filterAmpliconsOnLength, filterAmpliconsForNs)
@transform(collateAmplicons,
           regex('(.+)_unfiltered.fasta.gz'),
           add_inputs(r'\1_removed01.fasta.gz', r'\1_removed02.fasta.gz'),
           r'\1_discarded_summary.tsv')
def countPrimerMatchFailures(infiles, outfile):
    '''Parse cutadapt output to find sequences that failed due to a lack of
    matching primers. 
    '''
    amplicon_fasta, length_fasta, n_fasta = infiles
    n_fail_both = 0
    n_fail_forward = 0
    n_fail_reverse = 0

    s_matches = P.snip(amplicon_fasta, '_unfiltered.fasta.gz') + '_s_info.txt'
    as_matches = P.snip(amplicon_fasta, '_unfiltered.fasta.gz') + '_as_info.txt'

    # find sense matches for both forward and reverse primer
    s_forward = set()
    s_reverse = set()

    for line in IOTools.openFile(s_matches):
        line = line.split()
        # sequences for which there was no primer match
        if int(line[1]) == -1:
            n_fail_both += 1
            continue
        # cutadapt output can have varying number of fields...
        # ...search for primer 'name' instead (see above)
        elif 'forward' in line[2:]:
            s_forward.add(line[0])
        elif 'reverse' in line[2:]:
            s_reverse.add(line[0])
        else:
            raise ValueError('Unexpected line format: {}'.format(" ".join(line)))

    # find the number of sequences with forward but no reverse and vice versa
    n_fail_forward += len(s_reverse - s_forward)
    n_fail_reverse += len(s_forward - s_reverse)
        
    # find antisense matches for both forward and reverse primer
    as_forward = set()
    as_reverse = set()

    for line in IOTools.openFile(as_matches):
        line = line.split()
        if int(line[1]) == -1:
            continue
        elif 'forward' in line[2:]:
            as_forward.add(line[0])
        elif 'reverse' in line[2:]:
            as_reverse.add(line[0])
        else:
            raise ValueError('Unexpected line format: {}'.format(" ".join(line)))

    # find the number of sequences with forward but no reverse and vice versa
    n_fail_forward += len(as_reverse - as_forward)
    n_fail_reverse += len(as_forward - as_reverse)

    # find the number of sequences that failed on length
    n_fail_length = 0
    for line in FastaIterator.FastaIterator(IOTools.openFile(length_fasta)):
        n_fail_length += 1

    # find the number of sequences that failed because they contained N's
    n_fail_n = 0
    for line in FastaIterator.FastaIterator(IOTools.openFile(n_fasta)):
        n_fail_n +=1

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('Reason\tSequences\n')
        outf.write('No primer match\t' + str(n_fail_both) + '\n')
        outf.write('No forward primer match\t' + str(n_fail_forward) + '\n')
        outf.write('No reverse primer match\t' + str(n_fail_reverse) + '\n')
        outf.write('Amplicon too short\t' + str(n_fail_length) + '\n')
        outf.write('Amplicon contains Ns\t' + str(n_fail_n) + '\n')


@collate(countPrimerMatchFailures,
         regex('(.+)/.+/.+_discarded_summary.tsv'),
         r'\1/lost_sequence_summary.tsv')
def collatePrimerMatchFailures(infiles, outfile):
    '''Combine tables summarizing lost reads'''
    infiles = " ".join(infiles)
    statement = ("sleep 2m;"
                 " python {scriptsdir}/combine_tables.py"
                 " --add-file-prefix"
                 " --regex-filename='(.+)_amplicons_discarded_summary.tsv'"
                 " --log={outfile}.log"
                 " {infiles}"
                 " > {outfile}")
    cluster_options = '-l walltime=00:03:00'
    P.run()

    
###############################################################################
## Task: Run RDP Classifier
###############################################################################
@subdivide(filterUniqueSpecies,
           regex('(.+)/(.+)/(.+)_filtered04.fasta.gz'),
           r'\1/\2/*/\3_filtered04.fasta.gz')
def createRDPThresholdSubdirectories(infile, outfiles):
    '''Create a subdirectory for each of the RDP thresholds to be run'''
    outfile = os.path.basename(infile)
    out_parent_dir = os.path.dirname(infile)

    for threshold in PARAMS['mothur_cutoff'].split(','):
        directory = os.path.join(out_parent_dir, 'threshold_' + str(threshold))
        if not os.path.exists(directory):
            os.makedirs(directory)
            outf = os.path.join(directory, outfile)
            to_cluster = False
            statement = "cp -f {infile} {outf}"
            P.run()


@transform(createRDPThresholdSubdirectories,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)/(.+).fasta.gz'),
           add_inputs([r'\1_insilico_16S.dir/\1_reference.fasta.gz',
                       r'\1_insilico_16S.dir/\1_reference_taxonomy.txt.gz']),
           r'\1_insilico_16S.dir/\2/\3/\4.\1_reference_taxonomy.wang.taxonomy')
def runMothurRDP(infiles, outfile):
    '''Run mothur's version of the RDP classifier from the commandline'''
    input_fasta, reference_files = infiles
    ref_fasta, ref_taxonomy = reference_files

    # get cut off from directory name
    cut_off = os.path.basename(os.path.dirname(input_fasta))
    cut_off = cut_off.split('_')[1]

    # mothur can't handle gzipped files and has no option for redirecting
    # output... runnning mothur from each subdirectory
    working_dir = os.path.dirname(input_fasta)
    tmp_input_fasta = P.snip(input_fasta, '.gz')
    tmp_ref_fasta = os.path.join(working_dir, 
                                 P.snip(os.path.basename(ref_fasta), '.gz'))
    tmp_ref_tax = os.path.join(working_dir,
                               P.snip(os.path.basename(ref_taxonomy), '.gz'))

    statement = ("gunzip -c {input_fasta} > {tmp_input_fasta};"
                 " gunzip -c {ref_fasta} > {tmp_ref_fasta};"
                 " gunzip -c {ref_taxonomy} > {tmp_ref_tax}")
    to_cluster = False
    P.run()

    tmp_input_fasta = os.path.abspath(tmp_input_fasta)
    tmp_ref_fasta = os.path.abspath(tmp_ref_fasta)
    tmp_ref_tax = os.path.abspath(tmp_ref_tax)

    # set outfile
    outfile = P.snip(input_fasta, '.fasta.gz') + \
        P.snip(os.path.basename(ref_taxonomy), '.txt.gz') + \
        '.wang.taxonomy'

    # run mothur classify.seqs()
    statement = ("cd {working_dir};"
                 " mothur \"#classify.seqs("
                 "fasta={tmp_input_fasta},"
                 " template={tmp_ref_fasta},"
                 " taxonomy={tmp_ref_tax},"
                 " cutoff={cut_off},"
                 " processors={mothur_processors},"
                 " iters={mothur_iters})\";"
                 " cd -")

    cluster_options = '-l walltime=192:00:00,nodes=1:ppn={}'.format(PARAMS['mothur_processors']) 
    to_cluster = True

    P.run()

    os.unlink(tmp_input_fasta)
    os.unlink(tmp_ref_fasta)
    os.unlink(tmp_ref_tax)


###############################################################################
### Subtask: Summarize RDP Classification
###############################################################################
@transform(runMothurRDP,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)/(.+).wang.taxonomy'),
#           regex('(.+)_insilico_16S.dir/(V1_V2)/(threshold_30)/(.+).wang.taxonomy'),
           add_inputs(r'\1_insilico_16S.dir/\1_reference_taxonomy.txt.gz'),
           r'\1_insilico_16S.dir/\2/\3/tax_assignment_summary.txt.gz')
def summarizeMothurRDPClassification(infiles, outfile):
    '''Calculate the proportion of sequences that are correctly assigned
    at each taxonomic level'''
    classified_file, reference_file = infiles

    # set up the taonomy file for easy inputting into 
    tmpf = P.getTempFilename('.')
    statement = "zcat {reference_file} | sed 's/\\t/;/' | sed 's/;$//' > {tmpf}"
    to_cluster = False
    P.run()

    # read the taxonomy file into a dataframe...
    df = pd.read_table(tmpf, sep=';', header=None, index_col=0)
    df.columns = ['Kingdom', 'Phylum', 'Class',
                  'Order', 'Family', 'Genus', 'Species']
    df.index.name='SeqID'    
    
    # write the dataframe of original taxonomies to a flatfile for reference
    out_ref = P.snip(outfile, '_summary.txt.gz') + '_reference_table.txt.gz'
    df.to_csv(IOTools.openFile(out_ref, 'w'), sep='\t')

    os.unlink(tmpf)
    
    # write out summaries of taxonomic assignment success
    outf = IOTools.openFile(outfile, 'w')
    out_ambig = P.snip(outfile, '_summary.txt.gz') + '_ambig.txt.gz'
    out_ambig = IOTools.openFile(out_ambig, 'w')
    out_full = P.snip(outfile, '.txt.gz') + '_full.txt.gz'
    out_full = IOTools.openFile(out_full, 'w')

    # Subset the reference taxonomy so that it only contains the relvent seqs
    to_keep = [x.split()[0] for x in IOTools.openFile(classified_file)]
    to_keep = [re.sub('_V[1-9]*', '', x) for x in to_keep]
    df = df.loc[to_keep,:]
    L.warn('Succesfully subset reference dataframe to %i entries' % len(df.index))

    completed = 0
    for line in IOTools.openFile(classified_file):
        line = line.strip(';')
        seqID, taxon = line.split('\t')
        # get the orignal sequence identification number
        seqID = re.sub('_V[1-9]*', '', seqID)

        # discard the match threshold (either int or float)
        taxon = re.sub('\(\d+\.?\d*\)', '', taxon)
        tax_cl = taxon.split(';')[:-1] # there is a ';' at the end of the taxon line.    
        # check that there are indeed 7 levels in the taxonomic classification
        # output
        assert len(tax_cl) == 7, \
            'Ambiguous taxonomic classification output {}'.format(tax_cl)
    

        # fetch the lowest classification level that is not 'unclassified'
        tax_cl_level = None
        if tax_cl[0] == 'unclassified':
            tax_cl_level = -1 * len(tax_cl)
            tax_id = 'k__unclassified'
        else:
            for level, tax in enumerate(tax_cl[::-1], start=1):
                if tax != 'unclassified':
                    tax_cl_level = -1 * level
                    # GREENGENES: [] are removed
                    tax_id = re.sub('^\[|\]$', '', tax)
                    break

        # fetch the reference classification for this sequence
        # WARNING: index values are sometimes int... see above
        ref_cl = df.loc[seqID]

        # fetch the taxonomic level to which sequence was classified
        ref_cl_level = ref_cl.index[tax_cl_level]

        # find out whether it was classified correctly at this level
        # GREENGENES: [] are removed
        ref_id = re.sub('^\[|\]$', '', ref_cl[tax_cl_level])
        if ref_id == tax_id:
            out_cl = ref_cl_level
            # sanity check... consistency of parent taxon, if possible
            if tax_cl_level > len(ref_cl) * -1:
                try:
                    a = re.sub('^\[|\]$', '', ref_cl[tax_cl_level-1])
                    b = re.sub('^\[|\]$', '', tax_cl[tax_cl_level-1])
                except IndexError:
                    raise IndexError('Index Out of Bounds:\n %s \n %s \n %s' % \
                                     (str(tax_cl_level), str(ref_cl), tax_cl))

                if a != b:
                    warn = "Wrong tax assignment:\n {} {}".format(str(tax_cl), str(ref_cl))
                    L.warn(warn)
                    out_cl = ref_cl.index[tax_cl_level - 1] + '_error'
                    warn = re.sub('\n', '\t', warn)
                    out_ambig.write(warn + '\n')

        else:
            out_cl = ref_cl_level + '_error'

        outf.write(seqID + '\t' + out_cl + '\n')
        out_full.write('\t'.join(map(str, list(ref_cl))) + '\t' + out_cl + '\n')

        completed += 1
        if completed % 1000 == 0: L.warn('Completed %i reads for file %s' % \
                                         (completed, classified_file))

        
    outf.close()
    out_full.close()
    out_ambig.close()


@transform(summarizeMothurRDPClassification,
           suffix('.txt.gz'),
           '_table.txt')
def countMothurRDPAssignmentSuccess(infile, outfile):
    '''Summarize the number of sequences assigned correctly'''
    
    infile = P.snip(infile, '.txt.gz') + '_full.txt.gz'

    statement = ("printf 'Count\\tClassification\\n' > {outfile};"
                 " zcat {infile} | awk '{{print $NF}}' | sort | uniq -c |"
                 " awk '{{print $1\"\\t\"$2}}' >> {outfile}")
    to_cluster = False
    P.run()


@collate(countMothurRDPAssignmentSuccess,
         regex('(.+)/(.+)/(.+)/tax_assignment_summary_table.txt'),
         r'\1/rdp_\3_tax_assignment_summary_table.txt')
def collateMothurRDPAssignmentSuccess(infiles, outfile):

    headers = [x.split('/')[-3] for x in infiles]
    headers = ','.join(headers)
    
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 "  --columns=2"
                 "  --skip-titles"
                 "  --header-names={headers}"
                 "  --log={outfile}.log"
                 " {infiles}"
                 " > {outfile}")
    to_cluster = True
    P.run()


@jobs_limit(1, 'RGlobal')
@transform(collateMothurRDPAssignmentSuccess,
           suffix('_table.txt'),
           '.pdf')
def plotMothurRDPAssignmentSuccess(infile, outfile):
    '''Plot stacked bar chart of the proportion of reads correctly
    classified'''

    # get the rdp threshold
    threshold = P.snip(os.path.basename(infile),
                       '_tax_assignment_summary_table.txt')
    infile = os.path.abspath(infile)

    R('''
      rm(list=ls())
      plotBar <- function(df){
      df[is.na(df)] <- 0


      # Collapse all errors above family
      df_error <- df[rownames(df) %in% c('Order_error', 'Class_error', 'Phylum_error', 'Kingdom_error'),]
      other_error <- colSums(df_error)
      df <- rbind(df, other_error)
      rownames(df)[length(rownames(df))] <- 'Other_error'
      df <- df[!rownames(df) %in% c('Order_error', 'Class_error', 'Phylum_error', 'Kingdom_error'),]

      # Collapse all taxa above family
      df_other <- df[rownames(df) %in% c('Order', 'Class', 'Phylum', 'Kingdom'),]
      other <- colSums(df_other)
      df <- rbind(df, other)
      rownames(df)[length(rownames(df))] <- 'Other'
      df <- df[!rownames(df) %in% c('Order', 'Class', 'Phylum', 'Kingdom'),]

      # reorder dataframe
      o <- c('Species_error', 'Genus_error', 'Family_error', 'Other_error',
             'Species', 'Genus', 'Family', 'Other')
      df <- df[o,]

      # set fields for melting as a factor with the correct order
      df$order <- factor(rownames(df),
                         levels=c('Species_error', 'Genus_error', 'Family_error', 'Other_error',
                                  'Species', 'Genus', 'Family', 'Other'))

      # melt the dataframe 
      df.m <- melt(df, idvar=c('order'))

      my_cols=c(rev(brewer.pal(4, 'Reds')), "#2EB789", "#1CBCC2", "#26B1E6", "#001814")

      # plot 
      pl <- ggplot(df.m, aes(x=variable, y=value, fill=order)) + geom_bar(stat='identity')
      pl <- pl + theme_bw() + theme(panel.grid = element_blank())
      pl <- pl + scale_fill_manual(values=my_cols)
      pl <- pl + theme(legend.position="none")
      plot(pl)
      }''')

    # and plot...

    R('''
      require('ggplot2')
      require('reshape2')
      require('RColorBrewer')

      df <- read.table('{infile}',
                       sep='\t',
                       header=TRUE,
                       row.names=1,
                       stringsAsFactors=FALSE,
                       na.strings='na')
      pdf('{outfile}', height=2.5, width=3)
      plotBar(df)
      dev.off()'''.format(**locals()))   


@collate(collateMothurRDPAssignmentSuccess,
         regex('(.+)_insilico_16S.dir/(.+)_tax_assignment_summary_table.txt'),
         r'\1_insilico_tax_assignment_summary_table.load')
def loadMothurRDPAssignmentSuccess(infiles, outfile):
    '''
    '''
    table_name = P.snip(os.path.basename(outfile), '.load')

    to_cluster = False
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 "  --cat=RDP_Threshold"
                 "  --regex-filename='.+/rdp_threshold_(.+)_tax_assignment_summary_table.txt'"
                 "  --log={outfile}.log"
                 " {infiles} |"
                 " sed 's/na/0/g' |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 " > {outfile}")
    P.run()


@jobs_limit(1, 'RGlobal')
@follows(findAmpliconIntersection)
@transform(loadMothurRDPAssignmentSuccess,
           regex('(.+)_insilico_tax_assignment_summary_table.load'),
           add_inputs(r'\1_insilico_16S.dir/amplicons_to_keep.tsv'),
           r'\1_insilico_16S.dir/rdp_threshold_comparison.pdf')
def plotMothurRDPThresholdAssignmentSuccess(infiles, outfile):
    '''For each amplicon region, for each RDP classification threshold,
    plot the proportion of sequences successfully classified to species
    level.
    '''
    count_tab, total = infiles
    total = len(IOTools.openFile(total).readlines())

    database = P.snip(count_tab, '_insilico_tax_assignment_summary_table.load')
    table = P.snip(count_tab, '.load')
    statement = "SELECT * FROM {} WHERE bin == 'Species'".format(table)
    df = pd.read_sql(sql=statement, con=connect())
    df = df.drop('bin', axis=1)
    df = pd.melt(df, id_vars='RDP_Threshold')

    dfr = pandas2ri.py2ri(df)
    R('''rm(list=ls())''')
    R.assign('dfr', dfr)

    R('''
      require('ggplot2')
      pl <- ggplot(dfr, aes(x=RDP_Threshold,
                            y=value/{total}*100,
                            colour=variable))
      pl <- pl + geom_line()
      pl <- pl + theme_bw()
      pl <- pl + xlab('RDP Classification Threshold')
      pl <- pl + ylab('% Reads Correctly Assigned')
      pl <- pl + ggtitle('{database}\n({total} Reads)')

      pdf('{outfile}')
      plot(pl)
      dev.off()'''.format(**locals()))

    
###############################################################################
### Subtask: Create MetaCoder Fasta Files
###############################################################################
@follows(createRDPThresholdSubdirectories)
@transform(summarizeMothurRDPClassification,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)/tax_assignment_summary.txt.gz'),
           add_inputs(r'\1_insilico_16S.dir/\2/\3/\2_amplicons_filtered04.fasta.gz'),
           r'\1_insilico_16S.dir/\2/\3/\2_amplicons_annotated.fasta.gz')
def createMetaCoderFasta(infiles, outfile):
    '''Parse the input fasta, output any sequence that was not correctly
    identified with the correct taxonomy information in the fasta header
    '''

    in_tax, in_fasta = infiles
    in_tax = P.snip(in_tax, '.txt.gz') + '_full.txt.gz'

    # parse the taxonomy file, add anything that is not correct to species
    # do dictionary
    tax_dict = {}
    for line in IOTools.openFile(in_tax):
        line = line.strip()
        if line.endswith('Species'):
            continue

        # line is K,P,C,O,F,G,SeqID,Classification\n
        # generate taxonomy in format applicable for metacoder headers
        line = line.split()
        seq_id  = line.pop(-2)
        taxonomy = line[:-1]
        
        assert seq_id not in tax_dict.keys()
        tax_dict[seq_id] = taxonomy

        
    L.warn('Finished creating tax dict for sample %s' % in_fasta)

    with IOTools.openFile(outfile, 'w') as outf:
        # parse the fasta file, for those sequence headers in the tax_dict,
        # output the sequence with taxonomy added to header
        n = 0
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(in_fasta)):
            n += 1
            if n % 1000 == 0:
                L.warn('Processed %i sequences' % n)
            header = re.sub('_V\d.*', '', fasta.title)
            
            tax = tax_dict.get(header, False)
            if tax:
                outf.write('>' + fasta.title + '\t' + tax + '\n' \
                           + fasta.sequence + '\n')


@transform(filterUniqueSpecies,
           regex('(.+)_insilico_16S.dir/V1_V9/V1_V9_amplicons.filtered04.fasta.gz'),
           add_inputs(r'\1_insilico_16S.dir/\1_reference_taxonomy.txt.gz'),
           r'\1_insilico_16S.dir/\1_metacoder.fasta.gz')
def createMasterMetaCoderFasta(infiles, outfile):
    '''Create a metacoder fasta file for all sequences involved in the
    in silico analysis'''

    in_fasta, in_tax = infiles
    
    # parse the taxonomy file, add anything that is not correct to species
    # do dictionary
    i = 0
    tax_dict = {}
    for line in IOTools.openFile(in_tax):
        seq_id, taxonomy = line.split()
        # tax line ends with ;sequenceID;
        taxonomy = ';'.join(taxonomy.split(';')[:-2])

        tax_dict[seq_id] = taxonomy
        i += 1
        if i % 10000 == 0:
            L.warn('Added %i out of 93463 sequences to mapping_dict' % i)
        
    assert len(tax_dict.keys()) == i        
    L.warn('Finished creating tax dict for sample %s' % in_fasta)

    with IOTools.openFile(outfile, 'w') as outf:
        # parse the fasta file, for those sequence headers in the tax_dict,
        # output the sequence with taxonomy added to header
        n = 0
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(in_fasta)):
            n += 1
            if n % 1000 == 0:
                L.warn('Processed %i sequences' % n)
            header = re.sub('_V\d.*', '', fasta.title)
            
            tax = tax_dict.get(header, False)
            if tax:
                outf.write('>' + header + '\t' + tax + '\n' \
                           + fasta.sequence + '\n')
    
    


    
###############################################################################
## Task: Create In Silico OTUs
###############################################################################
@jobs_limit(1)
@transform(filterUniqueSpecies,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)_filtered04.fasta.gz'),
           r'\1_insilico_otus.dir')
def makeInSilicoOTUDirs(infile, outfile):
    '''Create directories for in silico OTU analysis'''
    if not os.path.exists(outfile):
        os.mkdir(outfile)    


@follows(makeInSilicoOTUDirs)
@subdivide(filterUniqueSpecies,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)_filtered04.fasta.gz'),
           r'\1_insilico_otus.dir/\2/dist*')
def makeInSilicoOTUSubdirs(infile, outfiles):
    '''Make a subdirectory for each difference threshold used for OTU
    identification.
    '''
    thresholds = map(str, PARAMS['usearch_thresholds'].split(','))
    thresholds = [x.zfill(2) for x in thresholds]

    outdir = re.sub('_16S', '_otus', os.path.dirname(infile)) 
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for threshold in thresholds:
        outfile = 'dist' + threshold
        outfile = os.path.join(outdir, outfile)
        if not os.path.exists(outfile):
            os.mkdir(outfile)


@follows(makeInSilicoOTUSubdirs)
@transform(filterUniqueSpecies,
           regex('(.+)_insilico_16S.dir/(.+)/(.+)_filtered04.fasta.gz'),
           r'\1_insilico_otus.dir/\2/\3_dereplicated.fasta.gz')
def dereplicateInSilicoAmplicons(infile, outfile):
    '''Derepilcate sequences that appear in all amplicon subsets'''

    outf_tmp = P.snip(outfile, '.gz')

    tmpf = P.getTempFilename('.')
    statement = ("gunzip -c {infile} > {tmpf};"
                 " usearch -derep_fulllength {tmpf} -fastaout {outf_tmp};"
                 " rm {tmpf};"
                 " gzip {outf_tmp}")
    P.run()


@transform(dereplicateInSilicoAmplicons,
           regex('(.+)/V1_V9/(.+)_dereplicated.fasta.gz'),
           r'\1/\2_sequence_order.txt')
def outputAmpliconOrder(infile, outfile):
    '''Output the order of sequences in the V1-V9 dereplicated amplicon
    file
    '''
    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            header = re.sub('_V[1-9]*', '', fasta.title)
            outf.write(header + '\n')
 
   
@follows(outputAmpliconOrder)
@transform(dereplicateInSilicoAmplicons,
           regex('(.+)/(.+)/(.+).fasta.gz'),
           add_inputs(r'\1/V1_V9_amplicons_sequence_order.txt'),
           r'\1/\2/\3_sorted.fasta.gz')
def sortAmplicons(infiles, outfile):
    '''Sort in silico amplicons for all subregions so they're in the
    same order as the V1-V9 amplicons
    '''
    in_fasta, out_order = infiles

    # there must be a better way to do this...
    seq_dict = {}
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(in_fasta)):
        header = re.sub('_V[1-9]*', '', fasta.title)
        assert header not in seq_dict.keys(), 'Duplicate headers in fasta file {}'.format(header)
        seq_dict[header] = fasta.sequence

    with IOTools.openFile(outfile, 'w') as outf:
        for seq_id in IOTools.openFile(out_order):
            seq_id = seq_id.strip()
            out_id = '>' + seq_id + ';size=1'
            try:
                out_seq = seq_dict[seq_id]
                outf.write(out_id + '\n' + out_seq + '\n')
            except KeyError:
                # This gives a lot of warnings!
                L.warn('%s missing from %s' % (seq_id, in_fasta))
             
@follows(sortAmplicons)
@transform(makeInSilicoOTUSubdirs,
           regex('(.+)/(.+)/dist(.+)'),
           add_inputs(r'\1/\2/\2_amplicons_dereplicated_sorted.fasta.gz'),
           r'\1/\2/dist\3/\2_otus.fa')
def createInSilicoOTUs(infiles, outfile):
    '''Cluster in silico amplicons to form OTUs'''

    indir, in_fasta = infiles
    tmpf = P.getTempFilename('.')
    out_up = P.snip(outfile, '.fa') + '_up.txt'

    # get radius %age from the directory name
    rd_pc = indir[-1] + '.0'

    statement = ("gunzip -c {in_fasta} > {tmpf};"
                 "usearch"
                 " -cluster_otus {tmpf}"
                 " -otus {outfile}"
                 " -uparseout {out_up}"
                 " -otu_radius_pct {rd_pc}"
                 " -uparse_break -999" # R.E. suggestion for removal of chimera check
                 " -relabel OTU_")
    cluster_options = '-l walltime=04:00:00,nodes=1:ppn=8'
    P.run()

    os.unlink(tmpf)


@transform(createInSilicoOTUs, suffix('.fa'), '.count')
def countInSilicoOTUs(infile, outfile):
    
    statement = ("echo OTU_Count > {outfile};"
                 "cat {infile} | grep '>' | wc -l >> {outfile}")
    P.run()


@collate(countInSilicoOTUs,
         regex('(.+)/(.+)/(.+)/.+_otus.count'),
         r'\1/\3_summary.tsv')
def collateInSilicoOTUs(infiles, outfile):
    '''Collate counts of OTUs across amplicons'''
    
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --cat=Amplicon"
                 " --regex-filename='.+/(.+)_otus.count'"
                 " --log={outfile}.log"
                 " {infiles} > {outfile}")
    P.run()


@collate(collateInSilicoOTUs,
         regex('(.+)_insilico_otus.dir/.+_summary.tsv'),
         r'\1_insilico_otus.dir/\1_otu_summary.tsv')
def collateInSilicoOTUs2(infiles, outfile):
    '''Collate counts of OTUs across distance'''
    
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --cat=Distance"
                 " --regex-filename='.+/(.+)_summary.tsv'"
                 " --log={outfile}.log"
                 " {infiles} > {outfile}")
    P.run()


@transform(collateInSilicoOTUs2,
           regex('.+/(.+).tsv'),
           r'\1.load')
def loadInSilicoOTUs(infile, outfile):
    '''load the summaries of per-base entropy'''
    table_name = P.snip(os.path.basename(outfile), '.load')
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 " > {outfile}")
    P.run()
    

@jobs_limit(1, 'RGlobal')
@transform(collateInSilicoOTUs2, suffix('.tsv'), '_linePlot.pdf')
def plotInSilicoOTUs(infile, outfile):
    '''Plot summaries of number of OTUs called per amplicon region'''

    R('''rm(list=ls())''')
    outf_bar = P.snip(outfile, '_linePlot.pdf') + '_barPlot.pdf'

    df = pd.read_table(infile, header=0)
    df = pd.melt(df, id_vars=['Distance', 'Amplicon',],
                 value_name='Number_of_OTUs')

    R.assign('df', df)
    
    R('''
      require('ggplot2')
      x <- ggplot(df, aes(x=factor(Distance),
                      y=Number_of_OTUs,
                      colour=Amplicon,
                      fill=Amplicon))
      x <- x  + geom_bar(position='dodge', 
                         stat='identity')
      x <- x + theme_bw()

      pdf('{outf_bar}')
      plot(x)
      dev.off()

      y <- ggplot(df, aes(x=factor(Amplicon),
                          y=Number_of_OTUs,
                          colour=Distance,
                          group=Distance))
      y <- y  + geom_line() + geom_point() + theme_bw()

      pdf('{outfile}')
      plot(y)
      dev.off()'''.format(**locals()))


###############################################################################
### Subtask: Classify OTUs and subsequently classify amplicons
###############################################################################
# @follows(runMothurRDP)
# @transform(createInSilicoOTUs,
#            regex('(.+)_insilico_otus.dir/(.+)/(.+)/(.+)_otus.fa'),
#            add_inputs(r'\1_insilico_16S.dir/\2/threshold_80/\4_amplicons_filtered03.\1_reference_taxonomy.wang.taxonomy'),
#            r'\1_insilico_otus.dir/\2/\3/\4_otus_taxonomic_id.tsv')
# def classifyOTUs(infiles, outfile):
#     '''Find the RDP classification for OTU at 80% confidence threshold
#     '''
 
#     otu_fasta, tax_assignment = infiles
#     otu_id = P.snip(otu_fasta, '.fa') + '_up.txt'


#     # parse taxonomic assignment into a table that can be loaded as dataframe
#     tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')
#     tmpf = P.getTempFilename('.')
#     statement = ("cat {tax_assignment} |"
#                  " python {tax_parser}"
#                  "  --remove-confidence-thresholds"
#                  "  --out-format=table"
#                  "  --log={outfile}.log"
#                  " > {tmpf}")
#     P.run()

#     df_tax = pd.read_table(tmpf, header=0, index_col=0)
#     # remove the suffix from the SequenceIDs
#     seq_ids = df_tax.index.tolist()
#     seq_ids = [re.sub('_V[1-9]*', '', x) for x in seq_ids]
#     df_tax.index = seq_ids

#     # read in detail of the OTUs and their representative sequence IDs
#     seqIDs = []
#     otuIDs = []
#     for line in IOTools.openFile(otu_id):
#         otuIDs.append(line[-1])
#         seqIDs.append(line[0].split(';')[0])

#     df_otu = {'OTU_ID' : pd.Series(otuIDs, index=seqIDs),}
#     df_otu = pd.DataFrame(df_otu)

#     df = df_otu.merge(df_tax, how='inner', left_index=True, right_index=True)

#     df.to_csv(outfile, sep='\t')

#     os.unlink(tmpf)


# @follows(filterAmpliconIntersection)
# @transform(createInSilicoOTUs,
#            regex('(.+)_insilico_otus.dir/(.+)/(.+)/(.+)_otus.fa'),
#            add_inputs(r'\1_insilico_16S.dir/\2/\2_amplicons_filtered03.fasta.gz'),
#            r'\1_insilico_otus.dir/\2/\3/\4_amplicon_otu_assignment.txt')
# def assignAmpliconsToOTUs(infiles, outfile):
#     '''Get the OTU for each in silico amplicon'''
#     otu_fasta, amplicon_fasta = infiles
    
#     tmp_amplicon_fasta = P.snip(amplicon_fasta, '.gz')

#     # determine the classification threshold used to generate OTUs (1-3)
#     id_threshold = os.path.basename(os.path.dirname(otu_fasta))
#     id_threshold = str((100 - int(id_threshold[-1]))/100.0)

#     statement = ("zcat {amplicon_fasta} > {tmp_amplicon_fasta};"
#                  " usearch"
#                  "  -usearch_global {tmp_amplicon_fasta}"
#                  " -db {otu_fasta}"
#                 " -strand plus"
#                  " -id {id_threshold}")
#     P.run()  


###############################################################################
## Task: Calculate Entropy Across Aligned V1-V9
###############################################################################
# Something is screwed with the dependencies in this pipeline. The symlink in
# subsampleGreengenes step is not maintaining timestamp. Attempts to update
# this with touch are failing due to differences between the clock (??) when
# using touch and using ruffus.
# To fix this, I'm temporarily restarting the pipeline from the this step
# having removed dependency...

@jobs_limit(1)

#@transform(filterAmpliconIntersection,
@transform('*_insilico_16S.dir/V1_V9/*_filtered03.fasta.gz',
           regex('(.+)_insilico_16S.dir/(V1_V9)/(.+)_filtered03.fasta.gz'),
           add_inputs(PARAMS['data_mock_ref_db']),
           r'\1_insilico_entropy.dir/aligned_sequence_baseline.fa')
def extractBaselineSequence(infiles, outfile):
    '''Pull out baseline e coli sequence along which to calculate
    entropy.
    '''
    o = os.path.dirname(outfile)
    if not os.path.exists(o):
        os.mkdir(o)

    out_dir, infile = infiles
    base_seq_id = PARAMS['entropy_base_sample']
    check = False
    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            if fasta.title.startswith(base_seq_id):
                outf.write('>BASE\n' + fasta.sequence + '\n')
                check = True
            else:
                continue
    assert check, "Couldn't find base sequence in ref db"


@follows(extractBaselineSequence)
#@transform(filterAmpliconIntersection,
@transform('*_insilico_16S.dir/V1_V9/*_filtered03.fasta.gz',
           regex('(.+)_insilico_16S.dir/V1_V9/(.+)_filtered03.fasta.gz'),
           add_inputs(r'\1_insilico_entropy.dir/aligned_sequence_baseline.fa'),
           r'\1_insilico_entropy.dir/aligned_sequence.fasta.gz')
def alignV1_V9Amplicons(infiles, outfile):
    '''Align the V1-V9 amplicons, plus the baseline sequence
    using muscle
    '''

    seq_fasta, base_fasta = infiles
    tmpf = P.getTempFilename('.')
    tmp_outf = P.snip(outfile, '.gz')

    to_cluster = False
    statement = ("cat {base_fasta} > {tmpf};"
                 " zcat {seq_fasta} >> {tmpf};"
                 " muscle "
                 "  -in {tmpf}"
                 "  -out {tmp_outf}"
#                 "  -maxhours 48"
                 " 2> {outfile}.log &&"
                 " gzip {tmp_outf}")
    P.run()


@transform(alignV1_V9Amplicons, suffix('.fasta.gz'), '_entropy.tsv')
def  calculateEntropy(infile, outfile):
    '''Calculate the shannon entropy for multiple sequence alignments'''

    entropy_scr = os.path.join(os.path.dirname(__file__),
                               'alignment2entropy.py')
    statement = ("python {entropy_scr}"
                 " -a {infile}"
                 " -b BASE"
                 " -o {outfile}")
    P.run()
           


@transform(alignV1_V9Amplicons, suffix('.fasta.gz'), '_smoothed_entropy.tsv')
def  calculateSmoothedEntropy(infile, outfile):
    '''Calculate the shannon entropy for multiple sequence alignments'''

    entropy_scr = os.path.join(os.path.dirname(__file__),
                               'alignment2entropy.py')
    statement = ("python {entropy_scr}"
                 " -a {infile}"
                 " -b BASE"
                 " -s 70"
                 " -o {outfile}")
    to_cluster=False
    P.run()


    
@jobs_limit(1)
@transform([calculateSmoothedEntropy, calculateEntropy],
           suffix('.tsv'),
           add_inputs(PARAMS['in_silico_amplicon_bed']),
           '_table.tsv')
def findAmpliconEntropy(infiles, outfile):
    '''Calculate summaries of the entropy across different in silico
    amplicons
    '''

    entropy_tab, amplicon_bed = infiles

    # read in the entropy table and extract entropies as list
    df = pd.read_table(entropy_tab, index_col=0, header=0)
    entropy = df['Entropy'].tolist()

    interval_dict = {}
    # read through bed intervals and extract 
    for line in IOTools.openFile(amplicon_bed):
        line = line.split()
        start, end, interval  = line[1:4] 
        sub_entropy = entropy[int(start):int(end)]
        tot = sum(sub_entropy)
        mn = np.mean(sub_entropy)
        md = np.median(sub_entropy)

        interval_dict[interval] = [tot, mn, md]

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('Amplicon\tTotalEntropy\tMeanEntropy\tMedianEntropy\n')
        for key, value in interval_dict.iteritems():
            outf.write(key + '\t' + '\t'.join(map(str, value)) + '\n')


@transform(findAmpliconEntropy,
           regex('(.+)_insilico_entropy.dir/aligned_sequence(.*)_entropy_table.tsv'),
           r'\1\2_insilico_entropy.load')
def loadAmpliconEntropy(infile, outfile):
    ''' '''
    table_name = P.snip(os.path.basename(outfile), '.load')
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 " > {outfile}")
    P.run()


# @jobs_limit(1, 'RGlobal')
# #@follows(loadInSilicoOTUs)
# @transform(loadAmpliconEntropy,
#            regex('(.+)_insilico_entropy.load'),
#            r'\1_insilico_entropy.dir/\1_total_entropy_vs_otu_dist01.pdf')
# def plotAmpliconEntropy_vs_OTUs(infile, outfile):
#     '''Plot the number of OTUs that are generated as a function of the
#     median/total entropy for a subregion'''
#     entropy_table = P.snip(infile, '.load', strip_path=True)
#     amplicon_table = P.snip(entropy_table, '_insilico_entropy') + '_otu_summary'

#     for dist in ('dist01', 'dist02', 'dist03'):
#         R('''rm(list=ls())''')
#         statement = ("SELECT a.*, b.OTU_Count FROM {entropy_table} AS a"
#                      " INNER JOIN {amplicon_table} AS b"
#                      " ON a.Amplicon == b.Amplicon"
#                      " WHERE Distance = '{dist}'".format(**locals()))
#         df = pd.read_sql(sql=statement, con=connect())

#         dfr = pandas2ri.py2ri(df)
#         R.assign('df', dfr)

#         for entropy in ('MeanEntropy', 'MedianEntropy', 'TotalEntropy'):
#             outf = P.snip(outfile, '_dist01.pdf') + '_' + dist + '.pdf'
#             if entropy == 'MeanEntropy':
#                 outf = re.sub('_total_', '_mean_', outf)
#                 print outf
#             elif entropy == 'MedianEntropy':
#                 outf = re.sub('_total_', '_median_', outf)
#                 print outf
#             else:
#                 print outf
#                 pass
                
#             R('''require(ggplot2)''')
#             R('''pl = ggplot(df, aes(x= %s, y=OTU_Count, label=Amplicon))''' % entropy)
#             R('''pl = pl + geom_point() + geom_text() + theme_bw()''')
#             R('''pdf('%s')''' % outf)
#             R('''plot(pl)''')
#             R('''dev.off()''')


@follows(loadAmpliconEntropy)
def full():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
