
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
# Metatask: ANALYSIS OF MOCK COMMUNITY
###############################################################################
###############################################################################
## Task: Run crossmatch to compare 16S reads to reference db
###############################################################################
###############################################################################
### Subtask: Obtain reads from PacBio ccs run output.
###############################################################################
@follows(mkdir('fasta_cleaned.dir'))
@merge(PARAMS['data_ccs_bams'].split(';'), 'fasta_cleaned.dir/mock16s.bam')
def mergeInputBams(infiles, outfile):
    '''Takes bamfiles output by new smartseq ccs pipeline sorts them 
    and merges them'''

    # sort the bamfiles
    tempfiles = []
    for f in infiles:
        tmpf = P.getTempFilename('.')
        statement = ("samtools sort"
                     " -T {tmpf}"
                     " -o {tmpf}.bam"
                     " -O bam"
                     " {f}")
        to_cluster = False
        P.run()
        tempfiles.append(tmpf)

    # merge sorted bamfiles, note that header is taken from in1.bam
    tempfiles = [x + '.bam' for x in tempfiles]
    tmpfiles = " ".join(tempfiles)
    statement = ("samtools merge -f {outfile} {tmpfiles};"
                 " samtools index {outfile}")
    to_cluster = False
    P.run()

    # remove the temporary files
    for f in tempfiles:
        os.unlink(f)
        os.unlink(P.snip(f, '.bam'))


@transform(mergeInputBams, suffix('.bam'), '.fasta.gz')
def convertInputBamToFasta(infile, outfile):
    '''Take input bamfile and convert to fasta'''
    statement = ("samtools bam2fq -n {infile} |"
                 " awk 'NR % 4==1 {{print $0}} NR % 4==2 {{print $0}}' |"
                 " sed 's/@/>/' |"
                 " gzip > {outfile}")
    P.run()


###############################################################################
### Subtask: Remove adapter sequences using cutadapt
###############################################################################
@transform(convertInputBamToFasta,
           regex('.+/(.+).fasta.gz'),
           r'fasta_cleaned.dir/\1_trimmed.fasta.gz')
def runCutAdapt(input_fasta, outfile):
    '''Trim primers from PacBio reads using cutadapt'''
    
    f_adapter = PARAMS['cutadapt_forward_adapter']
    r_adapter = PARAMS['cutadapt_reverse_adapter']
    r_adapter = P16S.reverseComplement(r_adapter)

    # create outfiles for different stages of adapter trimming...
    # info files to contain details of primer match location
    info_forward = P.snip(outfile, '.fasta.gz') + '_info_fwd.txt'
    info_forward_rc = P.snip(info_forward, '.txt') + '_rc.txt'
    info_reverse = P.snip(info_forward, '_fwd.txt') + '_rev.txt'
    # outfile to contain reads with trimmed forward primer
    out_forward = P.snip(outfile, '.fasta.gz') + '_forward.fasta.gz'
    # outfile to contain those reads that failed to find primers
    failed_fwd = P.snip(outfile, '_trimmed.fasta.gz') + '_failed_fwd.fasta'
    failed_rev = P.snip(outfile, '_trimmed.fasta.gz') + '_failed_rev.fasta'

    # temporary files for the untrimmed output in first step
    tmpf_rejected = 'ctmp_cutadapt_untrimmed.fasta'
    tmpf_rc = 'ctmp_cutadapt_reverse_complement.fasta'

    # various run parameters
    f_length = len(f_adapter)
    r_length = len(r_adapter)
    length = min([r_length, f_length])
    error_rate = float(PARAMS['in_silico_errors'])/length
    params = ("--error-rate={}"
              " --match-read-wildcards".format(error_rate))
    

    # Run cutadapt to trim the forward primer. 
    to_cluster = False
    statement = ("cutadapt"
                 " -g forward={f_adapter}"
                 " --info-file={info_forward}"
                 " --untrimmed-output={tmpf_rejected}"
                 " {params}"
                 " {input_fasta} 2> {outfile}.1.log |"
                 " gzip > {out_forward}".format(**locals()))
    P.run()

    # reverse complement sequences for which adapter wasn't found. 
    with IOTools.openFile(tmpf_rc, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(tmpf_rejected)):
            outf.write('>' + fasta.title + '\n')
            outf.write(P16S.reverseComplement(fasta.sequence) + '\n')

    # worried about read/write speed...
    assert os.path.exists(tmpf_rc), 'File is missing %s' % tmpf_rc

    # rerun cutadapt to trim the forward primer from the reverse complement
    to_cluster = False
    statement = ("cutadapt"
                 "  -g forward={f_adapter}"
                 "  --info-file={info_forward_rc}"
                 "  --untrimmed-output={failed_fwd}"
                 "  {params}"
                 "  {tmpf_rc} 2> {outfile}.2.log |"
                 " gzip >> {out_forward};"
#                 " checkpoint;" ??
                 " cutadapt"
                 "  -a reverse={r_adapter}"
                 "  --info-file={info_reverse}"
                 "  --untrimmed-output={failed_rev}"
                 "  {params}"
                 "  {out_forward} 2> {outfile}.3.log |"
                 " gzip > {outfile}".format(**locals()))
    P.run()
 

@transform(runCutAdapt, suffix('.fasta.gz'), '_filtered.fasta.gz')
def filterTrimmedOutputByLength(infile, outfile):
    '''Filter cutadapt cleaned fasta file to remove reads over and under
    a specified length
    '''

    min_len = PARAMS['filter_min_length']
    max_len = PARAMS['filter_max_length']

    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            if len(fasta.sequence) >= min_len and len(fasta.sequence) <= max_len:
                outf.write( '\n'.join([ '>' + fasta.title,
                                       fasta.sequence]) + '\n')
            else:
                continue


@transform(filterTrimmedOutputByLength,
           suffix('_trimmed_filtered.fasta.gz'),
           '_filtering_summary.tsv')
def summarizeMockSequenceFiltering(infile, outfile):
    '''Count the number of reads at each filtering step'''
    filtered_fasta = infile
    trimmed_fasta = P.snip(infile, '_filtered.fasta.gz') + '.fasta.gz'
    failed_fwd_fasta = P.snip(trimmed_fasta, '_trimmed.fasta.gz') + '_failed_fwd.fasta'
    failed_rev_fasta = P.snip(trimmed_fasta, '_trimmed.fasta.gz') + '_failed_rev.fasta'

    def _get_n(inf):
        n = 0
        for line in FastaIterator.FastaIterator(IOTools.openFile(inf)):
            n += 1
        return n

    filtered = _get_n(infile)
    trimmed = _get_n(trimmed_fasta)
    failed_trimming = trimmed - filtered
    failed_fwd = _get_n(failed_fwd_fasta)
    failed_rev = _get_n(failed_rev_fasta)

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('Passed Filtering' + '\t' + str(filtered) + '\n')
        outf.write('Failed Length Filtering' + '\t' + str(failed_trimming) + '\n')
        outf.write('Failed fwd Primer Match' + '\t' + str(failed_fwd) + '\n')
        outf.write('Failed rev Primer Match' + '\t' + str(failed_rev) + '\n')

###############################################################################
### Subtask: Map reads to refdb using cross-match.
###############################################################################
@follows(mkdir('fasta_cross_match.dir'))
@subdivide(filterTrimmedOutputByLength,
           regex('.+/(.+)_trimmed_filtered.fasta.gz'),
           r'fasta_cross_match.dir/\1_*.fasta.gz')
def splitFasta(infile, outfiles):
    '''Split the fasta file into n smaller files to speed up cross_match'''
    
    outfile_stub = P.snip(os.path.basename(infile),
                          '_trimmed_filtered.fasta.gz')
    outfile_stub = os.path.join('fasta_cross_match.dir', outfile_stub)
    multiSplitFasta(infile, outfile_stub, PARAMS['cross_match_subdivide'])



       
CROSS_MATCH_TARGETS = []
for run in  PARAMS['cross_match_parameter_options'].split(';'):
    # construct the outfile directory
    out_dir = 'cross_match' +  re.sub('\W+', '_', run)
    out_dir = re.sub('_discrep_lists_tags', '', out_dir)
    out_dir = os.path.join('fasta_cross_match.dir', out_dir, 'cross_match.sentinel')
    CROSS_MATCH_TARGETS.append(str(out_dir))


@follows(splitFasta)
@split(None, CROSS_MATCH_TARGETS)
def makeCrossMatchDirs(infile, outfile):
    '''Create multiple directories for cross_match runs.
    Create sentinel as dependency checker, it also
    contains the run parameters.'''

    for run in  PARAMS['cross_match_parameter_options'].split(';'):
        # construct the outfile directory
        out_dir = 'cross_match' +  re.sub('\W+', '_', run)
        out_dir = re.sub('_discrep_lists_tags', '', out_dir)
        out_dir = os.path.join('fasta_cross_match.dir', out_dir)
        outfile = os.path.join(out_dir, 'cross_match.sentinel')
        
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        with IOTools.openFile(outfile, 'w') as outf:
            outf.write(re.sub(',', ' ', run))


@combinatorics.product(splitFasta, formatter("fasta_cross_match.dir/(?P<STUB>.+).fasta.gz"),
                       makeCrossMatchDirs, formatter(".+/(?P<SUBDIR>.+)/(.+)"),
                       add_inputs(PARAMS['data_mock_ref_db']),
                       "fasta_cross_match.dir/{SUBDIR[1][0]}/{STUB[0][0]}_cross_match.txt")
def runCrossMatch(infiles, outfile):
    '''Run cross_match using a variety of different parameters, 
    which are specified in the config file
    '''

    # because of the way the input tsks are passed...
    fasta_out_dir_sentinel, ref_db = infiles
    fasta_file, out_dir_sentinel = fasta_out_dir_sentinel
    out_dir = os.path.dirname(out_dir_sentinel)

    # cross_match run parameters are in the sentinel file
    params = open(out_dir_sentinel).readline()

    # cross_match can't handle gzipped file format
    tmp_inf = P.getTempFilename(out_dir)

    statement = ("zcat {fasta_file} > {tmp_inf};"
                 " cross_match {params} {tmp_inf} {ref_db} > {outfile}")
    P.run()
    
    # cross_match doesn't give an option to rename logfiles
    shutil.move(tmp_inf + '.log', outfile + '.log')
    os.unlink(tmp_inf)


@collate(runCrossMatch,
         regex('(.+)/(.+)/(.+)_[0-9]+_cross_match.txt'),
         r'\1/\2/\3_\2.txt')
def collateCrossMatchResults(infiles, outfile):
    '''Concatenate all the crossmatch output into a single file
    '''
    infiles = ' '.join(infiles)
    to_cluster = False
    statement = ('cat {infiles} > {outfile}')
    P.run()


@transform(collateCrossMatchResults,
           suffix('.txt'),
           '_readsMapped.txt')
def countReadMappingSuccess(infile, outfile):
    '''Count the number of reads mapped'''
    mapped = set()
    for line in IOTools.openFile(infile):
        if line.startswith('ALIGNMENT'):
            mapped.add(line.split()[5])
    file_name = P.snip(infile, '.txt', strip_path=True)
    outf = IOTools.openFile(outfile, 'w')
    outf.write(file_name + '\t' + str(len(mapped)) + '\n')

    
@collate(countReadMappingSuccess,
         regex('(.+)/.+/.+.txt'),
         r'\1/cross_match_read_mapping_success.tsv')
def collateReadMappingSuccess(infiles, outfile):
    '''Summarize read mapping success''' 
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --cat=MappingParameters"
                 " --log={outfile}.log"
                 " {infiles} > {outfile}")
    to_cluster = False
    P.run()
    
    
@transform(collateCrossMatchResults,
           suffix('.txt'),
           add_inputs(PARAMS['data_mock_ref_db']),
           '_summary.txt')
def summarizeCrossMatch(infiles, outfile):
    '''parse crossmatch output and write summary tables.
    For output format see  www.phrap.org/phredphrap/phrap.html
    '''

    infile, ref_db = infiles

    # WARNING!
    if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
        L.warn('There is no sanity check for duplicate matches in'
               ' default cross_match run') 

    # get the maximum length of the sequences in the reference database
    # (informs output table size)
    max_len = 0
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(ref_db)):
        if len(fasta.sequence) > max_len:
            max_len = len(fasta.sequence)

    # parse cross match output
    out_dict = P16S.parseCrossMatchDiscrepList(infile, max_len, outfile)

    # pickle and dump out_dict
    out_pickle = P.snip(outfile, '.txt') + '.p'
    pickle.dump(out_dict, open(out_pickle, 'wb'))


@transform(summarizeCrossMatch,
           suffix('_summary.txt'),
           '_readMappingSuccess.txt')
def summarizeCrossMatchMappingSuccess(infile, outfile):
    '''read through the cross match output and count the number of reads
    successfully mapped to each reference sequence'''
    
    # a list of the target sequence IDs
    targets  = []
    header = True
    for line in IOTools.openFile(infile):
        if header:
            header = False
            assert line.split()[9] == 'TargetSeqID', "Problem with infile format"
            continue
        targets.append(line.split()[9])

    outf = IOTools.openFile(outfile, 'w')
    outf.write('TargetSequenceID\tReadsMapped\n')
    for seqID in set(targets):
        cnt = targets.count(seqID)
        outf.write(seqID + '\t' + str(cnt) + '\n')
    outf.close()

    
@follows(mkdir('mapped_read_lengths.dir'))
@subdivide(summarizeCrossMatch,
          regex('(.+)/.+_summary.txt'),
           r'\1/mapped_read_lengths.dir/*_mapped_reads.tsv')
def fetchAlignedReads(infile, outfile):

    reads_dict = collections.defaultdict(list)
    header = True
    for line in IOTools.openFile(infile):
        if header:
            header = False
            continue
        line = line.split()
        reads_dict[line[9]].append(line[4])

    out_dir = os.path.join(os.path.dirname(infile), 'mapped_read_lengths.dir')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for k, v in reads_dict.iteritems():
        outfile = os.path.join(out_dir, k + '_mapped_reads.tsv')
        with IOTools.openFile(outfile, 'w') as outf:
            outf.write('\n'.join(v))
        
    
    # length_dict = collections.defaultdict(list)
    # header = True
    # for line in IOTools.openFile(infile):
    #     if header:
    #         header = False
    #         continue
    #     line = line.split()
    #     # $6-QueryStart, $7-QueryEnd, $8-QueryTail,$10-TargetSeqId,
    #     q_start = int(line[5])
    #     q_end = int(line[6])
    #     q_tail = int(line[7])
    #     target_seq = line[9]
    #     target_seq_length = q_start + q_end + q_tail -1
    #     length_dict[target_seq].append(target_seq_length)

    # pickle.dump(length_dict, open(outfile, 'wb'))

@transform(fetchAlignedReads,
           suffix('s.tsv'),
           add_inputs(filterTrimmedOutputByLength),
           '_lengths.tsv')
def fetchAlignedReadLengths(infiles, outfile):
    '''Parse the cross_match aligned sequences and fetch the length
    distributions of reads aligned to each reference transform'''
    read_file, fasta_file = infiles

    # read the aligned reads into a list
    aligned_reads = []
    for line in IOTools.openFile(read_file):
        aligned_reads.append(line.strip())

    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
            if fasta.title in aligned_reads:
                outf.write(str(len(fasta.sequence)) + '\n')
            else:
                continue


    
@transform(summarizeCrossMatchMappingSuccess,
           regex('.+/(.+).txt'),
           r'\1.load')
def loadCrossMatchMappingSuccess(infile, outfile):
    ''' '''
    table_name = P.snip(infile, '.txt', strip_path=True)
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 " --table={table_name}"
                 " > {outfile}")
    to_cluster = False
    P.run()


@transform(summarizeCrossMatchMappingSuccess,
           suffix('.txt'),
           add_inputs(loadRefDBTaxSummary),
           '.png')
def plotCrossMatchMappingSuccess(infiles, outfile):
    '''plot number of reads successfully mapped to each reference 
    sequence
    '''
    map_table, ref_table = infiles
    ref_table = P.snip(os.path.basename(ref_table), '.load')
    ref_table = re.sub('-', '_', ref_table)

    df_map = pd.read_table(map_table, index_col=0)

    statement = 'SELECT StrainID,StrainName FROM {}'.format(ref_table)
    df_ref = pd.read_sql_query(statement,
                               con=connect(),
                               index_col='StrainID')
    sp = df_ref['StrainName'].tolist()
    sp = ['_'.join(x.split()[:2]) for x in sp]
    df_ref['SpName'] = sp

    df = df_map.join(df_ref)

    R('''rm(list=ls())''')
    R.assign('df', df)
    R('''png('{}')
         par(mar=c(15.1,4.1,4.1,2.1))
         barplot(df$ReadsMapped,
                 names.arg=df$SpName,
                 las=2)
         dev.off()'''.format(outfile))


@jobs_limit(1, 'RGlobal')
@transform(summarizeCrossMatch, suffix('.txt'), '_alignment_length_hist.png')
def plotCrossMatchAlignmentSummary(infile, outfile):
    '''Plot histogram of the alignment length and alignment
    score distributions.'''
    R('''rm(list=ls())''')

    df = pd.read_table(infile, sep='\t', index_col='QuerySeqID')
    
    # get the alignment lengths
    lengths = df['QueryEnd'] - df['QueryStart'] + 1
    lengths = rpy2.robjects.vectors.IntVector(lengths.tolist())
    R.assign('lengths', lengths )
    # plot 
    R('''
      png('{}')
      hist(lengths,
      breaks=100,
      main='Histogram of aligned length distributions')
      dev.off()
      '''.format(outfile))

    # get the SWA scores
    R('''rm(list=ls())''')
    outfile = P.snip(outfile, '_alignment_length_hist.png') + '_scores_hist.png'
    scores = rpy2.robjects.vectors.IntVector(df['Score'].tolist())
    R.assign('scores', scores)
    R('''
      png('{}')
      hist(scores,
      breaks=100,
      main='Hisgtogram of SW alignment scores')
      dev.off()
      '''.format(outfile))


@jobs_limit(1)
@transform(summarizeCrossMatch,
           regex('.+/(.+).txt'),
           r'\1.load')
def loadCrossMatchSummary(infile, outfile):
    '''Load the cross_match ALIGNMENT matches into database.
    Beware trying to access sqlite database in parallel
    '''

    table_name = P.snip(infile, '.txt', strip_path=True)
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 " --table={table_name}"
                 " --add-index=QuerySeqID"
                 " > {outfile}")

    P.run()


@jobs_limit(1)
@follows(plotCrossMatchAlignmentSummary)
@subdivide(summarizeCrossMatch,
           regex('(.+)/.+_summary.txt'),
           add_inputs(loadRefDBTaxSummary),
           r'\1/*_errorProfile.pdf')
def plotErrorHistograms(infiles, outfiles):
    '''Plot error histograms
    '''
    pickle_file, ref_table = infiles
    pickle_file = P.snip(pickle_file, '.txt') + '.p'
    ref_table = P.snip(ref_table, '.load', strip_path=True)
    ref_table = re.sub('-', '_', ref_table)
    outdir = os.path.dirname(pickle_file)

    # get a table containing strain name and strain id. 
    statement = ("SELECT StrainID,StrainName"
                 " FROM {ref_table}".format(**locals()))
    ref_df = pd.read_sql_query(sql=statement, con=connect(),
                               index_col='StrainID')

    # read in the table containing all mismatches (this is large!)
    df = pd.DataFrame(pickle.load(open(pickle_file, 'rb')))

    # iterate over every strain ID in the mismatch table
    for strain in set(df.columns.get_level_values(0)):
        # fetch strain name from database
        strain_name = str(ref_df.loc[strain, 'StrainName'])
        # subset only that strain from table
        df_sub = df[strain]

        # the number of reads mapping to this strain
        n_seq = str(len(df_sub.columns))

        # pandas counts empty strings as True
        df_sub = df_sub.applymap(lambda x: np.nan if x == '' else x)

        outfile_name = strain + '_errorProfile.pdf'
        outfile_name = os.path.join(outdir, outfile_name)
        
        R('''rm(list=ls())''')
        R('''pdf('%s', height=3.5)''' % outfile_name)
        R('''par(mfrow=c(1,3))''')
        for error in ('S', 'I', 'D'):
            if error == 'D':
                df_sub1 = df_sub.applymap(lambda x: np.nan
                                          if x in ['S', 'I'] else x)
            elif error == 'S':
                df_sub1 = df_sub.applymap(lambda x: np.nan
                                          if x in ['D', 'I'] else x)
            elif error == 'I':
                df_sub1 = df_sub.applymap(lambda x: np.nan
                                          if x in ['S', 'D'] else x)
            else:
                df_sub1 = df_sub.copy()
            df_sub1 = df_sub1.count(axis=1)
            df_sub1 = df_sub1.to_frame()
            df_sub1.reset_index(inplace=True)
            df_sub1.columns = ['Index', 'Count']

            # set to robjects
            x_axis = rpy2.robjects.vectors.IntVector(df_sub1['Index'].tolist())
            y_axis = rpy2.robjects.vectors.IntVector(df_sub1['Count'].tolist())

            R.assign('xa', x_axis)
            R.assign('ya', y_axis)

            R('''
              plot(xa,
                   ya/{n_seq}*100,
                   type='l',
                   main=paste('{strain_name}', '\n(', {n_seq}, ' sequences aligned)', sep=''),
                   xlab='Location',
                   ylab='% Error',
                   ylim=c(0,100))
              '''.format(**locals()))

        R('''dev.off()''')
        R('''rm(list=ls())''')
        L.info('Plots completed for strain {}, written to {}'.format(strain, outfile_name))


@jobs_limit(1)
@subdivide(summarizeCrossMatch,
           regex('(.+)/(.+)_summary.txt'),
           add_inputs(PARAMS['data_mock_ref_db'], loadRefDBTaxSummary),
           r'\1/\2_*_DeletionsVSHomopolymers.pdf')
def plotDeletionsVsHomopolymerRuns(infiles, outfiles):
    '''Plot the homopolymer runs in conjunction with the deletion
    frequency for alignments to each reference sequence'''
    
    pickle_file, ref_fasta, ref_table = infiles
    infile = pickle_file

    # get a table containing strain name and strain id. 
    ref_table = P.snip(ref_table, '.load', strip_path=True)
    ref_table = re.sub('-', '_', ref_table)
    statement = ("SELECT StrainID,StrainName"
                 " FROM {ref_table}".format(**locals()))
    ref_df = pd.read_sql_query(sql=statement, con=connect(),
                               index_col='StrainID')


    # create a dictionary out of reference sequences
    fasta_dict = {}
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(ref_fasta)):
        strain_id = fasta.title.split()[0]
        fasta_dict[strain_id] = fasta.sequence

    pickle_file = P.snip(pickle_file, '.txt') + '.p'
    df = pd.DataFrame(pickle.load(open(pickle_file, 'rb')))

    # iterate over every strain ID in the mismatch table
    for strain in set(df.columns.get_level_values(0)):
        outfile = P.snip(infile, '_summary.txt')
        outfile = outfile + '_' + strain + '_DeletionsVSHomopolymers.pdf'

        strain_name = str(ref_df.loc[strain, 'StrainName'])

        R('''rm(list=ls())''')
        # subset only that strain from table
        df_sub = df[strain]

        # the number of reads mapping to this strain
        n_seq = str(len(df_sub.columns))

        # Write a tmep matrix of the Deletions at each position in each
        # aligned read
        df_sub = df_sub.applymap(lambda x: 1 if x == 'D' else 0)
        tmpf_d = P.getTempFilename('.')
        df_sub.to_csv(tmpf_d, header=False)

        # write a tempfile of homopolymers for this strain
        strain_seq = fasta_dict[strain]
        homopolymers = {}
        for index, run in P16S.homopolymer_iterator(strain_seq,
                                                    index_pos='middle'):
            homopolymers[index] =(run[0], len(run))
            
        tmpf_h = P.getTempFilename('.')
        with open(tmpf_h, 'w') as outf:
            outf.write('index,base,length\n')
            for i in range(len(strain_seq)):
                if i in homopolymers.keys():
                    base, length = homopolymers[i]
                else:
                    base, length = '', '0'
                outf.write(','.join(map(str, [i, base, length])) + '\n')

        min_run_size = str(PARAMS['cross_match_homopolymer_length'])

        R('''
         require("ggplot2")
         # fetch the matrix of deletions
         df.del <- read.csv('{tmpf_d}',
                            header=FALSE,
                            stringsAsFactors=FALSE,
                            row.names=1)
         n.reads <- ncol(df.del)
         # summarize the number of deletions per position
         df.del <- rowSums(df.del)
         df.del <- as.data.frame(df.del)
         df.del <- df.del/n.reads*100
         # add upper panel and colour columns 
         df.del$panel = c('upper')
         df.del$base = c('deletion')
         # fix column names 
         colnames(df.del) <- c('height', 'panel', 'base')
         df.del$index <- rownames(df.del)

         # fetch the homopolymer runs
         df.hom <- read.csv('{tmpf_h}',
                            header=TRUE,
                            row.names=1,
                            stringsAsFactors=FALSE)
         # add lower panel and fix column names
         df.hom$panel = c('lower')
         colnames(df.hom) <- c('base', 'height', 'panel')
         df.hom <- df.hom[,c('height', 'panel', 'base')]
         df.hom$index <- rownames(df.hom)

         # remove homopolymer runs of < specified distance
         df.hom$height[df.hom$height < {min_run_size}] <- 0

         # stack
         df <- rbind(df.hom, df.del)
         # convert lower panel to negative
         df$height[df$panel=="lower"] <- -df$height[df$panel=="lower"]

         # assumes plot order
         cols <- c(NA, "dodgerblue4", "red3",
                   "lavenderblush4", "green2", "yellow2")

         #plot 
         pl <- ggplot(df, aes(x=as.numeric(index), 
                              y=height,
                              colour=base,
                              fill=base)) + geom_bar(stat="identity")
         pl <- pl + scale_colour_manual(values=cols) + scale_fill_manual(values=cols)
         pl <- pl + theme_bw() + theme(legend.position = 'none') + xlab('') + ylab('% Deletions')
         pl <- pl + ggtitle('{strain_name}')
         pdf('{outfile}')
         plot(pl, height=2.5, width=2.5)
         dev.off()
         '''.format(**locals()))
        L.info('Plots completed for strain {}'.format(strain))

@follows(plotErrorHistograms,
         plotDeletionsVsHomopolymerRuns)
def plotCM():
    pass


###############################################################################
## Task: Find read accuracy vs. number of passes
###############################################################################
###############################################################################
### Subtask: Parse CCS output for number of passs
###############################################################################
@follows(mkdir('mock_error_vs_passes.dir'))
@transform(PARAMS['data_ccs_bams'].split(';'),
           regex('.+/(.+).bam'),
           r'mock_error_vs_passes.dir/\1.tsv.gz')
def fetchMockCCSPassCounts(infile, outfile):
    '''For each read in the CCS bam files, output the number of passes'''

    samfile = pysam.AlignmentFile(infile, 'rb', check_sq=False)

    # a list of sequences for which there was no np tag
    out_failed = IOTools.openFile(P.snip(outfile, '.tsv.gz') + '_failed.tsv.gz',
                                  'w')
    with IOTools.openFile(outfile, 'w') as outf:
        # iterate over all reads in sam file, fetch # passes and read name
        for sam in samfile.fetch(until_eof=True):
            read_name = sam.qname
            try:
                n_passes = sam.get_tag('np')
            except KeyError:
                out_failed.write(read_name + '\n')
                continue
            outf.write(read_name + '\t'
                       + str(n_passes) + '\t'
                       + str(len(sam.query_sequence)) + '\n')

    out_failed.close()

            
@merge(fetchMockCCSPassCounts, 'mock_error_passes.load')
def loadMockCCSPassCounts(infiles, outfile):
    '''Collate the pass counts for each CCS read and load'''

    # Check that there are no duplicates...
    read_names = []
    for infile in infiles:
        for line in IOTools.openFile(infile):
            read_names.append(line.split()[0])
    assert len(read_names) == len(set(read_names)), 'Duplicate read names!'

    # load...
    table_name = P.snip(outfile, '.load', strip_path=True)
    infiles = ' '.join(infiles)
    statement = ("zcat {infiles} |"
                 " python {scriptsdir}/csv2db.py"
                 "  --table={table_name}"
                 "  --header-names=read_name,number_of_passes,read_length"
                 "  --add-index=read_name"
                 " > {outfile}")
    to_cluster = False

    P.run()


###############################################################################
### Subtask: Calculate the relationship between # passes and error rate... 
###############################################################################
@jobs_limit(1)
@transform(loadCrossMatchSummary,
           regex('(.+).load'),
           add_inputs(loadMockCCSPassCounts),
           r'mock_error_vs_passes.dir/\1_pass_summary.tsv.gz')
def summarizeMockErrorsVsPasses(infiles, outfile):
    '''Collate information on the number of errors vs. number of passes
    for the mock community data'''

    error_table, pass_table = [P.snip(x, '.load', strip_path=True) for x in infiles]
    
    statement = ("SELECT"
                 "  a.QuerySeqID,"
                 "  a.QueryStart,"
                 "  a.QueryEnd,"
                 "  a.QueryTail,"
                 "  a.TargetSeqID,"
                 "  a.TargetStart,"
                 "  a.TargetEnd,"
                 "  a.TargetTail,"
                 "  a.INS,"
                 "  a.DEL,"
                 "  a.SNP,"
                 "  b.number_of_passes,"
                 "  b.read_length"
                 " FROM {error_table} AS a"
                 " INNER JOIN {pass_table} AS b"
                 " ON a.QuerySeqID == b.read_name".format(**locals()))
    df = pd.read_sql(statement, connect(), index_col='QuerySeqID')
    df.to_csv(IOTools.openFile(outfile, 'w'), sep='\t')


###############################################################################
# Metatask: ANALYSIS OF IN VIVO DATASETS
###############################################################################
## Task: Filter in vivo data
###############################################################################
### Subtask: Process output of ccs pipeline 
###############################################################################
# samples J036[30-40] are stool samples J036[70-80] are saliva
@follows(mkdir('amp_cleaned.dir'))
@transform(os.path.join(PARAMS['amp_ccs_bam_dir'], '*.bam'),
           regex('.+/(J.+).bam'),
           r'amp_cleaned.dir/\1.bam')
def sortAMPBams(infile, outfile):
    '''merge output of the ccs pipeline'''

    prefix = P.snip(outfile, '.bam')

    cluster_options = '-l walltime=00:20:00'
    statement = ("samtools sort"
                 " -o {outfile}"
                 " -T {prefix}"
                 " -O bam"
                 " {infile}")
    P.run()


@follows(mkdir('mouse_cleaned.dir'))
@collate(os.path.join(PARAMS['mouse_ccs_bam_dir'], '*.bam'),
         regex('.+/(J.+?)(?:_2)?.bam'),
         r'mouse_cleaned.dir/\1.bam')
def mergeMouseBams(infiles, outfile):
    '''merge output of the ccs pipeline'''

    # sort the bamfiles
    tempfiles = []
    for f in infiles:
        tmpf = P.getTempFilename('.')
        cluster_options = '-l walltime=00:20:00'
        statement = ("samtools sort"
                     " -T {tmpf}"
                     " -o {tmpf}.bam"
                     " -O bam"
                     " {f}")
        P.run()
        tempfiles.append(tmpf)

    # merge sorted bamfiles, note that header is taken from in1.bam
    tempfiles = [x + '.bam' for x in tempfiles]
    tmpfiles = " ".join(tempfiles)
    cluster_options = '-l walltime=00:20:00'
    statement = ("samtools merge -f {outfile} {tmpfiles};"
                 " samtools index {outfile}")
    P.run()

    # remove the temporary files
    for f in tempfiles:
        os.unlink(f)
        os.unlink(P.snip(f, '.bam'))


@transform([sortAMPBams, mergeMouseBams], suffix('.bam'), '.fasta.gz')
def convertInVivoBamsToFasta(infile, outfile):
    '''Convert bam files to fasta files'''
    cluster_options = '-l walltime=00:10:00'
    statement = ("samtools bam2fq -n {infile} |"
                 " awk 'NR % 4==1 {{print $0}} NR % 4==2 {{print $0}}' |"
                 " sed 's/@/>/' |"
                 " gzip > {outfile}")
    to_cluster = False
    P.run()


@collate(convertInVivoBamsToFasta,
         regex('(.+)/.+fasta.gz'),
         r'\1/sample_read_depth.tsv')
def countInVivoSampleReadDepth(infiles, outfile):

    with IOTools.openFile(outfile,'w') as outf:
        outf.write('SampleID\tTotalReads\n')
        for infile in infiles:
            name = P.snip(infile, '.fasta.gz', strip_path=True)
            n_reads = 0
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
                n_reads += 1
            outf.write(name + '\t' + str(n_reads) + '\n')

###############################################################################
### Subtask: Remove adapter sequences from in vivo datasets  using cutadapt
###############################################################################
@transform(convertInVivoBamsToFasta,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/\2_trimmed.fasta.gz')
def runInVivoCutAdapt(input_fasta, outfile):
    '''Trim primers from PacBio reads using cutadapt'''
    
    f_adapter = PARAMS['cutadapt_forward_adapter']
    r_adapter = PARAMS['cutadapt_reverse_adapter']
    r_adapter = P16S.reverseComplement(r_adapter)

    # create outfiles for different stages of adapter trimming...
    # info files to contain details of primer match location
    info_forward = P.snip(outfile, '.fasta.gz') + '_info_fwd.txt'
    info_forward_rc = P.snip(info_forward, '.txt') + '_rc.txt'
    info_reverse = P.snip(info_forward, '_fwd.txt') + '_rev.txt'
    # outfile to contain reads with trimmed forward primer
    out_forward = P.snip(outfile, '.fasta.gz') + '_forward.fasta.gz'
    # outfile to contain those reads that failed to find primers
    failed_fwd = P.snip(outfile, '_trimmed.fasta.gz') + '_failed_fwd.fasta'
    failed_rev = P.snip(outfile, '_trimmed.fasta.gz') + '_failed_rev.fasta'

    # temporary files for the untrimmed output in first step
    tmpf_rejected = P.getTempFilename('.') + '.fasta'
    tmpf_rc = P.getTempFilename('.') + '.fasta'

    # various run parameters
    f_length = len(f_adapter)
    r_length = len(r_adapter)
    length = min([r_length, f_length])
    error_rate = float(PARAMS['cutadapt_errors'])/length
    params = ("--error-rate={}"
              " --match-read-wildcards".format(error_rate))
    

    # Run cutadapt to trim the forward primer. 
    to_cluster = False
    statement = ("cutadapt"
                 " -g forward={f_adapter}"
                 " --info-file={info_forward}"
                 " --untrimmed-output={tmpf_rejected}"
                 " {params}"
                 " {input_fasta} 2> {outfile}.1.log |"
                 " gzip > {out_forward}".format(**locals()))
    P.run()

    # reverse complement sequences for which adapter wasn't found. 
    with IOTools.openFile(tmpf_rc, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(tmpf_rejected)):
            outf.write('>' + fasta.title + '\n')
            outf.write(P16S.reverseComplement(fasta.sequence) + '\n')

    # worried about read/write speed...
    assert os.path.exists(tmpf_rc), 'File is missing %s' % tmpf_rc

    # rerun cutadapt to trim the forward primer from the reverse complement
    to_cluster = False
    statement = ("cutadapt"
                 "  -g forward={f_adapter}"
                 "  --info-file={info_forward_rc}"
                 "  --untrimmed-output={failed_fwd}"
                 "  {params}"
                 "  {tmpf_rc} 2> {outfile}.2.log |"
                 " gzip >> {out_forward};"
#                 " checkpoint;" ??
                 " cutadapt"
                 "  -a reverse={r_adapter}"
                 "  --info-file={info_reverse}"
                 "  --untrimmed-output={failed_rev}"
                 "  {params}"
                 "  {out_forward} 2> {outfile}.3.log |"
                 " gzip > {outfile}".format(**locals()))
    P.run()

    os.unlink(tmpf_rejected)
    os.unlink(tmpf_rc)


@transform(runInVivoCutAdapt, suffix('.fasta.gz'), '_filtered.fasta.gz')
def filterInVivoReadsByLength(infile, outfile):
    '''Filter cutadapt cleaned fasta file to remove reads over and under
    a specified length
    '''

    min_len = PARAMS['filter_min_length']
    max_len = PARAMS['filter_max_length']

    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            if len(fasta.sequence) >= min_len and len(fasta.sequence) <= max_len:
                outf.write( '\n'.join([ '>' + fasta.title,
                                       fasta.sequence]) + '\n')
            else:
                continue


@collate(filterInVivoReadsByLength,
         regex('(.+)/.+_trimmed_filtered.fasta.gz'),
         r'\1/sample_filtering_summary.tsv')
def summarizeInVivoFiltering(infiles, outfile):

    def _get_n(inf):
        n = 0
        for line in FastaIterator.FastaIterator(IOTools.openFile(inf)):
            n += 1
        return n

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('SampleID\tFilteredReads\tTrimmedReads\tFailedTrim\tFailedFwd\tFailedRv\n')

        for infile in infiles:
            name = os.path.basename(infile)
            name = name.split('_')[0]
            
            filtered_fasta = infile
            trimmed_fasta = P.snip(infile, '_filtered.fasta.gz') + '.fasta.gz'
            failed_fwd_fasta = P.snip(trimmed_fasta, '_trimmed.fasta.gz') + '_failed_fwd.fasta'
            failed_rev_fasta = P.snip(trimmed_fasta, '_trimmed.fasta.gz') + '_failed_rev.fasta'

            filtered = _get_n(infile)
            trimmed = _get_n(trimmed_fasta)
            failed_trimming = trimmed - filtered
            failed_fwd = _get_n(failed_fwd_fasta)
            failed_rev = _get_n(failed_rev_fasta)

            outf.write('\t'.join([name,
                                  str(filtered),
                                  str(trimmed),
                                  str(failed_trimming),
                                  str(failed_fwd),
                                  str(failed_rev)]) + '\n')
            

@collate([summarizeInVivoFiltering, countInVivoSampleReadDepth],
         regex('(.+)/.+tsv'),
         r'\1/sample_processing_summary.tsv')
def summarizeInVivoSampleProcessing(infiles, outfile):
    infiles = ' '.join(infiles)

    statement = ("python {scriptsdir}/combine_tables.py"
                 " --log={outfile}.log"
                 " {infiles}"
                 " > {outfile}")
    to_cluster = False
    P.run()

            
@transform(filterInVivoReadsByLength,
           suffix('.fasta.gz'),
           '_dereplicated.fasta.gz')
def dereplicateInVivoReads(infile, outfile):
    '''Derepilcate sequences that appear in all amplicon subsets'''

    outf_tmp = P.snip(outfile, '.gz')

    tmpf = P.getTempFilename('.')
    cluster_options = '-l walltime=00:20:00'
    statement = ("gunzip -c {infile} > {tmpf};"
                 " usearch"
                 "  -derep_fulllength {tmpf}"
                 "  -fastaout {outf_tmp}"
                 "  -sizeout;"
                 " rm {tmpf};"
                 " gzip {outf_tmp}")
    P.run()


@collate(filterInVivoReadsByLength,
         regex('(.+)_cleaned.dir/.+_trimmed_filtered.fasta.gz'),
         r'\1_cleaned.dir/\1_pooled.fasta.gz')
def poolInVivoReads(infiles, outfile):
    '''Pool reads across samples'''
    infiles = ' '.join(infiles)
    to_cluster = False
    statement = "zcat {infiles} | gzip > {outfile}"
    P.run()


@transform(poolInVivoReads,
           suffix('.fasta.gz'),
           '_dereplicated.fasta.gz')
def dereplicatePooledReads(infile, outfile):
    outf_tmp = P.snip(outfile, '.gz')

    tmpf = P.getTempFilename('.')
    cluster_options = '-l walltime=00:20:00'
    statement = ("gunzip -c {infile} > {tmpf};"
                 " usearch"
                 "  -derep_fulllength {tmpf}"
                 "  -fastaout {outf_tmp}"
                 "  -sizeout;"
                 " rm {tmpf};"
                 " gzip {outf_tmp}")
    P.run()



###############################################################################
### Subtask: Extract the AMP V1-V3 stool samples
###############################################################################
V1_V3_Targets =PARAMS['in_vivo_v1v3_fastas'].split(',')

@follows(mkdir('ampV1V3_cleaned.dir'))
@transform(V1_V3_Targets,
           regex('.+/AMP1-(.+?)-.+.clean.fasta'),
           r'ampV1V3_cleaned.dir/\1_trimmed_filtered.fasta.gz')
def fetchAMPV1_V3Reads(infile, outfile):
    '''Retreive the V1-V3 Illumina sequence data'''
    
    statement = "cat {infile} | gzip > {outfile}"
    P.run()


###############################################################################
### Subtask: Downsample AMP V1-V3 or AMP V1-V9 reads
###############################################################################
@follows(fetchAMPV1_V3Reads)
@transform(filterInVivoReadsByLength,
           regex('(.+)/(J03630|J03631|J03632|J03644)_trimmed_filtered.fasta.gz'),
           add_inputs(r'ampV1V3_cleaned.dir/\2_trimmed_filtered.fasta.gz'),
           r'\1/\2ds_trimmed_filtered.fasta.gz')
def downsampleAMPReads_v1v9(infiles, outfile):
    v9, v3 = infiles

    n_v9 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9))])
    n_v3 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3))])

    if n_v9 > n_v3:
        x = 0
        pptn = float(n_v3)/float(n_v9)
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9)):
                n = random.random()
                if n <= pptn:
                    x += 1
                    outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
                else:
                    continue
        L.warn("For sample %s there are %i V1-V9 reads %i V1-V3 reads "
               "and %i reads in the output V9 file" % (os.path.basename(v9), n_v9, n_v3, x))
    else:
        shutil.copyfile(v9, outfile)


@follows(filterInVivoReadsByLength)
@transform(fetchAMPV1_V3Reads,
           regex('(.+)/(.+)_trimmed_filtered.fasta.gz'),
           add_inputs(r'amp_cleaned.dir/\2_trimmed_filtered.fasta.gz'),
           r'\1/\2ds_trimmed_filtered.fasta.gz')
def downsampleAMPReads_v1v3(infiles, outfile):
    v3, v9 = infiles

    n_v3 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3))])
    n_v9 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9))])

    if n_v3 > n_v9:
        x = 0
        pptn = float(n_v9)/float(n_v3)
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3)):
                n = random.random()
                if n <= pptn:
                    x += 1
                    outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
                else:
                    continue
        L.warn("For sample %s there are %i V1-V3 reads %i V1-V9 reads "
               "and %i reads in the output V3 file" % (os.path.basename(v3), n_v3, n_v9, x))
    else:
        shutil.copyfile(v3, outfile)


@merge([downsampleAMPReads_v1v3, downsampleAMPReads_v1v9],
       'amp_cleaned.dir/downsampling_summary.tsv')
def summarizeDownsampledReadDepths(infiles, outfile):
    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('SampleID\tReadDepth\n')
        for infile in infiles:
            name_prefix = os.path.basename(os.path.dirname(infile)).split('_')[0]
            name_suffix = os.path.basename(infile).split('_')[0]
            name = name_prefix + '_' + name_suffix

            n = 0
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
                n += 1
            outf.write(name + '\t' + str(n) + '\n')

###############################################################################
## Task: generate genus level OTUs for AMP1 and mouse stool data
###############################################################################
#sort dereplicated fasta files by size.
@follows(mkdir('amp_genera_otus.dir'))
@follows(mkdir('ampV1V3_genera_otus.dir'))
@follows(mkdir('mouse_genera_otus.dir'))
@collate([filterInVivoReadsByLength, fetchAMPV1_V3Reads],
         regex('(.+)_cleaned.dir/(J03630|J03631|J03632|J03644|J00916|J00925|J00937)'
               '_trimmed_filtered.fasta.gz'),
         r'\1_genera_otus.dir/\1_pooled.fasta.gz')
def poolAMPReads_usearch(infiles, outfile):
    '''Pool the V1-V9 reads in a manner compatible with USEARCH'''
    outf = IOTools.openFile(outfile, 'w')
    for infile in infiles:
        sample_id = P.snip(os.path.basename(infile), '_trimmed_filtered.fasta.gz')
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            title = fasta.title + ';barcodelabel=' + sample_id
            outf.write('>' + title + '\n' + fasta.sequence + '\n')
    outf.close()


@follows(mkdir('ampds_genera_otus.dir'))
@follows(mkdir('ampV1V3ds_genera_otus.dir'))
@collate([downsampleAMPReads_v1v9, downsampleAMPReads_v1v3],
         regex('(.+)_cleaned.dir/(J03630ds|J03631ds|J03632ds|J03644ds)'
               '_trimmed_filtered.fasta.gz'),
         r'\1ds_genera_otus.dir/\1ds_pooled.fasta.gz')
def poolDownsampledAMPReads_usearch(infiles, outfile):
    '''Pool the V1-V9 reads in a manner compatible with USEARCH'''
    outf = IOTools.openFile(outfile, 'w')
    for infile in infiles:
        sample_id = P.snip(os.path.basename(infile), '_trimmed_filtered.fasta.gz')
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            title = fasta.title + ';barcodelabel=' + sample_id
            outf.write('>' + title + '\n' + fasta.sequence + '\n')
    outf.close()


@transform([poolAMPReads_usearch, poolDownsampledAMPReads_usearch],
            suffix('.fasta.gz'), '_dereplicated.fasta.gz')
def dereplicateAMPReads(infile, outfile):
    '''Dereplicate the AMPV1-V3 reads'''
    outf_tmp = P.snip(outfile, '.gz')

    tmpf = P.getTempFilename('.')
    cluster_options = '-l walltime=00:20:00'
    statement = ("gunzip -c {infile} > {tmpf};"
                 " usearch"
                 "  -derep_fulllength {tmpf}"
                 "  -fastaout {outf_tmp}"
                 "  -sizeout;"
                 " rm {tmpf};"
                 " gzip {outf_tmp}")
    P.run()

    
@transform(dereplicateAMPReads,
           suffix('.fasta.gz'),
           '_sorted.fasta.gz')
def sortAMPReads(infile, outfile):
    '''Sort the faecal reads for the AMP V1-V3 and V1-V9 samples'''
    tmpf = P.getTempFilename('.')
    outf_tmpf = P.snip(outfile, '.gz')
    cluster_options = '-l walltime=00:20:00'
    statement = ("gunzip -c {infile} > {tmpf};"

                 " usearch"
                 "  -sortbysize {tmpf}"
                 "  -fastaout {outf_tmpf}"
                 "  -minsize {in_vivo_min_reps};"
                 " gzip {outf_tmpf}")
    P.run()


@subdivide(sortAMPReads,
           regex('(.+)/(.+)_pooled_dereplicated_sorted.fasta.gz'),
           r'\1/distance_*/\2.fasta')
def generateAMPOTUSubdirectories(infile, outfiles):
    '''Generate subdirectories for the AMP OTU clustering thresholds'''
    directories = PARAMS['in_vivo_otu_dist'].split(',')

    for directory in directories:
        out_dir = os.path.abspath(os.path.dirname(infile))
        out_dir = os.path.join(out_dir, 'distance_' + directory)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        outfile = os.path.join(out_dir,
                               P.snip(os.path.basename(infile),
                                      '_pooled_dereplicated_sorted.fasta.gz')
                               + '.fasta')
        
        statement = " zcat {infile} > {outfile}"
        to_cluster = False
        P.run()


@transform(generateAMPOTUSubdirectories, suffix('.fasta'), '_otus.fasta')
def makeAMPOTUs(infile, outfile):
    '''Create OTUs using usearch'''
    out_up = P.snip(outfile, '.fasta') + '_up.txt'   
    rd_pc = os.path.dirname(infile)[-3:]

    statement = ("usearch"
                 " -cluster_otus {infile}"
                 " -otus {outfile}"
                 " -uparseout {out_up}"
                 " -otu_radius_pct {rd_pc}"
                 " -relabel OTU_")
    cluster_options = '-l nodes=1:ppn=4'
    P.run()


###############################################################################
### Subtask: Reassign reads to AMP OTUs
###############################################################################
@transform(makeAMPOTUs,
           regex('(.+)_genera_otus.dir/(.+)/(.+)_otus.fasta'),
           add_inputs(r'\1_genera_otus.dir/\1_pooled.fasta.gz'),
           r'\1_genera_otus.dir/\2/\3_\2_otu_assignment.txt')
#           regex('amp(.*)_genera_otus.dir/(.+)/(.+)_otus.fasta'),
#           add_inputs(r'amp\1_genera_otus.dir/amp\1_pooled.fasta.gz'),
#           r'amp\1_genera_otus.dir/\2/\3_\2_otu_assignment.txt')
def assignAMPReadsToOTUs(infiles, outfile):
    'Re-assign all reads to OTUs at the specified similarity threshold'
    otu_file, fasta_file = infiles 
    tmpf = P.getTempFilename('.')

    sim = float(os.path.dirname(outfile)[-3:])
    sim = (100 - sim) / 100.0
    
    statement = (" zcat {fasta_file} |"
                 # its necessary to remove the ' ' in fasta headers
                 " sed 's/ /_/g' > {tmpf};"
                 " usearch "
                 "  -usearch_global {tmpf}"
                 "  -db {otu_file}"
                 "  -strand plus"
                 "  -id {sim}"
                 "  -uc {outfile}"
                 " > {outfile}")
    cluster_options = '-l nodes=1:ppn=4,walltime=00:20:00'
    P.run()


@transform(assignAMPReadsToOTUs,
           suffix('_assignment.txt'),
           '_countTable.txt')
def countAMPOTUAssignment(infile, outfile):
    '''Use RE scripts to generate OTU count tables'''
    
    statement = ("python {in_vivo_uparse_scriptsdir}/uc2otutab.py"
                 " {infile} > {outfile}")
    cluster_options = '-l walltime=00:20:00'
    P.run()


###############################################################################
### Subtask: Parse the greengenes databases into a form for use with mothur
###############################################################################
@follows(mkdir('mouse_classified.dir'))
@follows(mkdir('amp_classified.dir'))
@follows(mkdir('ampV1V3_classified.dir'))
@follows(mkdir('ampds_classified.dir'))
@follows(mkdir('ampV1V3ds_classified.dir'))
@split(PARAMS['in_vivo_ref_taxonomy'],
       ['mouse_classified.dir/99_otu_taxonomy.txt',
        'amp_classified.dir/99_otu_taxonomy.txt', 
        'ampV1V3_classified.dir/99_otu_taxonomy.txt',
        'ampds_classified.dir/99_otu_taxonomy.txt', 
        'ampV1V3ds_classified.dir/99_otu_taxonomy.txt'])
def parseInVivoTaxonomyFiles(infile, outfiles):
    '''Parse greengenes into an acceptable format for use with mothur's
    RDP classifier
    '''
    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

    for outfile in outfiles:
        statement = ("cat {infile} |"
                     " python {tax_parser}"
                     "  --log={outfile}.log"
                     "  --tax-split='; '"
                     "  --drop-empty"
                     "  -x ID"
                     "  --remove-whitespace"
                     "  --clean-greengenes"
                     " > {outfile}")
        to_cluster = False
        P.run()


@follows(parseInVivoTaxonomyFiles)
@split(PARAMS['in_vivo_ref_sequence'],
       ['mouse_classified.dir/99_otu.fasta',
        'amp_classified.dir/99_otu.fasta', 
        'ampV1V3_classified.dir/99_otu.fasta',
        'ampds_classified.dir/99_otu.fasta', 
        'ampV1V3ds_classified.dir/99_otu.fasta'])
def parseInVivoFastaFiles(infile, outfiles):
    '''Parse greengenes fasta file'''

    for outfile in outfiles:
        tax_file = os.path.dirname(outfile) + '/99_otu_taxonomy.txt'
        # get a list of the sequence IDs to be retained
        IDs = []
        for line in IOTools.openFile(tax_file):
            IDs.append(line.split()[0])

        # parse the fasta file
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
                if 'ID' + fasta.title in IDs:
                    outf.write('>ID' + fasta.title + '\n' + fasta.sequence + '\n')
                else:
                    continue

###############################################################################
### Subtask: Assign taxonomy to AMP OTUs
###############################################################################
@follows(parseInVivoTaxonomyFiles)
@follows(parseInVivoFastaFiles)
@transform(makeAMPOTUs,
           regex('amp(.*)_genera_otus.dir/(.+)/(.+)_otus.fasta'),
           add_inputs(r'amp\1_classified.dir/99_otu.fasta',
                      r'amp\1_classified.dir/99_otu_taxonomy.txt'),
           r'amp\1_genera_otus.dir/\2/\3_otus.99_otu_taxonomy.wang.taxonomy')
def runAMPMothurRDP(infiles, outfile):
    '''Classify AMP OTUs against greengenes'''
    input_fasta, ref_fasta, ref_taxonomy = infiles

    ref_fasta = os.path.abspath(ref_fasta)
    ref_taxonomy = os.path.abspath(ref_taxonomy)
    input_fasta=os.path.abspath(input_fasta)
    
    # for safety keep parallel mothur runs in separate directories
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))
    working_dir = os.path.dirname(outfile)

    statement = ("cd {working_dir};"
                 " ln -s {ref_fasta} .;"
                 " ln -s {ref_taxonomy} .;"
                 " mothur \"#classify.seqs("
                 "fasta={input_fasta},"
                 " template=99_otu.fasta,"
                 " taxonomy=99_otu_taxonomy.txt,"
                 " cutoff={in_vivo_cut_off},"
                 " iters=75,"
                 " processors={in_vivo_mothur_processors})\""
                 " cd -")

    cluster_options = '-l walltime=08:00:00,mem=12Gb,nodes=1:ppn={}'.format(PARAMS['in_vivo_mothur_processors'])
    to_cluster = False
    P.run()


@transform(runAMPMothurRDP,
           regex('(.+)/(.+)/(.+)_otus.99_otu_taxonomy.wang.taxonomy'),
           r'\1/\2/\3_\2_otu_taxonomy.txt')
def fetchAMPReadTaxonomy(infile, outfile):
    '''Parse output of mothur taxonomic assignment'''

    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

    statement = ("cat {infile} |"
                 " python {tax_parser}"
                 "  --log={outfile}.log"
                 "  --out-format=table"
                 "  --remove-confidence-thresholds"
                 " > {outfile}")
    cluster_options = '-l walltime=00:20:00'
    to_cluster = False
    P.run()


###############################################################################
### Subtask: Collate AMP genus counts with taxonomic identification
###############################################################################
@collate([fetchAMPReadTaxonomy, countAMPOTUAssignment],
         regex('(.+)/(.+)/(.+)(_countTable.txt|_taxonomy.txt)'),
         r'\1/\2/\3_genus_counts.txt')
def combineAMPReadCountsWithTaxonomy(infiles, outfile):
    '''Combine the taxonomic assignment and count tables, then take just
    the genera counts'''
    infiles = " ".join(infiles)

    statement = ("python {scriptsdir}/combine_tables.py"
                 "  --log={outfile}.log"
                 "  {infiles} "
#                 " cut -f7,9,10,11,12" # the genus and the four samples
                 " > {outfile}")
    P.run()


@collate(combineAMPReadCountsWithTaxonomy,
         regex('(.+)/(.+)/(amp|ampV1V3)_(.+)_otu_genus_counts.txt'),
         r'amp_genera_otus.dir/\2/amp_\4_genus_correlations.txt')
def combineV1V3andV1V9GenusCounts(infiles, outfile):
    '''Join the genus counts for the V1-V3 and V1-V9 OTU analysis'''
    for f in infiles:
        if os.path.basename(f).startswith('ampV1V3_'):
            v3 = f
        elif os.path.basename(f).startswith('amp_'):
            v9 = f
        else:
            raise Exception('Unrecognized file')

    # load the V1-V3 data
    v3 = pd.read_table(v3)
    # select only Genus J03630 J03631 J03632 J03644
    v3 = v3[['Genus', 'J03630', 'J03631', 'J03632', 'J03644']]
    # get the sum of replicate rows
    v3 = v3.groupby('Genus').sum()
    
            
    # same for the V1-V9 data
    v9 = pd.read_table(v9)
    v9 = v9[['Genus', 'J03630', 'J03631', 'J03632', 'J03644']]
    v9 = v9.groupby('Genus').sum()        

    # perform full outer join on the tables
    df = v3.merge(v9, how='outer', left_index=True,
                  right_index=True, suffixes=('_V1_V3', '_V1_V9'))

    df.to_csv(outfile, sep='\t')


@collate(combineAMPReadCountsWithTaxonomy,
         regex('(.+)/(.+)/(ampds|ampV1V3ds)_(.+)_otu_genus_counts.txt'),
         r'ampds_genera_otus.dir/\2/ampds_\4_genus_correlations.txt')
def combineV1V3andV1V9DownsampledGenusCounts(infiles, outfile):
    '''Join the genus counts for the V1-V3 and V1-V9 OTU analysis'''
    for f in infiles:
        if os.path.basename(f).startswith('ampV1V3ds_'):
            v3 = f
        elif os.path.basename(f).startswith('ampds_'):
            v9 = f
        else:
            raise Exception('Unrecognized file')

    # load the V1-V3 data
    v3 = pd.read_table(v3)
    # select only Genus J03630ds J03631ds J03632ds J03644ds
    v3 = v3[['Genus', 'J03630ds', 'J03631ds', 'J03632ds', 'J03644ds']]
    # get the sum of replicate rows
    v3 = v3.groupby('Genus').sum()
    
            
    # same for the V1-V9 data
    v9 = pd.read_table(v9)
    v9 = v9[['Genus', 'J03630ds', 'J03631ds', 'J03632ds', 'J03644ds']]
    v9 = v9.groupby('Genus').sum()        

    # perform full outer join on the tables
    df = v3.merge(v9, how='outer', left_index=True,
                  right_index=True, suffixes=('_V1_V3', '_V1_V9'))

    df.to_csv(outfile, sep='\t')


###############################################################################
## Task: Search dereplicated reads against the greengenes database
###############################################################################
@transform(PARAMS['green_genes_full'],
           regex('.+/(.+).fasta.gz'),
           r'\1.udb')
def indexGreenGenes(infile, outfile):
    '''Convert GreenGenes to a UDB database file'''
    tmpf = P.snip(os.path.basename(outfile)) + '.fasta'
    statement = ("zcat {infile} > {tmpf};"
                 " sleep 2m;"
                 " usearch -makeudb_usearch {tmpf}"
                 "  -output {outfile}"
                 " -threads 20"
                 " &> {outfile}.log")
    cluster_options = '-l walltime=24:00:00,mem=60Gb,nodes=1:ppn=20'
    P.run()
    os.unlink(tmpf)


@follows(mkdir('ampV1V3_genera_otus.dir'), 
         mkdir('ampV1V3_genera_otus.dir/gg_blast'),
         mkdir('amp_genera_otus.dir'),
         mkdir('amp_genera_otus.dir/gg_blast'))
@subdivide(sortAMPReads,
           regex('(.+)/(.+)_pooled_dereplicated_sorted.fasta.gz'),
           r'\1/gg_blast/distance*/amplicons.fasta')
def makeUSEARCHDirectories(infile, outfiles):
    '''Create subdirectories for searching dereplicated sequences
    against greengenes'''
    working_dir = os.path.dirname(infile)

    for dist in PARAMS['in_vivo_otu_dist'].split(','):
        sub_dir = 'distance_' + dist
        sub_dir = os.path.join(working_dir, 'gg_blast', sub_dir)
        if not os.path.exists(sub_dir):
            os.mkdir(sub_dir)
        outfile = os.path.join(sub_dir, 'amplicons.fasta')
        statement = "gunzip -c {infile} > {outfile}"
        cluster_options = '-l walltime=00:01:00'
        P.run()


@transform(makeUSEARCHDirectories,
           suffix('.fasta'),
           add_inputs(indexGreenGenes),
           '_gg_matches.b6')
def searchAMPSequencesAgainstGreenGenes(infiles, outfile):
    '''Match AMP sequences against GreenGenes database using Usearch'''
    amplicon_fasta, ref_db = infiles

    # fetch distance
    dist = float(os.path.dirname(outfile).split('_').pop())
    dist = (100-dist)/100.0

    # output sam file
    sam_out = P.snip(outfile, '.b6') + '.sam'
    
    statement = ("usearch -usearch_local {amplicon_fasta}"
                 " -db {ref_db}"
                 " -id {dist}"
                 " -maxaccepts 1" # terminate as soon as one positive hit made
                 " -maxrejects 0"
                 " -strand both"
                 " -output_no_hits"
                 " -query_cov {dist}"
                 " -blast6out {outfile}"
                 " -samout {sam_out}"
                 " -threads 16"
                 " &> {outfile}.log")
    cluster_options = '-l walltime=24:00:00,mem=60Gb,nodes=1:ppn=16'
    P.run()

    
###############################################################################
## Task: Assign all in vivo reads to a particular genus
###############################################################################
### Subtask: Classify all in vivo reads using RDP
###############################################################################
### need the classifed V1V3 reads to enter at this point. 
@follows(fetchAMPV1_V3Reads)
@follows(parseInVivoTaxonomyFiles)
@follows(parseInVivoFastaFiles)
@transform([filterInVivoReadsByLength, fetchAMPV1_V3Reads],
           regex('(.+)_cleaned.dir/(.+)_trimmed_filtered.fasta.gz'),
           add_inputs([r'\1_classified.dir/99_otu.fasta', # reference fasta
                       r'\1_classified.dir/99_otu_taxonomy.txt']), # reference taxonomy
           r'\1_classified.dir/\2/\2.99_otu_taxonomy.wang.taxonomy')
def runInVivoMothurRDP(infiles, outfile):
    '''Classify in vivo reads against greengenes'''
    input_fasta, reference_files = infiles
    ref_fasta, ref_taxonomy = reference_files

    ref_fasta = os.path.abspath(ref_fasta)
    ref_taxonomy = os.path.abspath(ref_taxonomy)

    # for safety keep parallel mothur runs in separate directories
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))
    working_dir = os.path.dirname(outfile)

    fasta = P.snip(outfile, '.99_otu_taxonomy.wang.taxonomy') + '.fasta'
    fasta_file = os.path.basename(fasta)

    statement = (" zcat {input_fasta} > {fasta};"
                 " cd {working_dir};"
                 " ln -s {ref_fasta} .;"
                 " ln -s {ref_taxonomy} .;"
                 " mothur \"#classify.seqs("
                 "fasta={fasta_file},"
                 " template=99_otu.fasta,"
                 " taxonomy=99_otu_taxonomy.txt,"
                 " cutoff={in_vivo_cut_off},"
                 " iters=75,"
                 " processors={in_vivo_mothur_processors})\""
                 " cd -")

    cluster_options = '-l walltime=08:00:00,mem=12Gb,nodes=1:ppn={}'.format(PARAMS['in_vivo_mothur_processors'])
    to_cluster = True
    P.run()


@transform(runInVivoMothurRDP,
           suffix('.99_otu_taxonomy.wang.taxonomy'),
           '_taxonomic_assignment.tsv.gz')
def fetchInVivoReadTaxonomy(infile, outfile):
    '''Parse output of mothur taxonomic assignment'''

    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

    statement = ("cat {infile} |"
                 " python {tax_parser}"
                 "  --log={outfile}.log"
                 "  --out-format=table"
                 "  --remove-confidence-thresholds |"
                 " gzip > {outfile}")
    cluster_options = '-l walltime=00:20:00'
    P.run()


###############################################################################
### Subtask: Extract reads identified as belonging to genera of interest
###############################################################################
@follows(filterInVivoReadsByLength)
@subdivide(fetchInVivoReadTaxonomy,
           # taking only the stool samples
           regex('(.+)_classified.dir/(J03630|J03631|J03632|J03644)/.+_taxonomic_assignment.tsv.gz'),
           add_inputs(r'\1_cleaned.dir/\2_trimmed_filtered.fasta.gz'),
           r'\1_classified.dir/\2_*fasta.gz')
def filterGenusSpecificReads(infiles, outfiles):
    '''Filter those reads classified as belonging to genera of interest
    '''
    # a table containing read taxonomic classification
    # and a fasta file containing reads
    tax_table, fasta_file = infiles
    tax_df = pd.read_table(tax_table, index_col=0, compression='gzip')
    
    # outfile stub
    out_dir = os.path.dirname(os.path.dirname(tax_table))
    out_prefix = os.path.basename(os.path.dirname(tax_table))
    out_stub = os.path.join(out_dir, out_prefix + '_')
    
    # a list of those genera to be extracted
    genera_list = PARAMS['in_vivo_genera']
    genera_list = [x.lower() for x in genera_list.split(',')]

    for genus in genera_list:
        L.info('Processing genus %s for sample %s' % (genus, out_prefix))
        outfile = out_stub + genus + '.fasta.gz'
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
                # the V1-V3 data have ':' in title which mothur replaces with
                # underscores. Mothur also only takes first title element.
                title = fasta.title.split()[0]
                title = re.sub(':', '_', title)
                g = tax_df.loc[title, 'Genus']
                if g.lower() == genus:
                    outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
                else:
                    continue 


@transform(filterGenusSpecificReads,
           regex('ampV1V3_classified.dir/(.+)_(.+).fasta.gz'),
           add_inputs(r'amp_classified.dir/\1_\2.fasta.gz'),
           r'ampV1V3_classified.dir/\1ds_\2.fasta.gz')
def downSampleGenusSpecificReads_v1v3(infiles, outfile):
    '''Downsample the V1-V3 samples, if necessary'''
    v3, v9 = infiles

    n_v3 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3))])
    n_v9 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9))])

    if n_v3 > n_v9:
        x = 0
        pptn = float(n_v9)/float(n_v3)
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3)):
                n = random.random()
                if n <= pptn:
                    x += 1
                    outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
                else:
                    continue
        L.warn("For sample %s there are %i V1-V3 reads %i V1-V9 reads "
               "and %i reads in the output V3 file" % (os.path.basename(v3), n_v3, n_v9, x))
    else:
        shutil.copyfile(v3, outfile)


@transform(filterGenusSpecificReads,
           regex('amp_classified.dir/(.+)_(.+).fasta.gz'),
           add_inputs(r'ampV1V3_classified.dir/\1_\2.fasta.gz'),
           r'amp_classified.dir/\1ds_\2.fasta.gz')
def downSampleGenusSpecificReads_v1v9(infiles, outfile):
    '''Downsample the V1-V9 samples, if necessary'''
    v9, v3 = infiles

    n_v9 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9))])
    n_v3 = len([fasta.title for fasta in FastaIterator.FastaIterator(IOTools.openFile(v3))])

    if n_v9 > n_v3:
        x = 0
        pptn = float(n_v3)/float(n_v9)
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(v9)):
                n = random.random()
                if n <= pptn:
                    x += 1
                    outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
                else:
                    continue
        L.warn("For sample %s there are %i V1-V9 reads %i V1-V3 reads "
               "and %i reads in the output V9 file" % (os.path.basename(v9), n_v9, n_v3, x))
    else:
        shutil.copyfile(v9, outfile)

                
@collate([filterGenusSpecificReads,
          downSampleGenusSpecificReads_v1v3,
          downSampleGenusSpecificReads_v1v9],
         regex('amp(.*)_classified.dir/J\d{5}(.*)_(.+).fasta.gz'),
         r'amp\1_classified.dir/amp\1\2_\3.fasta')
def poolGenusReads(infiles, outfile):
    '''Pool the reads in a manner compatible with USEARCH'''
    outf = IOTools.openFile(outfile, 'w')
    for infile in infiles:
        sample_id = os.path.basename(infile)
        sample_id = sample_id.split('_')[0]
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
            title = fasta.title + ';barcodelabel=' + sample_id
            outf.write('>' + title + '\n' + fasta.sequence + '\n')
    outf.close()


@transform(poolGenusReads,
           suffix('.fasta'),
           '_dereplicated.fasta')
def dereplicateGenusReads(infile, outfile):
    cluster_options = '-l walltime=00:40:00'
    statement = (" usearch"
                 "  -derep_fulllength {infile}"
                 "  -fastaout {outfile}"
                 "  -sizeout")
    to_cluster = False
    P.run()


#################################################################################
## Task: Calculate species/strain level OTUs for each genus
#################################################################################
@follows(mkdir('amp_species_otus.dir'))
@follows(mkdir('ampV1V3_species_otus.dir'))
@transform(dereplicateGenusReads,
           regex('amp(.*)_classified.dir/amp(.*)_(.+)_dereplicated.fasta'),
           r'amp\1_species_otus.dir/amp\2_\3.dir/amp\2_\3_sorted.fasta')
def sortGenusAmpliconsBySize(infile, outfile):
    '''Sort the bacteroides sequencs by the number of duplicates discard 
    sequences with fewer than n replicates'''

    # create the subdirectory
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    cluster_options = '-l walltime=00:20:00'
    statement = ("usearch"
                 "  -sortbysize {infile}"
                 "  -fastaout {outfile}"
                  "  -minsize {in_vivo_min_reps}")
    to_cluster = False
    P.run()
    
           
@subdivide(sortGenusAmpliconsBySize,
           regex('(.+)/(.+).fasta'),
           r'\1/distance_*/\2.fasta')
def generateOTUSubdirectories(infile, outfiles):
    '''Generate subdirectories for OTU clustering using usearch'''

    directories = PARAMS['in_vivo_otu_dist'].split(',')

    for directory in directories:
        out_dir = os.path.abspath(os.path.dirname(infile))
        out_dir = os.path.join(out_dir, 'distance_' + directory)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        
        outfile = os.path.join(out_dir, os.path.basename(infile))
        
        statement = "cat {infile} > {outfile}"
        to_cluster = False
        P.run()
           

@transform(generateOTUSubdirectories,
           suffix('_sorted.fasta'),
           '_otus.fasta')
def makeGenusSpecificOTUs(infile, outfile):
    '''Create OTUs using usearch'''

    out_up = P.snip(outfile, '.fasta') + '_up.txt'   
    rd_pc = os.path.dirname(infile)[-3:]

    # a sanity check, after downsampling and removing reads with fewer
    # than three replicates some files are empty.
    n_reads = IOTools.openFile(infile).readlines()
    if len(n_reads) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("usearch"
                     " -cluster_otus {infile}"
                     " -otus {outfile}"
                     " -uparseout {out_up}"
                     " -otu_radius_pct {rd_pc}"
                     " -relabel OTU_")
        cluster_options = '-l walltime=00:20:00'
        to_cluster = False
        P.run()


@follows(poolGenusReads)
@transform(makeGenusSpecificOTUs,
           regex('amp(.*)_species_otus.dir/(.+)_(.+).dir/(.+)/(.+)_otus.fasta'),
           add_inputs(r'amp\1_classified.dir/\2_\3.fasta'),
           r'amp\1_species_otus.dir/\2_\3.dir/\4/\5_otu_assignment.txt')
def assignGeneraReadsToOTUs(infiles, outfile):
    '''Re-assign all reads to OTUs at a specified similarity threshold'''
    otu_file, fasta_file = infiles
    tmpf = P.getTempFilename('.')
    
    sim = float(os.path.dirname(outfile)[-3:])
    sim = (100 - sim) / 100.0

    # a sanity check, after downsampling and removing reads with fewer
    # than three replicates some files are empty.
    n_reads = IOTools.openFile(otu_file).readlines()
    if len(n_reads) == 0:
        open(outfile, 'w').close()
    else:
        statement = (" cat {fasta_file} |"
                     # its necessary to remove the ' ' in fasta headeres
                     " sed 's/ /_/g' > {tmpf};"
                     " usearch "
                     "  -usearch_global {tmpf}"
                     "  -db {otu_file}"
                     "  -strand plus"
                     "  -id {sim}"
                     "  -uc {outfile}"
                     " > {outfile}")
        cluster_options = '-l walltime=00:40:00'
        P.run()


@transform(assignGeneraReadsToOTUs,
           suffix('_assignment.txt'),
           '_countTable.txt')
def countGeneraOTUAssignment(infile, outfile):
    '''Use RE scripts to generate OTU count tables'''
    
    # a sanity check, after downsampling and removing reads with fewer
    # than three replicates some files are empty.
    n_reads = IOTools.openFile(infile).readlines()
    if len(n_reads) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("python {in_vivo_uparse_scriptsdir}/uc2otutab.py"
                     " {infile} > {outfile}")
        cluster_options = '-l walltime=00:20:00'
        P.run()


###############################################################################
## Subtask: Retrieve sequences that were successfully assigned to OTUs. 
###############################################################################
# parse the amp_stool_otus.dir/distance_1.0/*_otu_assignment.txt file
# if the read is assigned. then add read title to dict list, where key is barcode label
# parse the four bacteroides fasta files in amp_classified.dir, 
# if the fasta title is in the relevent dict list, then write sequence to outfile. 

@follows(filterGenusSpecificReads)
@transform(assignGeneraReadsToOTUs,
           regex('amp_species_otus.dir/amp(.*)_(.+).dir/(.+)/(.+)_assignment.txt'),
           # can't use glob
           add_inputs(r'amp_classified.dir/J03630\1_\2.fasta.gz',
                      r'amp_classified.dir/J03631\1_\2.fasta.gz',
                      r'amp_classified.dir/J03632\1_\2.fasta.gz',
                      r'amp_classified.dir/J03644\1_\2.fasta.gz'),
           r'amp_species_otus.dir/amp\1_\2.dir/\3/\4_assigned.fasta.gz')
def fetchOTUAssignedReads(infiles, outfile):
    '''Fetch reads that were successfully assigned to OTUs'''
    infile = infiles[0]
    fasta_files = infiles[1:]

    ff = os.path.basename(fasta_files[0]).split('_')[0]
    if ff.endswith('ds'):
        sample_dict = {'J03630ds': [], 'J03631ds': [], 'J03632ds': [], 'J03644ds': []}        
    else:    
        sample_dict = {'J03630': [], 'J03631': [], 'J03632': [], 'J03644': []}
    for line in open(infile):
        if line.endswith('*\n'):
            continue
        else:
            title = line.split()[8]
            title, sample = title.split(';')
            sample = sample.split('=')[1]
            sample_dict[sample].append(title)


    outf = IOTools.openFile(outfile, 'w')
    for fasta_file in fasta_files:
        fasta_id = os.path.basename(fasta_file).split('_')[0]

        if fasta_id in sample_dict.keys():
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
                if fasta.title.split()[0] in sample_dict[fasta_id]:
                    outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')
                else:
                    continue
        else:
            continue

        
###############################################################################
### Subtask: Align reads for each genus and create entropy plots
###############################################################################
@transform(fetchOTUAssignedReads,
           regex('(.+)/(.+)_otu_assigned.fasta.gz'),
           add_inputs(r'\1/\2_otus.fasta'),
           r'\1/\2_otu_assigned_base.fasta')
def addBaseOTUSequence(infiles, outfile):
    '''Add a single sequence representative of the first OTU'''
    full_fasta, otu_fasta = infiles

    # fetch the fasta sequence for the first OTU identified
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(otu_fasta)):
        base = fasta.sequence
        break

    if not base:
        IOTools.openFile(outfile, 'w').close()
    else:
        with IOTools.openFile(outfile, 'w') as outf:
            for fasta in FastaIterator.FastaIterator(IOTools.openFile(full_fasta)):
                outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
            outf.write('>BASE\n' + base + '\n')
        

@transform(addBaseOTUSequence, suffix('_base.fasta'), '_aligned.fasta')
def alignOTUAssignedReads(infile, outfile):
    '''Create muscle alignment of OTU assigned reads'''

    # There may be empty files that will cause muscle to fail
    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        statement = (" muscle -in {infile} -out {outfile} -maxhours 4"
                     " &> {outfile}.log")
        cluster_options = '-l walltime=04:00:00,mem=12Gb,nodes=1:ppn=4'
        P.run()

        
@transform(alignOTUAssignedReads, suffix('.fasta'), '_entropy.tsv')
def  calculateEntropyForEachGenus(infile, outfile):
    '''Calculate the shannon entropy for multiple sequence alignments'''

    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        entropy_scr = os.path.join(os.path.dirname(__file__),
                                   'alignment2entropy.py')
        statement = ("python {entropy_scr}"
                     " -a {infile}"
                     " -b BASE"
                     " -s 75"
                     " -o {outfile}")
        cluster_options = '-l walltime=00:20:00'
        P.run()


@jobs_limit(1, 'RGlobal')
@transform(calculateEntropyForEachGenus,
           regex('(.+)/(.+)/(.+)_otu_assigned_aligned_entropy.tsv'),
           r'\1/\3_\2.png')
def plotEntropyForEachGenus(infile, outfile):
    '''Plot the entropy for each alignment of reads belonging to a 
    particular genus. Intervals are added on to match primer coordinates
    on bacteroides vulgatus'''

    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:

        # fetch titles
        title = os.path.basename(infile).split('_')[1]
        
        R('''rm(list=ls())''')
        df = pd.read_table(infile)
        df_m = pd.melt(df, id_vars=['Base_Position',])
        R.assign('df', df_m)

        R('''
        require('ggplot2')
        require('scales')
        
        # normalize the entropy values
        df$value_n <- df$value / max(df$value, na.rm=T)
        
        # region title offset
        n = 0.01
        # region title alpha
        ra=0.8
        # region title colour
        rc='dark grey'
        # panel alpha
        pa=0.2
        #panel colour
        pc='grey'
        #amplicon colour
        ac = 'dark red'
        
        plt <- ggplot(df)
        plt <- plt + geom_line(aes(x=Base_Position, y=value_n, fill=NULL))
        plt <- plt + scale_colour_manual(values=c('black',))
        plt <- plt + theme_bw() + ylim(-0.05, 1.05)
        
        # add individual bars
        #V1
        plt <- plt + annotate('rect', xmin=31, xmax=110, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=31, xend=110, y=1, yend=1, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=9 + (110-31)/2, y=1+n, label='V1', cex=1)
        #V2
        plt <- plt + annotate('rect', xmin=126, xmax=258, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=126, xend=258, y=1+n, yend=1+n, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=126 + (258-126)/2, y=1+n+n, label='V2', cex=1)
        #V3
        plt <- plt + annotate('rect', xmin=365, xmax=519, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=365, xend=519, y=1, yend=1, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=365 + (519-365)/2, y=1+n, label='V3', cex=1)
        #V4
        plt <- plt + annotate('rect', xmin=537, xmax=689, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=537, xend=689, y=1+n, yend=1+n, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=537 + (689-537)/2, y=1+n+n, label='V4', cex=1)
        #V5
        plt <- plt + annotate('rect', xmin=801, xmax=910, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=801, xend=910, y=1, yend=1, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=801 + (910-801)/2, y=1+n, label='V5', cex=1)
        #V6
        plt <- plt + annotate('rect', xmin=986, xmax=1059, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=986, xend=1059, y=1+n, yend=1+n, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=986 + (1059-986)/2, y=1+n+n, label='V6', cex=1)
        #V7
        plt <- plt + annotate('rect', xmin=1114, xmax=1176, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=1114, xend=1176, y=1, yend=1, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=1114 + (1176-1114)/2, y=1+n, label='V7', cex=1)
        #V8
        plt <- plt + annotate('rect', xmin=1240, xmax=1369, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=1240, xend=1369, y=1+n, yend=1+n, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=1240 + (1369-1240)/2, y=1+n+n, label='V8', cex=1)
        #V9
        plt <- plt + annotate('rect', xmin=1394, xmax=1482, ymin=0, ymax=1, fill=pc, alpha=pa)
        plt <- plt + annotate('segment', x=1394, xend=1482, y=1, yend=1, colour=rc, alpha=ra)
        #plt <- plt + annotate('text', x=1394 + (1482-1394)/2, y=1+n, label='V9', cex=1)
        
        # Add amplicon regions
        #V1-V3
        plt <- plt + annotate('segment', x=31, xend=522, y=-0.01-n, yend=-0.01-n, colour=ac, alpha=ra)
        #plt <- plt + annotate('text', x=31 + (522-31)/2, y=-0.01-2*n, label='V1-V3', colour=ac, cex=3)
        #V1-V9
        plt <- plt + annotate('segment', x=31, xend=1482, y=-0.02-2*n, yend=-0.02-2*n, colour=ac, alpha=ra)
        #plt <- plt + annotate('text', x=31 + (1482-31)/2, y=-0.02-3*n, label='V1-V9', colour=ac, cex=3)
        
        plt <- plt + theme(panel.grid = element_blank())
        plt <- plt + theme(legend.position='none', axis.title=element_blank())
        plt <- plt + ggtitle('{title}')
        
        # plot to outfile
        png('{outfile}', width=2.5, height=2.0, units='in', res=300)
        par(mar=c(0,0,0,0))
        plot(plt)
        dev.off()
        '''.format(**locals()))
    

# This doesn't work because ImageMagick is only installed on the head node. 
# @jobs_limit(1)
# @collate(plotEntropyForEachGenus,
#          regex('(.+)/.+/amp(.*)_(.+)_distance_(.+).png'),
#          r'\1/amp\2_distance_\4.pdf')
# def collateEntropyPlots(infiles, outfile):
#     '''Merge entropy plots into a single pdf'''
#     infiles = " ".join(infiles)
#     statement = "convert {infiles} {outfile}"
#     to_cluster = False
#     P.run()

###############################################################################
# Metatask: COMPARE IN VIVO ANALYSIS RESULTS TO AMP1 DATA
###############################################################################
## Task: Summarize mWGS data
###############################################################################
### Subtask: Parse the mWGS data for species of interest
###############################################################################
@follows(mkdir('rtg_mwgs_data.dir'))
@subdivide(os.path.join(PARAMS['amp_rtg_mwgs'], '*.tsv'),
           regex('.+/(.+).tsv'),
           r'rtg_mwgs_data.dir/\1_*_str.tsv')
def parseRTGStrainOutput(infile, outfiles):
    '''Parse the RTG output into text format to retrieve strain abundances'''
    df = pd.read_table(infile)
    
    for genus in PARAMS['in_vivo_genera'].split(','):
        outfile = P.snip(infile, '.tsv', strip_path=True)
        outfile = outfile + '_' + genus + '_str.tsv'
        outfile = os.path.join('rtg_mwgs_data.dir', outfile)
        
        # NB. This only works if genera start with uppercase!
        df_s = df[df['rank'] == 'strain']
        df_s = df_s[df_s['taxonomy-name'].str.contains(genus)]
        df_s = df_s[['taxonomy-name', 'abundance']]
        df_s.to_csv(outfile, sep='\t', index=False)

        
@follows(mkdir('rtg_mwgs_data.dir'))
@subdivide(os.path.join(PARAMS['amp_rtg_mwgs'], '*.tsv'),
           regex('.+/(.+).tsv'),
           r'rtg_mwgs_data.dir/\1_*_species.tsv')
def parseRTGSpeciesOutput(infile, outfiles):
    '''Parse the RTG output into text format to retrieve strain abundances'''
    df = pd.read_table(infile)
    
    for genus in PARAMS['in_vivo_genera'].split(','):
        outfile = P.snip(infile, '.tsv', strip_path=True)
        outfile = outfile + '_' + genus + '_species.tsv'
        outfile = os.path.join('rtg_mwgs_data.dir', outfile)
        
        # NB. This only works if genera start with uppercase!
        df_s = df[df['rank'] == 'species']
        df_s = df_s[df_s['taxonomy-name'].str.contains(genus)]
        df_s = df_s[['taxonomy-name', 'abundance']]
        df_s.to_csv(outfile, sep='\t', index=False)


@collate([parseRTGSpeciesOutput, parseRTGStrainOutput],
         regex('(.+)/.+_(.+)_(str|species).tsv'),
         r'\1/\2_\3_abundance.tsv')
def collateRTGOutput(infiles, outfile):
    '''Combine mWGS data for different samples'''
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --log={outfile}.log"
                 " --add-file-prefix"
                 " --regex-filename='(.+).tsv'"
                 " {infiles} |"
                 " sed 's/-/_/g' |" # WARNING! This is affecting scientific numbers. 
                 " sed 's/_abundance//g' > {outfile}")
    cluster_options = '-l walltime=00:01:00'
    to_cluster = False
    P.run()


@transform(collateRTGOutput,
           regex('(.+)/(.+)_str_abundance.tsv'),
           add_inputs(r'\1/\2_species_abundance.tsv'),
           r'\1/\2_strain_abundance.tsv')
def addMissingSpeciesToRTGStrainData(infiles, outfile):
    '''There are unidentified species that don't appear in the
    strain-level data.
    '''
    strain_file, species_file = infiles
    df_str = pd.read_table(strain_file, index_col=0)
    df_sp = pd.read_table(species_file, index_col=0) 
    
    # a list of the strains in the strain level table
    strain_list = df_str.index.tolist()
    
    for row in df_sp.iterrows():
        # there are sometimes [] in the ncbi taxon ids.
        strain_list = [re.sub('\[|\]', '', x) for x in strain_list]
        name = re.sub('\[|\]', '', row[1].name)

        # check whether the species matches any of the strain names        
        present = [re.search(row[1].name, x) for x in strain_list]
        present = sum([x != None for x in present])
        if present:
            continue
        else:
            out_series = row[1]
            out_series.index = [re.sub('_species', '_str', x) for x in out_series.index]
            df_str = df_str.append(out_series)
            
    df_str.to_csv(outfile, sep='\t')


@transform(addMissingSpeciesToRTGStrainData,
           regex('(.+)/(.+).tsv'),
           r'\1/\2.load')
def loadRTGStrainOutput(infile, outfile):
    table_name = P.snip(infile, '.tsv', strip_path=True)
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 " --table={table_name}"
                 " --allow-empty-file"
                 " > {outfile}")
    cluster_options = '-l walltime=00:01:00'
    to_cluster = False
    P.run()


@transform(collateRTGOutput,
           regex('(.+)/(.+)_species_abundance.tsv'),
           r'\1/\2_species.load')
def loadRTGSpeciesOutput(infile, outfile):
    table_name = P.snip(infile, '.tsv', strip_path=True)
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 " --table={table_name}"
                 " --allow-empty-file"
                 " > {outfile}")
    cluster_options = '-l walltime=00:01:00'
    to_cluster = False
    P.run()


@follows(loadRTGStrainOutput,
         loadRTGSpeciesOutput)
def loadmWGSStrainData():
    pass


###############################################################################
## Task: Search in vivo OTUs against RTG database 
###############################################################################
### Subtask: Create a udb file from the RTG fasta
###############################################################################
@follows(mkdir('rtg_amp_comparison.dir'))
@files(PARAMS['data_rtg_fasta'], 'rtg_amp_comparison.dir/rtg_db.udb')
def indexRTGDb(infile, outfile):
    '''Convert RTG fasta into a UDB database file'''
    statement = ("usearch -makeudb_usearch {infile}"
                 " -output {outfile}"
                 " -threads 10")
    to_cluster = False
    P.run()


@follows(indexRTGDb,
         mkdir('amp_rtg_comparison.dir'),
         mkdir('ampV1V3_rtg_comparison.dir'))
@jobs_limit(1)
@transform(makeGenusSpecificOTUs,
           regex('amp(.*)_species_otus.dir/amp(.*)_(.+).dir/distance_(.+)/(.+)_otus.fasta'),
           r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/sentinel.txt')
def makeUSEARCHSubDirectories(infile, outfile):
    '''Make the subdirectories for USEARCH mapping to rtg database
    This is a hack to avoid errors when creating the same directory
    multiple times in parallel. 
    '''
    sub_dir = os.path.dirname(outfile)
    parent_dir = os.path.dirname(sub_dir)
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)
    open(outfile, 'w').close()


###############################################################################
### Subtask: Run USEARCH
###############################################################################
@follows(makeUSEARCHSubDirectories)
@transform(makeGenusSpecificOTUs,
           regex('amp(.*)_species_otus.dir/amp(.*)_(.+).dir/distance_(.+)/(.+)_(bacteroides|faecalibacterium)_otus.fasta'),
           add_inputs('rtg_amp_comparison.dir/rtg_db.udb'),
           r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_\6_otu_matches.b6')
def searchOTUsAgainstRTG(infiles, outfile):
    '''Match OTUs against RTG database using USEARCH'''
    # with the exception of amp_bacteroides/distance_0.1 this finishes in < 24hrs

    otu_fasta, reference_udb = infiles

    # fetch the distance at which to threshold 'hits'
    dist = float(os.path.dirname(outfile).split('_').pop())
    dist = (100 - dist)/100.0

    # output sam file
    sam_out = P.snip(outfile, '.b6') + '.sam'

    if len(IOTools.openFile(otu_fasta).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("usearch -usearch_local {otu_fasta}"
                     " -db {reference_udb}"
                     " -id {dist}"
                     " -maxaccepts 0"
                     " -maxrejects 0"
                     " -strand both"
                     " -top_hits_only"
                     " -maxhits 50"
                     " -query_cov {dist}"
                     " -blast6out {outfile}"
                     " -samout {sam_out}"
                     " -threads 16"
                     " &> {outfile}.log")
        cluster_options = '-l walltime=24:00:00,mem=110Gb,nodes=1:ppn=16'
        P.run()


@transform(searchOTUsAgainstRTG, suffix('.b6'), '_gi2taxid.tsv')
def getGi2TaxId(infile, outfile):
    '''For each of the gi numbers matching an OTU, fetch the taxonomic
    id from NCBI dump used to create the RTG fasta file. 

    The RTG fasta file is found in /data/weinstocklab/DBS/
    references-species-all-2015-04-10/references/species-all-2015-04-10-clean.fasta
    In the parent dir, there is a subdir called species-all. In this there is a
    file called taxonomy_lookup.tsv, which is a 2-column tsv. The first column is
    the taxID and the second column is the fasta header. There is also a second 
    file (taxonomy.tsv), which is a 4-column tsv. The first column is the taxID,
    the second column is the parentID, the third the rank and the fourth the name. 
    '''

    # fetch the taxonomy map file for the rtg_fasta headers
    rtg_fasta = PARAMS['data_rtg_fasta']
    rtg_dir = os.path.dirname(rtg_fasta)
    gi2tax_map = os.path.join(rtg_dir, 'species-all', 'taxonomy_lookup.tsv')
               
    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("rm -f {outfile};"
                     " cat {infile} |"
                     " cut -f2 |"
                     # stop at first match
                     " while read X; do grep -m 1 $X {gi2tax_map}"
                     " >> {outfile};"
                     " done")
        cluster_options = '-l walltime=00:02:00'
        P.run()


@transform(getGi2TaxId, suffix('_gi2taxid.tsv'), '_taxid2strain.tsv')
def getTaxIdStrains(infile, outfile):
    '''Fetch strain information based on tax id.'''

    # fetch the taxonomy file for the rtg_fasta headers
    rtg_fasta = PARAMS['data_rtg_fasta']
    rtg_dir = os.path.dirname(rtg_fasta)
    tax2strain_map = os.path.join(rtg_dir, 'species-all', 'taxonomy.tsv')

    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("rm -f {outfile};"
                     " cat {infile} |"
                     " cut -f1 |"
                     # stop at first match
                     " while read X; do grep -m 1 $X {tax2strain_map}"
                     " >> {outfile};"
                     " done")
        cluster_options = '-l walltime=00:02:00'
        P.run()
    
#@jobs_limit(1)
@transform(getTaxIdStrains, suffix('_taxid2strain.tsv'), '_strain2species.tsv')
def getStrainSpecies(infile, outfile):
    '''Fetch species information based on strain id'''

    # fetch the taxonomy file for the rtg_fasta headers
    rtg_fasta = PARAMS['data_rtg_fasta']
    rtg_dir = os.path.dirname(rtg_fasta)
    tax2strain_map = os.path.join(rtg_dir, 'species-all', 'taxonomy.tsv')

    if len(IOTools.openFile(infile).readlines()) == 0:
        open(outfile, 'w').close()
    else:
        statement = ("rm -f {outfile};"
                     " sleep 10;"
                     " cat {infile} |"
                     " cut -f2 |"
                     # stop at first match
                     " while read X;"
                     "  do grep -m 1 \"^$X\\s\" {tax2strain_map} |"
                     " grep species"
                     " >> {outfile};"
                     " sleep 5;"
                     " done")
        cluster_options = '-l walltime=00:40:00'
        to_cluster = False
        P.run()


# @jobs_limit(1)
# @follows(getStrainSpecies, getTaxIdStrains, getGi2TaxId, searchOTUsAgainstRTG)
# @transform(makeGenusSpecificOTUs,
#            regex('amp(.*)_species_otus.dir/amp(.*)_(.+).dir/distance_(.+)/(.+)_otus.fasta'),
#            add_inputs(r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches.b6',
#                       r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_gi2taxid.tsv',
#                       r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_strain2species.tsv',
#                       r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_taxid2strain.tsv'),
#            r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_classification.tsv')
@jobs_limit(1)
@follows(getStrainSpecies, getTaxIdStrains, getGi2TaxId, searchOTUsAgainstRTG)
@transform(makeGenusSpecificOTUs,
           regex('amp(.*)_species_otus.dir/amp(.*)_(bacteroides).dir/distance_(.+)/(.+)_otus.fasta'),
           add_inputs(r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches.b6',
                      r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_gi2taxid.tsv',
                      r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_strain2species.tsv',
                      r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_matches_taxid2strain.tsv'),
           r'amp\1_rtg_comparison.dir/amp\2_\3.dir/distance_\4/\5_otu_classification.tsv')
def fetchOTUClassificationSummary(infiles, outfile):
    '''Fetch OTU classification, based on RTG database'''
    otu_fasta, otu_b6, gi2taxid, strain2species, taxid2strain = infiles

    # fetch the OTU IDs from the fasta file
    otu_list = []
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(otu_fasta)):
        otu_list.append(fasta.title)

    ## Fetch dictionaries
    # from the b6 file... {(OTU, genbank_id):  identity,}
    b6_dict = {}
    b6_otus = set()
    for line in IOTools.openFile(otu_b6):
        line = line.split()
        b6_dict[(line[0], line[1])] = line[2]
        b6_otus.add(line[0])
    # from the gi2taxID file {genbank_id: taxid,}
    gi2taxid_dict = {}
    for line in IOTools.openFile(gi2taxid):
        line = line.split()
        gi2taxid_dict[line[1]] = line[0]
    # from the taxid2strain file {tax_id: [species_id, class_level, strain_name]}
    tax2strain_dict = collections.defaultdict(list)
    for line in IOTools.openFile(taxid2strain):
        line = line.strip().split('\t')
        tax2strain_dict[line[0]] = [line[1], line[2], line[3]]
    # from the strain2species file {species_id, species_name}
    strain2species_dict = {}
    for line in IOTools.openFile(strain2species):
        line = line.strip().split('\t')
        strain2species_dict[line[0]] = line[3]

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('\t'.join(['OTU', 'GenBankID', 'Identity', 'Species', 'Strain']) + '\n')
        for otu in otu_list:
            # need to catch those OTUs that didn't match anything in RTG
            if otu not in b6_otus:
                   outf.write('\t'.join([otu, 'None', 'None', 'None', 'None']) + '\n')
            for otu_id, genbank_id in b6_dict.keys():
                if otu == otu_id:
                    identity = b6_dict[(otu_id, genbank_id)]
                    tax_id = gi2taxid_dict[genbank_id]
                    species_id, class_level, strain_name = tax2strain_dict[tax_id]
                    if class_level == 'species':
                        species_name = strain_name
                        strain_name = 'unclassified'
                    else:
                        try:
                            species_name = strain2species_dict[species_id]
                        except KeyError:
                            species_name = 'unclassified'
                    outf.write('\t'.join([otu_id, genbank_id, identity, species_name, strain_name]) + '\n')
                else:
                    continue
        
###############################################################################
### Subtask: Collate mWGS and amplicon data
###############################################################################
@subdivide(fetchOTUClassificationSummary,
           regex('amp(.*)_rtg_comparison.dir/(.+).dir/distance_(.+)/(.+)_otu_classification.tsv'),
           add_inputs(r'amp\1_species_otus.dir/\2.dir/distance_\3/\2_otu_countTable.txt'),
           [r'amp\1_rtg_comparison.dir/\2.dir/distance_\3/\4_otu_strain_classification.tsv',
            r'amp\1_rtg_comparison.dir/\2.dir/distance_\3/\4_otu_species_classification.tsv'])
def collapseOTUClassificationSummary(infiles, outfiles):
    '''Collapse species-level OTU classifications by species or strain'''
    in_class, in_count = infiles
    out_strain, out_species = outfiles

    # some files may be empty
    if len(IOTools.openFile(in_class).readlines()) == 0 or \
       len(IOTools.openFile(in_count).readlines()) == 0:
        IOTools.openFile(out_strain, 'w').close()
        IOTools.openFile(out_species, 'w').close()
    else:
        df_class = pd.read_table(in_class)

        # drop duplicates to get a strain-level table
        df_strain = df_class.drop('GenBankID', axis=1)
        df_strain = df_strain.drop_duplicates()
        # drop duplicates to get a species-level table
        df_species = df_strain.drop(['Strain', 'Identity'], axis=1)
        df_species = df_species.drop_duplicates()
        df_species.set_index('OTU', inplace=True)
        df_strain = df_strain.drop(['Identity', 'Species'], axis=1)
        df_strain.set_index('OTU', inplace=True)
        
        # read in count table and replace headers
        df_count = pd.read_table(in_count, index_col=0)
        x = df_count.columns.tolist()
        for i, j in enumerate(x):
            if j == 'J03630' or j == 'J03630ds':
                x[i] = 'Scott'
            elif j == 'J03631' or j == 'J03631ds':
                x[i] = 'IronHorse'
            elif j == 'J03632' or j == 'J03632ds':
                x[i] = 'Commencal'
            elif j == 'J03644' or j == 'J03644ds':
                x[i] = 'Breezer'
            else:
                raise Exception('Unrecognized Sample ID: %s' % x[i])
        df_count.columns = x
    
        # merge tables
        df_species = df_species.join(df_count, how='outer')
        df_strain = df_strain.join(df_count, how='outer')
        df_species.to_csv(out_species, sep='\t', index_label='OTU')
        df_strain.to_csv(out_strain, sep='\t', index_label='OTU')
           


@jobs_limit(1) # the subdirectories are created during this task
@follows(loadmWGSStrainData)
@collate(collapseOTUClassificationSummary,
         regex('amp.*_rtg_comparison.dir/amp(.*?)(ds_|_)(.+).dir/'
               'distance_(.+)/(.+)_otu_(species|strain)_classification.tsv'), # or strain
         r'rtg_mwgs_data.dir/\3.dir/distance_\4/\3\2\6_abundance.tsv')
def collateOTUClassificationSummary(infiles, outfile):
    '''Combine V1V9, V1V3, and OTU abundance counts'''

    sub_dir = os.path.dirname(outfile)
    parent_dir = os.path.dirname(sub_dir)
    if not os.path.exists(parent_dir):
        os.mkdir(parent_dir)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)

    assert len(infiles) == 2, 'Expected two infiles: %s' % ' '.join(infiles)
    for infile in infiles:
        if os.path.basename(infile).startswith('amp_') or \
           os.path.basename(infile).startswith('ampds_'):
            in_v1v9 = infile
        elif os.path.basename(infile).startswith('ampV1V3_') or \
             os.path.basename(infile).startswith('ampV1V3ds_'):
            in_v1v3 = infile
        else:
            raise Exception('Unrecognised infile %s'  % infile)

    # fetch a dataframe containing the mWGS data from sqlite db
    species = os.path.basename(in_v1v9).split('_')[1]
    species = species[0].upper() + species[1:]
    if in_v1v9.endswith('strain_classification.tsv'):
        table = species + '_strain_abundance'
        tax = 'Strain'
    else:
        table = species + '_species_abundance'
        tax = 'Species'
    statement = ("SELECT * FROM {table}".format(**locals()))
    df_mwgs = pd.read_sql_query(sql=statement, con=connect(), index_col='taxonomy_name')
    x = df_mwgs.columns.tolist()
    x = [i.split('_')[0] + '_MGWS' for i in x]
    df_mwgs.columns = x
    
    # fetch dataframe containing the v1v3 data
    try:
        df_v1v3 = pd.read_table(in_v1v3, index_col=0)
        df_v1v3.index.name = 'V1V3'
        df_v1v3.reset_index(inplace=True)
        df_v1v3.set_index(tax, inplace=True)
        df_v1v3.index = ['None_v1v3' if x == 'None' else x for x in df_v1v3.index.tolist()]
        x = df_v1v3.columns.tolist()
        for i, j in enumerate(x):
            if j != 'V1V3':
                x[i] = j + '_V1V3'
        df_v1v3.columns = x
    except ValueError:
        df_v1v3 = pd.DataFrame(index=[tax,], columns=['V1V3',])
        
    # fetch a dataframe containing the v1v9 data
    try:
        df_v1v9 = pd.read_table(in_v1v9, index_col=0)
        df_v1v9.index.name = 'V1V9'
        df_v1v9.reset_index(inplace=True)
        df_v1v9.set_index(tax, inplace=True)
        df_v1v9.drop_duplicates(inplace=True)
        df_v1v9.index = ['None_v1v9' if x == 'None' else x for x in df_v1v9.index.tolist()]
        x = df_v1v9.columns.tolist()
        for i, j in enumerate(x):
            if j != 'V1V9':
                x[i] = j + '_V1V9'
        df_v1v9.columns = x
    except ValueError:
        df_v1v9 = pd.DataFrame(index=[tax,], columns=['V1V9',])

        
    # merge amplicon datasets and mwgs datasets
    df_amplicon = df_v1v9.join(df_v1v3, lsuffix='_v1v9' , rsuffix='_v1v3', how='outer')
    df = df_mwgs.join(df_amplicon, how='outer')
    
    # tidy up table
    df['V1V9'] = df['V1V9'].replace(np.nan, 'Missing')
    df['V1V3'] = df['V1V3'].replace(np.nan, 'Missing')
    df = df.replace(np.nan, 0)
    df.index = ['Missing' if x.startswith('None_') else x for x in df.index.tolist()]
    df.index.name = 'MWGS'
    
    df.to_csv(outfile, sep='\t')


@transform(collateOTUClassificationSummary,
           regex('rtg_mwgs_data.dir/bacteroides.dir/distance_1.0/(.+).tsv'),
           r'\1.load')
def loadOTUClassificationSummary(infile, outfile):
    ''' '''
    table_name = P.snip(infile, '.tsv', strip_path=True)
    statement = ("cat {infile} |"
                 " python {scriptsdir}/csv2db.py"
                 " --table={table_name}"
                 " --allow-empty-file"
                 " > {outfile}")
    to_cluster = False
    P.run()


@jobs_limit(1)
@transform(collateOTUClassificationSummary, suffix('.tsv'), '.png')
def plotOTUClassificationSuccess(infile, outfile):
    '''Plot a stacked bar chart to show how well OTUs correspond to 
    RTG mWGS data'''

    R('''rm(list=ls())''')
    def _summarizeColumn(df, column, genus):
        '''Summarize how well OTUs match entries in RTG output'''
        out_dict = {'02_Matches_multiple_species_in_RTG': 0,
                    '04_Not_present_in_RTG': 0,
                    '03_Matches_single_species_in_RTG': 0,
                    '01_Matches_wrong_genus_in_RTG': 0}
        
        series = df[column].value_counts()
        #print(series)
        for idx in series.index:
            if idx == 'Missing':
                continue
            elif series.loc[idx] > 1:
                #print idx
                #print series.loc[idx]
                matches_wrong_genus = False
                for i in df[df[column]==idx].index:
                    if not i.startswith(genus):
                        out_dict['01_Matches_wrong_genus_in_RTG'] += 1
                        #print 'Matches genus other than %s' % genus
                        matches_wrong_genus = True
                        break
                if not matches_wrong_genus:
                    out_dict['02_Matches_multiple_species_in_RTG'] += 1
                    #print 'Matches_multiple_species_in_RTG'
            else:
                #print idx
                assert len(df[df[column] == idx].index) == 1, \
                    'OTU should match only one taxon!'
                taxon = df[df[column]==idx].index[0]
                if taxon == 'Missing':
                    out_dict['04_Not_present_in_RTG'] += 1
                    #print 'Not_present_in_RTG'
                elif not taxon.startswith(genus):
                    out_dict['01_Matches_wrong_genus_in_RTG'] += 1
                        #print 'Matches genus other than %s' % genus
                else:
                    out_dict['03_Matches_single_species_in_RTG'] += 1
                    #print 'Matches_single_species_in_RTG'
                
        return out_dict

    # fetch the genus name
    genus = os.path.basename(infile).split('_')[0]
    genus = genus[0].upper() + genus[1:]
    # HACK... this will fail if the genus ends in ds... but I can't fix the names.
    if genus.endswith('ds'):
        genus = genus[:-2]

    # fetch a dataframe of the relationship between OTUs and RTG output
    df = pd.read_table(infile, index_col=0)
    df = df[['V1V9', 'V1V3']]

    # strip [] from the taxon ids in the index
    df.index = [re.sub('\[|\]', '', x) for x in df.index.tolist()]

    # create summary dataframe
    dt = {}
    for column in ['V1V9', 'V1V3']:
        dt[column] = _summarizeColumn(df, column, genus)

    df = pd.DataFrame(dt)
    # if this is strain rather than species... fix the rownames
    if re.search('_strain_', os.path.basename(infile)):
        df.index = [re.sub('_species_', '_strain_', x) for x in df.index.tolist()]

    # save the dataframe as a summary table
    outf = P.snip(outfile, '.png') + '_table.tsv'
    df.to_csv(outf, sep='\t', float_format='{:f}'.format, encoding='utf-8')

    df.reset_index(inplace=True)
    df = pd.melt(df, id_vars='index')

    # plot
    R.assign('df', df)
    R('''
    require("ggplot2")
    cols = c("darkred", "firebrick2", "lightblue3", "lightblue1")
    pl <- ggplot(df, aes(x=variable, y=value, fill=index)) + geom_bar(stat='identity')
    pl <- pl + scale_fill_manual(values=cols)
    pl <- pl + theme_bw() + theme(panel.grid = element_blank())
    pl <- pl + xlab('Amplicon Region') + ylab('Number of OTUs')
    pl <- pl + ggtitle("{genus}")
    png("{outfile}")
    plot(pl)
    dev.off()
    '''.format(**locals()))


@subdivide(plotOTUClassificationSuccess,
           regex('(.+).png'),
           r'\1_summary_*_table.tsv')
def summarizeOTUAssignment(infile, outfiles):
    '''Produce summary tables of the % otus correctly/incorrectly
    matched to RTG database'''
    infile = P.snip(infile, '.png') + '_table.tsv'
    
    # specify rows that fall into different categories
    d = {}
    d['error'] = ['01_Matches_wrong_genus_in_RTG', '02_Matches_multiple_species_in_RTG',]
    d['novel'] = ['04_Not_present_in_RTG',]
    d['correct'] = ['03_Matches_single_species_in_RTG',]
    if re.search('_strain_', os.path.basename(infile)):
        d['error'] = [re.sub('_species_', '_strain_', x) for x in d['error']]
        d['correct'] = [re.sub('_species_', '_strain_', x) for x in d['correct']]

    df = pd.read_table(infile, index_col=0)

    for cat in ['error', 'novel', 'correct']:
        outfile = P.snip(infile, '_table.tsv') + '_summary_' + cat + '_table.tsv'
        v1v3_t = float(sum(df['V1V3']))
        v1v3_s = float(sum(df.loc[d[cat], 'V1V3']))
        if v1v3_t == 0:
            v1v3 = 0.0
        else:
            v1v3 = (v1v3_s/v1v3_t)*100.0

        v1v9_t = float(sum(df['V1V9']))
        v1v9_s = float(sum(df.loc[d[cat], 'V1V9']))
        if v1v9_t == 0:
            v1v9 = 0
        else:
            v1v9 = (v1v9_s/v1v9_t)*100.0

        outf = IOTools.openFile(outfile, 'w')
        outf.write('\t'.join(map(str, [v1v3, v1v9, v1v3_t, v1v9_t])) + '\n')
#        outf.write(str(v1v3) + '\t' + str(v1v9) + '\n')
        outf.close()
    

@collate(summarizeOTUAssignment,
         regex('(.+)/.+.dir/distance(.+)/.+?(_strain_|ds_strain_|_species_|ds_species_)abundance_summary_(correct|error|novel)_table.tsv'),
         r'\1/distance\2\3\4_table.tsv')
def collateOTUAssignmentSummary(infiles, outfile):
    '''Combine summaries of the correct, novel, erroneous OTUs'''
    infiles = " ".join(infiles)
    statement = ("python {scriptsdir}/combine_tables.py"
                 " --cat=Genus"
                 #" --header-names=Genus,V1V3,V1V9"
                 " --no-titles"
                 " --regex-filename='.+/(.+?)_.+tsv'"
                 " --log={outfile}.log"
                 " {infiles} > {outfile}")
    to_cluster = False
    P.run()

        
@jobs_limit(1)
@collate(plotOTUClassificationSuccess,
         regex('(.+)/.+.dir/distance(.+)/'
               '.+?(_strain_|ds_strain_|_species_|ds_species_)abundance.png'),
         r'\1/distance\2\3plots.pdf')
def collateOTUClassificationSuccessPlots(infiles, outfile):
    '''Merge plots into a single pdf'''
    infiles = " ".join(infiles)
    statement = "convert {infiles} {outfile}"
    to_cluster=False
    P.run()


###############################################################################
# Metatask: MAP READS BACK ONTO OTU SEQUENCE
###############################################################################
@jobs_limit(8)
@follows(poolAMPReads_usearch, poolDownsampledAMPReads_usearch) #, poolGenusReads)
@subdivide(assignAMPReadsToOTUs, #assignGeneraReadsToOTUs],
           regex('(.+)/(.+)/(.+?)_(.+)_otu_assignment.txt'),
           add_inputs(r'\1/\3_pooled.fasta.gz'),
           r'\1/\2/otu_mapping/\3_\4_*.fasta.gz')
def splitReadsMappedToOTUs(infiles, outfiles):
    '''Separate reads mapped to different OTUs by parsing pooled fasta'''
    otu_file, fasta_file = infiles
    
    # generate outfile prefix
    out_stub = os.path.join(os.path.dirname(otu_file), 'otu_mapping')
    if not os.path.exists(out_stub):
        os.mkdir(out_stub)
    out_prefix = P.snip(otu_file, 'otu_assignment.txt', strip_path=True)
    out_prefix = os.path.join(out_stub, out_prefix)

    # create a dictionary where keys are assigned reads and value is otu
    otus = set()
    read_dict = {}
    for line in IOTools.openFile(otu_file):
        if line.strip().endswith('*'):
            continue
        else:
            line = line.split()
            otu = line[-1]
            read = line[-2]
            read_dict[read] = otu
            otus.add(otu)

    # generate a dictionary of outfile handles
    outfile_dict = {}
    for otu in otus:
        outf = out_prefix + otu + '.fasta.gz'
        outfile_dict[otu] = IOTools.openFile(outf, 'w')


    # iterate over the fasta files and write relevant reads to their respective
    # outfiles
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        fasta.title = re.sub('\s', '_', fasta.title)
        if fasta.title in read_dict.keys():
            outf = outfile_dict[read_dict[fasta.title]]
            outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
        else:
            continue

    # close the outfiles
    for outf in outfile_dict.values():
        outf.close()

        
@subdivide(makeAMPOTUs,
           regex('(.+)/(.+)/(.+)_otus.fasta'),
           r'\1/\2/otu_mapping/\3_\2_*_template.fasta')
def makeOTUTemplates(infile, outfile):
    '''Split each otu sequence into a separate fasta file'''

    # generate an outfile prefix
    out_dist = os.path.basename(os.path.dirname(infile))
    out_stub = os.path.join(os.path.dirname(infile), 'otu_mapping')
    assert os.path.exists(out_stub), 'otu_mapping directory is missing'
    out_prefix = P.snip(infile, 'otus.fasta', strip_path=True)
    out_prefix = out_prefix + out_dist + '_'
    out_prefix = os.path.join(out_stub, out_prefix)
    
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        outf = out_prefix + fasta.title + '_template.fasta'
        with IOTools.openFile(outf, 'w') as o:
            o.write('>' + fasta.title + '\n' + fasta.sequence + '\n')


@follows(poolGenusReads)
@subdivide(assignGeneraReadsToOTUs,
           regex('amp(.*)_species_otus.dir/(.+)/(.+)_otu_assignment.txt'),
           add_inputs(r'amp\1_classified.dir/\3.fasta'),
           r'amp\1_species_otus.dir/\2/otu_mapping/\3_*.fasta.gz')
def splitReadsMappedToSpeciesSpecificOTUs(infiles, outfiles):
    otu_file, fasta_file = infiles
    
    # generate outfile prefix
    out_stub = os.path.join(os.path.dirname(otu_file), 'otu_mapping')
    if not os.path.exists(out_stub):
        os.mkdir(out_stub)
    out_prefix = P.snip(otu_file, 'otu_assignment.txt', strip_path=True)
    out_prefix = os.path.join(out_stub, out_prefix)

    # create a dictionary where keys are assigned reads and value is otu
    otus = set()
    read_dict = {}
    for line in IOTools.openFile(otu_file):
        if line.strip().endswith('*'):
            continue
        else:
            line = line.split()
            otu = line[-1]
            read = line[-2]
            read_dict[read] = otu
            otus.add(otu)

    # generate a dictionary of outfile handles
    outfile_dict = {}
    for otu in otus:
        outf = out_prefix + otu + '.fasta.gz'
        outfile_dict[otu] = IOTools.openFile(outf, 'w')


    # iterate over the fasta files and write relevant reads to their respective
    # outfiles
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        fasta.title = re.sub('\s', '_', fasta.title)        
        if fasta.title in read_dict.keys():
            outf = outfile_dict[read_dict[fasta.title]]
            outf.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
        else:
            continue

    # close the outfiles
    for outf in outfile_dict.values():
        outf.close()


@subdivide(makeGenusSpecificOTUs,
           regex('(.+)/(.+)_otus.fasta'),
           r'\1/otu_mapping/\2_*_template.fasta')
def makeSpeciesOTUTemplates(infile, outfiles):
    '''Split each otu sequence into a separate fasta file'''

    # generate an outfile prefix
    out_stub = os.path.join(os.path.dirname(infile), 'otu_mapping')
    assert os.path.exists(out_stub), 'otu_mapping directory is missing'
    out_prefix = P.snip(infile, 'otus.fasta', strip_path=True)
    out_prefix = os.path.join(out_stub, out_prefix)
    
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        outf = out_prefix + fasta.title + '_template.fasta'
        with IOTools.openFile(outf, 'w') as o:
            o.write('>' + fasta.title + '\n' + fasta.sequence + '\n')
    

###############################################################################
## Task: Run Cross Match
###############################################################################
@follows(makeOTUTemplates, makeSpeciesOTUTemplates)
@transform([splitReadsMappedToOTUs, splitReadsMappedToSpeciesSpecificOTUs],
           regex('(.+)/(.+).fasta.gz'),
           add_inputs(r'\1/\2_template.fasta'),
           r'\1/\2_cross_match.txt')
def mapOTUReadsToTemplate(infiles, outfile):
    '''Run CrossMatch to map reads onto their OTU template sequence'''
    query_fasta, template_fasta = infiles

    # create a temporay unzipped query file
    out_dir = os.path.dirname(outfile)
    tmp_inf = P.getTempFilename(out_dir)

    # The cross match parameters are set as:
    # mask_level = 0 (only report the best alignment for each query)
    # min_score = 430/1300 (approx. 90% of the length of the template reads)
    if re.search('V1V3', os.path.basename(query_fasta)):
        min_score = str(430)
    else:
        min_score = str(1300)

    statement = ("zcat {query_fasta} > {tmp_inf};"
                 " cross_match"
                 "  -discrep_lists"
                 "  -tags"
                 "  -minscore {min_score}"
                 "  -masklevel 0"
                 " {tmp_inf}"
                 " {template_fasta}"
                 " > {outfile}")
    P.run()

    shutil.move(tmp_inf + '.log', outfile + '.log')
    os.unlink(tmp_inf)


@transform(mapOTUReadsToTemplate,
           regex('(.+)/(.+)_cross_match.txt'),
           add_inputs(r'\1/\2_template.fasta'),
           r'\1/\2_cross_match_summary.txt')
def summarizeReadToOTUMapping(infiles, outfile):
    '''Create a summary dataframe out of the cross match results'''
    infile, ref_db = infiles

    # WARNING!
    if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
        L.warn('There is no sanity check for duplicate matches in'
               ' default cross_match run') 

    # get the maximum length of the sequences in the reference database
    # (informs output table size)
    max_len = 0
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(ref_db)):
        if len(fasta.sequence) > max_len:
            max_len = len(fasta.sequence)

    # parse cross match output
    out_dict = P16S.parseCrossMatchDiscrepList(infile, max_len, outfile)

    # pickle and dump out_dict
    out_pickle = P.snip(outfile, '.txt') + '.p'
    pickle.dump(out_dict, open(out_pickle, 'wb'))


@jobs_limit(1)
@subdivide(summarizeReadToOTUMapping,
           regex('(.+)/(.+)_summary.txt'),
#           regex('(amp_genera_otus.dir/distance_1.0/otu_mapping)/(.+)_summary.txt'),
           r'\1/\2_errorProfile.pdf')
def plotOTUErrorHistograms(infile, outfile):
    '''Plot error histograms
    '''

    pickle_file = P.snip(infile, '.txt') + '.p'
    outdir = os.path.dirname(pickle_file)

    # read in the table containing all mismatches (this is large!)
    df = pd.DataFrame(pickle.load(open(pickle_file, 'rb')))

    # fetch the OTU_ID
    otu_id = P.snip(infile, '_cross_match_summary.txt', strip_path=True)
    otu_id = '_'.join(otu_id.split('_')[-2:])

    # select mismatch dataframe
    df = df[otu_id]
    n_seq = str(len(df.columns))
    # pandas counts empty strings as True
    df = df.applymap(lambda x: np.nan if x == '' else x)

    # plot substitutions, insertions, deletions
    R('''rm(list=ls())''')
    R('''pdf('%s', height=3.5)''' % outfile)
    R('''par(mfrow=c(1,3))''')
    for error in ('S', 'I', 'D'):
        if error == 'D':
            df_sub1 = df.applymap(lambda x: np.nan
                                  if x in ['S', 'I'] else x)
        elif error == 'S':
            df_sub1 = df.applymap(lambda x: np.nan
                                  if x in ['D', 'I'] else x)
        elif error == 'I':
            df_sub1 = df.applymap(lambda x: np.nan
                                  if x in ['S', 'D'] else x)
        else:
            df_sub1 = df.copy()
        df_sub1 = df_sub1.count(axis=1)
        df_sub1 = df_sub1.to_frame()
        df_sub1.reset_index(inplace=True)
        df_sub1.columns = ['Index', 'Count']

        # set to robjects
        x_axis = rpy2.robjects.vectors.IntVector(df_sub1['Index'].tolist())
        y_axis = rpy2.robjects.vectors.IntVector(df_sub1['Count'].tolist())

        R.assign('xa', x_axis)
        R.assign('ya', y_axis)

        R('''
          plot(xa,
               ya/{n_seq}*100,
               type='l',
               main=paste('{otu_id}', '\n(', {n_seq}, ' sequences aligned)', sep=''),
               xlab='Location',
               ylab='% Error',
               ylim=c(0,100))
        '''.format(**locals()))

    R('''dev.off()''')
    R('''rm(list=ls())''')
    L.info('Plots completed for otu {}, written to {}'.format(otu_id, outfile))

    
@collate(summarizeReadToOTUMapping,
           regex('(.+)/(.+)_summary.txt'),
           r'\1/otu_error_summary.tsv')
def collateOTUErrors(infiles, outfile):
    '''Create a dataframe summarizing the position of errors in reads
    aligned to each OTU'''

    first = True
    for otu_file in infiles:
        pickle_file = P.snip(otu_file, '.txt') + '.p'
        df = pd.DataFrame(pickle.load(open(pickle_file, 'rb')))
        if first:
            df_i = pd.DataFrame(index=range(0, len(df.index)))
            df_s = pd.DataFrame(index=range(0, len(df.index)))
            df_d = pd.DataFrame(index=range(0, len(df.index)))
            df_a = pd.DataFrame(index=range(0, len(df.index)))
            first = False
            
        # there is only one otu per dataframe
        otu_id = list(set(df.columns.get_level_values(0)))[0]
        df = df[otu_id]
        df = df.applymap(lambda x: np.nan if x == '' else x)

        for error in ('S', 'I', 'D'):
            if error == 'D':
                df_sub1 = df.applymap(lambda x: np.nan
                                      if x in ['S', 'I'] else x)
                df_sub1 = df_sub1.count(axis=1)
                df_sub1 = df_sub1.to_frame()
                df_sub1.columns = [otu_id,]
                df_d = df_d.join(df_sub1)
            elif error == 'S':
                df_sub1 = df.applymap(lambda x: np.nan
                                      if x in ['D', 'I'] else x)
                df_sub1 = df_sub1.count(axis=1)
                df_sub1 = df_sub1.to_frame()
                df_sub1.columns = [otu_id,]
                df_s = df_s.join(df_sub1)
            elif error == 'I':
                df_sub1 = df.applymap(lambda x: np.nan
                                      if x in ['S', 'D'] else x)
                df_sub1 = df_sub1.count(axis=1)
                df_sub1 = df_sub1.to_frame()
                df_sub1.columns = [otu_id,]
                df_i = df_i.join(df_sub1)
            else:
                df_sub1 = df.copy()
                df_sub1 = df_sub1.count(axis=1)
                df_sub1 = df_sub1.to_frame()
                df_sub1.columns = [otu_id,]
                df_a = df_a.join(df_sub1)
                
    outf_sub = P.snip(outfile, '.tsv') + '_substitutions.tsv'
    df_s.to_csv(outf_sub, sep='\t')
    outf_ins = P.snip(outfile, '.tsv') + '_insertions.tsv'
    df_i.to_csv(outf_ins, sep='\t')
    outf_del = P.snip(outfile, '.tsv') + '_deletions.tsv'
    df_d.to_csv(outf_del, sep='\t')
    df_a.to_csv(outfile, sep='\t')

            
@jobs_limit(1)
@transform(plotOTUErrorHistograms, suffix('.pdf'), '.png')
def convertOTUErrorHistogramsToPng(infile, outfile):
    '''For the purpose of collating'''
    statement = "convert {infile} {outfile}"
    to_cluster = False
    P.run()
    
@jobs_limit(1)
@collate(convertOTUErrorHistogramsToPng,
         regex('(.+)/(.+)_OTU_.*_cross_match_errorProfile.png'),
         r'\1/\2_errorProfile_summary.pdf')
def collateEntropyPlots(infiles, outfile):
    '''Merge entropy plots into a single pdf'''
    infiles = " ".join(infiles)
    statement = "convert {infiles} {outfile}"
    to_cluster = False
    P.run()


# @transform(fetchOTUStrainIDsFromNCBI,)
#            regex('amp(.*)_rtg_comparison.dir/(.+)/(.+)/(.+)_otus.b8'),
#            r'amp\1_rtg_comparison.dir/\2/\3/\3_\2_otu_classification.tsv')
# def classifyRTGMatchedOTUs(infile, outfile):
#     '''Combine OTU alignment table with the Genbank taxonomic IDs'''
#     df = pd.read_table(infile, header=None)
#     df = df[[0, 1, 2, 12]]
#     df.columns = ['OTU_ID', 'genBank_ID', 'Identity', 'Strain_ID']
#     df.to_csv(outfile, index=False, sep='\t')
    

# @transform(classifyRTGMatchedOTUs,
#            regex('.+/(.+).tsv'),
#            r'\1.load')
# def loadRTGMatchedOTUs(infile, outfile):
#     table_name = P.snip(infile, '.tsv', strip_path=True)
#     to_cluster = False
#     statement = ("cat {infile} |"
#                  " python {scriptsdir}/csv2db.py"
#                  " --table={table_name}"
#                  " > {outfile}")
#     P.run()



# ###############################################################################
# ### Subtask: Classify the bacteroides strains against mothur RDP
# ###############################################################################
# @subdivide(makeBacteroidesOTUs,
#            regex('(.+)/(.+)/(.+).fasta'),
#            r'\1/\2/*/\3.fasta')
# def createBacteroidesRDPSubdirectories(infile, outfiles):
#     '''Create subdirectory for each of the RDP thresholds to be run'''
#     outfile = os.path.basename(infile)
#     out_parent_dir = os.path.dirname(infile)

#     for threshold in PARAMS['mothur_cutoff'].split(','):
#         directory = os.path.join(out_parent_dir, 'threshold_' + str(threshold))
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         outf = os.path.join(directory, outfile)
#         to_cluster = False
#         statement = "cp -f {infile} {outf}"
#         P.run()


# @subdivide(makeBacteroidesOTUs,
#            regex('(.+)/(.+)/(.+).fasta'),
#            add_inputs(PARAMS['in_vivo_rdp_tax']),
#            r'\1/\2/rdp_taxonomy.tax')
# def createBacteroidesRDPTaxonomyFiles(infiles, outfile):
#     '''Parse the RDP taxonomy file so that the sequence ID is in species'''
#     otu_file, tax_file = infiles

#     tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

#     statement = ("cat {tax_file} |"
#                  " python {tax_parser}"
#                  "  --log={outfile}.log"
#                  "  --replace-species"
#                  " > {outfile}")
#     to_cluster = False
#     P.run()
    

# @transform(createBacteroidesRDPSubdirectories,
#            regex('(.+)/(.+)/(.+)/(.+).fasta'),
#            add_inputs(PARAMS['in_vivo_rdp_fasta'],
#                       r'\1/\2/rdp_taxonomy.tax'),
#            r'\1/\2/\3/\4.rdp_taxonomy.wang.taxonomy')
# def runMothurRDP_bacteroides(infiles, outfile):
#     '''Run mothur's version of the RDP classifier from the commandline'''
#     input_fasta, ref_fasta, ref_taxonomy = infiles

#     # get cut off from directory name
#     cut_off = os.path.basename(os.path.dirname(input_fasta))
#     cut_off = cut_off.split('_')[1]

#     working_dir = os.path.dirname(input_fasta)

#     tmp_input_fasta = os.path.abspath(input_fasta)
#     tmp_ref_fasta = os.path.abspath(ref_fasta)
#     tmp_ref_tax = os.path.abspath(ref_taxonomy)

#     # set outfile
#     outfile = P.snip(input_fasta, '.fasta') + \
#         P.snip(os.path.basename(ref_taxonomy), '.tax') + \
#         '.wang.taxonomy'

#     # run mothur classify.seqs()
#     statement = ("cd {working_dir};"
#                  " mothur \"#classify.seqs("
#                  "fasta={tmp_input_fasta},"
#                  " template={tmp_ref_fasta},"
#                  " taxonomy={tmp_ref_tax},"
#                  " cutoff={cut_off},"
#                  " processors={mothur_processors},"
#                  " iters={mothur_iters})\";"
#                  " cd -")

#     cluster_options = '-l nodes=1:ppn=2' 
#     to_cluster = True

#     P.run()


# @transform(runMothurRDP_bacteroides,
#            suffix('s.rdp_taxonomy.wang.taxonomy'),
#            '_taxonomic_asignment.txt')
# def parseMothurRDP_bacteroides(infile, outfile):
#     '''Convert taxonomy file to table'''

#     tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')
#     statement = ("cat {infile} |"
#                  "python {tax_parser}"
#                  " --remove-confidence-thresholds"
#                  " --out-format=table"
#                  " --log={outfile}.log"
#                  " > {outfile}")
#     cluster_options = '-l walltime=00:02:00'
#     P.run()


# ###############################################################################
# ### Subtask: Collate mWGS and amplicon taxonomic identification tables
# ###############################################################################
# @follows(loadRTGOutput)
# @follows(loadRTGMatchedOTUs)
# @transform(classifyRTGMatchedOTUs,
#            suffix('_classification.tsv'),
#            add_inputs(loadRTGOutput),
#            '_mWGS_match.tsv')
# def compareMWGSAndAmpliconClassification(infiles, outfile):
#     '''Collate information on the strains detected by mWGS vs amplicon
#     sequencing.
#     '''
#     amplicon_table, mwgs_table = infiles
#     amplicon_table = P.snip(amplicon_table, '.tsv', strip_path=True)
#     mwgs_table = P.snip(mwgs_table, '.load', strip_path=True)
#     amplicon_table = re.sub('\.', '_', amplicon_table)
#     mwgs_table = re.sub('\.', '_', mwgs_table)
    
#     statement = ("SELECT DISTINCT OTU_ID,taxonomy_name,Strain_ID"
#                  "  FROM {mwgs_table}"
#                  "  LEFT OUTER JOIN {amplicon_table}"
#                  "   ON taxonomy_name = Strain_ID"
#                  " UNION ALL"
#                  " SELECT DISTINCT OTU_ID,taxonomy_name,Strain_ID"
#                  "  FROM {amplicon_table}"
#                  "  LEFT OUTER JOIN {mwgs_table}"
#                  "   ON taxonomy_name = Strain_ID".format(**locals()))
#     df = pd.read_sql(sql=statement, con=connect())

#     df.to_csv(outfile, index=False, sep='\t')


# @follows(mkdir('rtg_amp_comparison.dir'))
# @collate(compareMWGSAndAmpliconClassification,
#          regex('.+/amp.*_stool_(.+).tsv'),
#          r'rtg_amp_comparison.dir/\1.tsv')
# def compareV3vsV9mWGSMatch(infiles, outfile):
#     '''Join the tables for V1-V3 and V1-V9 AMP mWGS comparison'''
#     for infile in infiles:
#         if os.path.basename(infile).startswith('ampV1V3'):
#             v1v3 = infile
#         elif os.path.basename(infile).startswith('amp_stool'):
#             v1v9 = infile
#         else:
#             raise Exception('Unrecognised infile')

#     assert v1v3, 'missing V1-V3 file'
#     assert v1v9, 'missing V1-V9 file'

#     df_v3 = pd.read_table(v1v3)
#     df_v3.drop(['Strain_ID',], inplace=True, axis=1)
#     df_v3.set_index('taxonomy_name', inplace=True)
#     df_v3.columns = ['V1_V3',]
    
#     df_v9 = pd.read_table(v1v9)
#     df_v9.drop(['Strain_ID',], inplace=True, axis=1)
#     df_v9.set_index('taxonomy_name', inplace=True)
#     df_v9.columns = ['V1_V9',]

#     df_out = df_v3.join(df_v9)
#     df_out.reset_index(inplace=True)
#     df_out.drop_duplicates(inplace=True)
#     df_out.to_csv(outfile, sep='\t', index=False)


# @jobs_limit(1)
# @transform(compareV3vsV9mWGSMatch,
#            regex('(.+)/bacteroides_(.+)_otu_mWGS_match.tsv'),
#            add_inputs(r'ampV1V3_stool_otus.dir/\2/ampV1V3_stool_bacteroides_otus.fasta',
#                       r'amp_stool_otus.dir/\2/amp_stool_bacteroides_otus.fasta'),
#            r'\1/bacteroides_\2_otu_mWGS_summary.tsv')
# def createV3vsV9mWGSSummary(infiles, outfile):
#     '''Assumes match file has columns taxonomy_name, V1_V3, V1_V9
#     '''
#     match_file, v1v3_fasta, v1v9_fasta = infiles

#     # fetch the full list of OTUs output by USEARCH
#     def fetch_otus(fasta_file):
#         l = []
#         for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
#             l.append(fasta.title)
#         return l

#     v1v3_otus = fetch_otus(v1v3_fasta)
#     v1v9_otus = fetch_otus(v1v9_fasta)

#     def summarizeOTUs(df, column, all_otus):
#         # extract the OTUs from column
#         otu_list = [x for x in set(df[column]) if not pd.isnull(x)]

#         # set up summary dicts
#         amplicon_strain_dict = collections.defaultdict(list)
#         amplicon_otu_dict = {}

#         for otu in otu_list:
#             df_sub = df[df[column] == otu]
#             # record the number of OTUs that map to a particular strain
#             for i in df_sub.index:
#                 # OTUs that don't match mWGS strains are currently NaN
#                 if pd.isnull(i):
#                     amplicon_strain_dict['Missing'].append(otu)
#                 else:
#                     amplicon_strain_dict[i].append(otu)
#                     # record whether a an OTU matches one, or multiple strains
#                 if len(df_sub.index) > 1:
#                     amplicon_otu_dict[otu] = 'Multiple'
#                 elif pd.isnull(df_sub.index.tolist()[0]):
#                     amplicon_otu_dict[otu] = 'Error'
#                 else:
#                     amplicon_otu_dict[otu] = 'Single'
#         # check for strains have multiple matching OTUs
#         for strain, otus in amplicon_strain_dict.items():
#             if strain == 'Missing':
#                 continue
#             if len(otus) > 1:
#                 for otu in otus:
#                     amplicon_otu_dict[otu] = 'Partial'
#                 else:
#                     continue
#         # finally, if there are OTUs not detected by mWGS
#         for otu in all_otus:
#             if otu not in amplicon_otu_dict.keys():
#                 amplicon_otu_dict[otu] = 'Novel'
                
#         return amplicon_otu_dict
                                                                                                                
#     # summarize the otus for V1-V3 & V1-V9
#     df = pd.read_table(match_file)
#     df.set_index('taxonomy_name', inplace=True)
#     v1v3 = collections.Counter(summarizeOTUs(df, 'V1_V3', v1v3_otus).values())
#     df_v1v3 = pd.DataFrame(v1v3.items(), columns=['Status', 'V1_V3'])
#     df_v1v3.set_index('Status', inplace=True)

#     v1v9 = collections.Counter(summarizeOTUs(df, 'V1_V9', v1v9_otus).values())
#     df_v1v9 = pd.DataFrame(v1v9.items(), columns=['Status', 'V1_V9'])
#     df_v1v9.set_index('Status', inplace=True)

#     df_out = df_v1v3.join(df_v1v9)
#     df_out.to_csv(outfile, sep='\t')


###############################################################################
# Metatask: SEARCH FOR FULL LENGTH READS THAT ARE NOT IN GG DATABASE
###############################################################################
## Task: Parse the full greengenes database. 
###############################################################################
@follows(mkdir('mouse_novel_amplicons.dir'))
@follows(mkdir('amp_novel_amplicons.dir'))
@originate(['mouse_novel_amplicons.dir/greengenes_taxonomy.txt',
            'amp_novel_amplicons.dir/greengenes_taxonomy.txt'])
def parseFullGGTaxonomyFiles(outfile):
    '''Parse the taxonomy file for the full greengenes database for use
    with mothurs RDP classifier '''
    infile = PARAMS['in_vivo_full_taxonomy']
    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

    statement = ("zcat {infile} |"
                 " python {tax_parser}"
                 "  --log={outfile}.log"
                 "  --tax-split='; '" # this also removes prefixes!
                 "  --replace-species"
#                 "  --drop-empty"
                 "  -x ID"
                 "  --remove-whitespace"
                 "  --clean-greengenes"
                 " > {outfile}")
    to_cluster = False
    P.run()

        
@follows(parseFullGGTaxonomyFiles)
@originate(['mouse_novel_amplicons.dir/greengenes.fasta',
            'amp_novel_amplicons.dir/greengenes.fasta'])
def parseFullGGFastaFiles(outfile):
    '''Parse the greengenes reference fasta file to include only those
    sequences that are in the parsed taxonomy file'''
    infile = PARAMS['in_vivo_full_sequence']
    statement = ("zcat {infile} | sed 's/>/>ID/' > {outfile}")
    to_cluster = False
    P.run()
    

###############################################################################
### Subtask: Assign taxonomy to  OTUs/dereplicated sequences
###############################################################################
@follows(parseFullGGTaxonomyFiles)
@follows(parseFullGGFastaFiles)
@transform(makeAMPOTUs,
           regex('(amp|mouse)_genera_otus.dir/(.+)/(.+).fasta'),
           add_inputs(r'\1_novel_amplicons.dir/greengenes.fasta',
                      r'\1_novel_amplicons.dir/greengenes_taxonomy.txt'),
           r'\1_novel_amplicons.dir/\2/\3.greengenes_taxonomy.wang.taxonomy')
def runMothurRDP_otu_denovo_discovery(infiles, outfile):
    '''Classify AMP OTUs against greengenes'''
    input_fasta, ref_fasta, ref_taxonomy = infiles
    fasta = os.path.basename(input_fasta)
    
    ref_fasta = os.path.abspath(ref_fasta)
    ref_taxonomy = os.path.abspath(ref_taxonomy)
    input_fasta=os.path.abspath(input_fasta)
    in_vivo_cut_off = PARAMS['in_vivo_cut_off']
    in_vivo_mothur_processors = PARAMS['in_vivo_mothur_processors']
    
    # for safety keep parallel mothur runs in separate directories
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))
    working_dir = os.path.dirname(outfile)

    statement = ("cd {working_dir};"
                 " ln -s {ref_fasta} .;"
                 " ln -s {ref_taxonomy} .;"
                 " ln -s {input_fasta} .;"
                 " mothur \"#classify.seqs("
                 "fasta={fasta},"
                 " template=greengenes.fasta,"
                 " taxonomy=greengenes_taxonomy.txt,"
                 " cutoff={in_vivo_cut_off},"
                 " iters=75,"
                 " processors={in_vivo_mothur_processors})\""
                 " cd -")

    cluster_options = '-l walltime=20:00:00,mem=50Gb,nodes=1:ppn={}'.format(PARAMS['in_vivo_mothur_processors'])
    P.run()


@follows(parseFullGGTaxonomyFiles)
@follows(parseFullGGFastaFiles)
@follows(mkdir('amp_full_novel_amplicons.dir'))
@follows(mkdir('mouse_full_novel_amplicons.dir'))
@transform(generateAMPOTUSubdirectories,
           regex('(amp|mouse)_genera_otus.dir/(.+)/(.+).fasta'),
           add_inputs(r'\1_novel_amplicons.dir/greengenes.fasta',
                      r'\1_novel_amplicons.dir/greengenes_taxonomy.txt'),
           r'\1_full_novel_amplicons.dir/\2/\3.greengenes_taxonomy.wang.taxonomy')
def runMothurRDP_full_denovo_discovery(infiles, outfile):
    '''Classify AMP OTUs against greengenes'''
    input_fasta, ref_fasta, ref_taxonomy = infiles
    fasta = os.path.basename(input_fasta)
    
    ref_fasta = os.path.abspath(ref_fasta)
    ref_taxonomy = os.path.abspath(ref_taxonomy)
    input_fasta=os.path.abspath(input_fasta)
    in_vivo_cut_off = PARAMS['in_vivo_cut_off']
    in_vivo_mothur_processors = PARAMS['in_vivo_mothur_processors']
    
    # for safety keep parallel mothur runs in separate directories
    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))
    working_dir = os.path.dirname(outfile)

    statement = ("cd {working_dir};"
                 " ln -s {ref_fasta} .;"
                 " ln -s {ref_taxonomy} .;"
                 " ln -s {input_fasta} .;"
                 " mothur \"#classify.seqs("
                 "fasta={fasta},"
                 " template=greengenes.fasta,"
                 " taxonomy=greengenes_taxonomy.txt,"
                 " cutoff={in_vivo_cut_off},"
                 " iters=75,"
                 " processors={in_vivo_mothur_processors})\""
                 " cd -")

    cluster_options = '-l walltime=20:00:00,mem=50Gb,nodes=1:ppn={}'.format(PARAMS['in_vivo_mothur_processors'])
    P.run()


@follows(runMothurRDP_full_denovo_discovery,
         runMothurRDP_otu_denovo_discovery)
def run_rdp_denovo_discovery():
    pass


@transform([runMothurRDP_otu_denovo_discovery, runMothurRDP_full_denovo_discovery],
           regex('(.+)/(.+)/(.+).greengenes_taxonomy.wang.taxonomy'),
           r'\1/\2/\3_otu_taxonomy.txt')
def fetchReadTaxonomy_denovo_discovery(infile, outfile):
    '''Parse output of mothur taxonomic assignment'''
    tax_parser = os.path.join(os.path.dirname(__file__), 'tax2tax.py')

    statement = ("cat {infile} |"
                 " python {tax_parser}"
                 "  --log={outfile}.log"
                 "  --out-format=table |"
                 " grep unclassified"
                 " > {outfile}")
    cluster_options = '-l walltime=00:20:00'
    to_cluster = False
    P.run()

###############################################################################
### Subtask: Local alignment of OTUs/unique sequences to full gg database. 
###############################################################################
@follows(run_rdp_denovo_discovery)
@transform(makeAMPOTUs,
           regex('(amp|mouse)_genera_otus.dir/(.+)/(.+).fasta'),
           add_inputs(indexGreenGenes),
           r'\1_novel_amplicons.dir/\2/\3_gg_matches.b6')
def usearchGreenGenes_otu_denovo_discovery(infiles, outfile):
    amplicon_fasta, ref_db = infiles

    # fetch distance
    dist = float(os.path.dirname(outfile).split('_').pop())
    dist = (100-dist)/100.0

    # output sam file
    sam_out = P.snip(outfile, '.b6') + '.sam'
    
    statement = ("usearch -usearch_local {amplicon_fasta}"
                 " -db {ref_db}"
                 " -id {dist}"
                 " -maxaccepts 1" # terminate as soon as one positive hit made
                 " -maxrejects 0"
                 " -strand both"
                 " -output_no_hits"
                 " -query_cov {dist}"
                 " -blast6out {outfile}"
                 " -samout {sam_out}"
                 " -threads 16"
                 " &> {outfile}.log")
    cluster_options = '-l walltime=24:00:00,mem=60Gb,nodes=1:ppn=16'
    P.run()


@follows(run_rdp_denovo_discovery)
@transform(generateAMPOTUSubdirectories,
           regex('(amp|mouse)_genera_otus.dir/(.+)/(.+).fasta'),
           add_inputs(indexGreenGenes),
           r'\1_full_novel_amplicons.dir/\2/\3_gg_matches.b6')
def usearchGreenGenes_full_denovo_discovery(infiles, outfile):
    amplicon_fasta, ref_db = infiles

    # fetch distance
    dist = float(os.path.dirname(outfile).split('_').pop())
    dist = (100-dist)/100.0

    # output sam file
    sam_out = P.snip(outfile, '.b6') + '.sam'
    
    statement = ("usearch -usearch_local {amplicon_fasta}"
                 " -db {ref_db}"
                 " -id {dist}"
                 " -maxaccepts 1" # terminate as soon as one positive hit made
                 " -maxrejects 0"
                 " -strand both"
                 " -output_no_hits"
                 " -query_cov {dist}"
                 " -blast6out {outfile}"
                 " -samout {sam_out}"
                 " -threads 16"
                 " &> {outfile}.log")
    cluster_options = '-l walltime=24:00:00,mem=60Gb,nodes=1:ppn=16'
    P.run()


@follows(usearchGreenGenes_full_denovo_discovery,
         usearchGreenGenes_otu_denovo_discovery)
def run_usearch_denovo_discovery():
    pass

###############################################################################
### Subtask: Collate RDP and USEARCH de novo discovery results
###############################################################################
@jobs_limit(1, 'RGlobal')
@follows(run_usearch_denovo_discovery)
@transform([runMothurRDP_full_denovo_discovery,
            runMothurRDP_otu_denovo_discovery],
           regex('(.+)/(.+).greengenes_taxonomy.wang.taxonomy'),
           add_inputs(r'\1/\2_gg_matches.b6'),
           r'\1/\2_classification_summary.tsv')
def collateDeNovoDiscoveryResults(infiles, outfile):
    ''' '''
    in_rdp, in_b6 = infiles
    match_dict = collections.defaultdict(list)

    # summarize RDP results
    for line in open(in_rdp):
        line = line.split()
        tax = line[1].strip(';')
        tax = tax.split(';')
        match_dict[line[0]].append( tax.pop() )

    # summarize USEARCH results
    for line in open(in_b6):
        line = line.split()
        classification = line[1]
        classification = re.sub('\*', 'unclassified', classification)
        classification = classification + '(' + line[2] + ')'
        match_dict[line[0]].append(classification)

    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('Query_ID\tRDP_Classification\tUSEARCH_Classification\n')
        for k, v in match_dict.items():
            outf.write(k + '\t' + '\t'.join(v) + '\n')

    # summarize classification
    rdp = 0
    usearch = 0
    both = 0
    uncl = re.compile('unclassified')
    for match in match_dict.values():
        if uncl.match(match[0]) and uncl.match(match[1]):
            both += 1
        elif uncl.match(match[0]):
            rdp += 1
        elif uncl.match(match[1]):
            usearch += 1
        else:
            continue

    # plot venn diagrams
    R.assign('rdp', rdp)
    R.assign('usearch', usearch)
    R.assign('both', both)

    outf_pdf = P.snip(outfile, '.tsv') + '.pdf'
    R('''
    require("VennDiagram")
    pdf('{outf_pdf}')
    draw.pairwise.venn(rdp+both, usearch+both, both, 
                           category=c('RDP', 'USEARCH'))
    dev.off()'''.format(**locals()))

###############################################################################
# Metatask: QUANTIFY MOCK COMMUNITY V1-V3 COVERAGE
###############################################################################
## Task: Map V1V3_reads to mock reference database using CrossMatch
###############################################################################
@follows(mkdir("fasta_cross_match_v1v3.dir"))
@transform(os.path.join(PARAMS['data_mock_v1v3'],
                        'Diaz_16S_Mock_community_Pool*L001.clean.fasta.gz'),
           regex('.+/.+_(Pool.)_.+'),
           r'fasta_cross_match_v1v3.dir/\1/Mock\1.fasta.gz')
def fetchMockV1V3Data(infile, outfile):
    '''Retrieve the Illumina V1-V3 sequence data for the mock community'''
    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    statement = "zcat {infile}  | gzip > {outfile}"
    P.run()


@subdivide(fetchMockV1V3Data,
           regex('(.+)/(.+).fasta.gz'),
           r'\1/\2_*.fasta.gz')
def splitV1V3Fasta(infile, outfiles):
    '''Split the fastas for faster crossmatch mapping'''
    outfile_stub = P.snip(infile, '.fasta.gz')
    multiSplitFasta(infile, outfile_stub, 10)


@transform(splitV1V3Fasta,
           suffix('.fasta.gz'),
           add_inputs(PARAMS['data_mock_ref_db']),
           '_cross_match.txt')
def runCrossMatch_v1v3(infiles, outfile):
    '''Run cross match to map the Illumina V1V3 data to mock reference'''
    fasta_file, ref_db = infiles

    out_dir = os.path.dirname(outfile)
    # cross_match can't handle gzipped file format
    tmp_inf = P.getTempFilename(out_dir)

    
    # The minscore for mapping V1-V9 is approx 1/2 sequence length
    statement = ("zcat {fasta_file} > {tmp_inf};"
                 " cross_match"
                 "  -discrep_lists"
                 "  -tags"
                 "  -minscore 250"
                 "  -masklevel 0"
                 "  -penalty -9"
                 " {tmp_inf}"
                 " {ref_db}"
                 " > {outfile}")
    P.run()

    # cross_match doesn't give an option to rename logfiles
    shutil.move(tmp_inf + '.log', outfile + '.log')
    os.unlink(tmp_inf)

    
@collate(runCrossMatch_v1v3,
         regex('(.+)/(.+)/(.+)_[0-9]+_cross_match.txt'),
         r'\1/\3.txt')
def collateCrossMatchResults_v1v3(infiles, outfile):
    '''Concatenate all the crossmatch output into a single file
    '''
    infiles = ' '.join(infiles)
    to_cluster = False
    statement = ('cat {infiles} > {outfile}')
    P.run()


@transform(collateCrossMatchResults_v1v3,
           suffix('.txt'),
           '_readsMapped.txt')
def countReadMappingSuccess_v1v3(infile, outfile):
    '''Count the number of reads mapped'''
    mapped = set()
    for line in IOTools.openFile(infile):
        if line.startswith('ALIGNMENT'):
            mapped.add(line.split()[5])
    file_name = P.snip(infile, '.txt', strip_path=True)
    outf = IOTools.openFile(outfile, 'w')
    outf.write(file_name + '\t' + str(len(mapped)) + '\n')


@transform(collateCrossMatchResults_v1v3,
           suffix('.txt'),
           add_inputs(PARAMS['data_mock_ref_db']),
           '_summary.txt')
def summarizeCrossMatch_v1v3(infiles, outfile):
    '''parse crossmatch output and write summary tables.
    For output format see  www.phrap.org/phredphrap/phrap.html
    '''

    infile, ref_db = infiles

    # WARNING!
    if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
        L.warn('There is no sanity check for duplicate matches in'
               ' default cross_match run') 

    # get the maximum length of the sequences in the reference database
    # (informs output table size)
    max_len = 0
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(ref_db)):
        if len(fasta.sequence) > max_len:
            max_len = len(fasta.sequence)

    # parse cross match output
    out_dict = P16S.parseCrossMatchDiscrepList(infile, max_len, outfile)

    # pickle and dump out_dict
    out_pickle = P.snip(outfile, '.txt') + '.p'
    pickle.dump(out_dict, open(out_pickle, 'wb'))


@transform(summarizeCrossMatch_v1v3,
           suffix('_summary.txt'),
           '_readMappingSuccess.txt')
def summarizeCrossMatchMappingSuccess_v1v3(infile, outfile):
    '''read through the cross match output and count the number of reads
    successfully mapped to each reference sequence'''
    
    # a list of the target sequence IDs
    targets  = []
    header = True
    for line in IOTools.openFile(infile):
        if header:
            header = False
            assert line.split()[9] == 'TargetSeqID', "Problem with infile format"
            continue
        targets.append(line.split()[9])

    outf = IOTools.openFile(outfile, 'w')
    outf.write('TargetSequenceID\tReadsMapped\n')
    for seqID in set(targets):
        cnt = targets.count(seqID)
        outf.write(seqID + '\t' + str(cnt) + '\n')
    outf.close()


@jobs_limit(1, 'RGlobal')
@transform(summarizeCrossMatchMappingSuccess_v1v3,
           suffix('.txt'),
           add_inputs(loadRefDBTaxSummary),
           '.png')
def plotCrossMatchMappingSuccess_v1v3(infiles, outfile):
    '''plot number of reads successfully mapped to each reference 
    sequence
    '''
    map_table, ref_table = infiles
    ref_table = P.snip(os.path.basename(ref_table), '.load')
    ref_table = re.sub('-', '_', ref_table)

    df_map = pd.read_table(map_table, index_col=0)

    statement = 'SELECT StrainID,StrainName FROM {}'.format(ref_table)
    df_ref = pd.read_sql_query(statement,
                               con=connect(),
                               index_col='StrainID')
    sp = df_ref['StrainName'].tolist()
    sp = ['_'.join(x.split()[:2]) for x in sp]
    df_ref['SpName'] = sp

    df = df_map.join(df_ref)

    R('''rm(list=ls())''')
    R.assign('df', df)
    R('''png('{}')
         par(mar=c(15.1,4.1,4.1,2.1))
         barplot(df$ReadsMapped,
                 names.arg=df$SpName,
                 las=2)
         dev.off()'''.format(outfile))


###############################################################################
# Metatask: RELATING ERROR TO NUMBER OF PASSES
###############################################################################
## Task: Calculate number of passes for each read
###############################################################################
@files(os.path.join(PARAMS['data_mock_fastq_dir'],
                    'filtered_subread_summary.csv'),
           'error_vs_passes/read_pass_number.tsv')
def calculateNumberOfPasses(infile, outfile):
    """
    For every sequence read in pacbio run, find the number of times the 
    insert read was passed.
    """
    df = pd.read_table(infile, sep=',')

    # extract column containing hole numbers
    holes = df['HoleNumber']
    # count occurrences of each hole number in column
    summary = holes.value_counts()
    summary.name = 'NumberOfPasses'

    summary.to_csv(outfile, header=True, index_label='HoleNumber', sep='\t')
    

@transform(calculateNumberOfPasses, regex('.+/(.+).tsv'), r'\1.load')
def loadNumberOfPasses(infile, outfile):

    table_name = P.snip(os.path.basename(infile), '.tsv')

    to_cluster=False
    statement = ("cat %(infile)s |"
                 " python %(scriptsdir)s/csv2db.py"
                 "  --table=%(table_name)s"
                 " > %(outfile)s")
    P.run()

###############################################################################
@follows(summarizeCrossMatchMappingSuccess, 
         plotCrossMatchAlignmentSummary,
         loadCrossMatchSummary, 
         plotErrorHistograms)
def runMockCommunityAnalysis():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
