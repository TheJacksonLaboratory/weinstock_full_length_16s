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
========================
pipeline_pb_strains.py

A pipeline to characterize strain isolates based on PacBio 16s ccs

:Author: Jethro Johnson
========================

A pipeline to characterize strain isolates based on PacBio 16s ccs

INPUTS:

OUTPUT:

TO RUN:

MODULES: 
module load usearch
module load samtools
"""

import os
import sys
import re
import shutil
import collections
import itertools
import pickle
import gzip
import logging as L
import pandas as pd
import numpy as np
import rpy2.robjects as robjects

from Bio import SeqIO
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FI

import ruffus
from ruffus import *

import CGATPipelines.Pipeline as P


###############################################################################
# Read essential configuration parameters from pipeline.ini
###############################################################################
# location of config file
PARAMS = P.getParameters(["%s/pipeline.ini" % os.path.splitext(__file__)[0],
                          "pipeline.ini"],)


# location of input fasta files
try:
    DATADIR = PARAMS['input']
except NameError:
    DATADIR = "."

    
###############################################################################
# Utility functions
###############################################################################
def connect():
    """Connect to default sqlite database, returning a database handle
    that can be used for database IO. See pandas.DataFrame.to_sql() and
    and pandas.read_sql_query()
    """
    dbh = sqlite3.connect(PARAMS['database']) # name of database is in INI
    
    return dbh


###############################################################################
# Main Pipeline
###############################################################################
## Section I: Reorient sequences and collate technincal replicates
###############################################################################
@follows(mkdir('01_reoriented_fastas.dir'))
@transform(os.path.join(DATADIR, '*.fasta'),
           regex('.+/(.+)\.fasta'),
           r'01_reoriented_fastas.dir/\1.fasta')
def reorientFastas(infile, outfile):
    '''PacBio sequences have no particular orientation. Use mothur
    to ensure 5'->3' is V1->V9'''


    align_file = P.snip(infile, '.fasta') + '.align'
    out_file = P.snip(infile, '.fasta') + '.ng.fasta'
    
    statement=("mothur \"#align.seqs(candidate=%(infile)s,"
               " template=%(database_silva_gold)s, flip=t, threshold=0.50);"
               " degap.seqs(fasta=%(align_file)s)\"")
    P.run()

    shutil.move(os.path.abspath(out_file), outfile)


@follows(mkdir('02_renamed_fastas.dir'))
@subdivide(reorientFastas,
           regex('.+/(.+).fasta'),
           add_inputs(PARAMS['database_sample_names']),
           r'02_renamed_fastas.dir/\1_*.fasta')
def renameFastas(infiles, outfile):
    '''Rename isolate fasta files to include the original sample names'''

    fasta, sample_names = infiles

    # create a dictionary of sample names
    d = {}
    header = True
    for line in open(sample_names):
        if header:
            header = False
            continue
        line = line.split()
        isolate = line[0]
        sample_name = line[2]
        sample_name = re.sub('-|\.', '_', sample_name)
        
        assert isolate not in d.keys()
        d[isolate] = sample_name

    # generate the outfile name
    isolate = P.snip(fasta, '.fasta', strip_path=True)
    sample_name = d[isolate]
    out_file = os.path.join('02_renamed_fastas.dir',
                            isolate + '_' + sample_name + '.fasta')

    # iterate over the fasta file and change the headers to the isolate id
    n = 0
    with open(out_file, 'w') as outf:
        for fasta in SeqIO.parse(open(fasta), 'fasta'):
            n += 1
            outf.write('>' + isolate + '_' + str(n) + '\n' \
                       + str(fasta.seq) + '\n')


@follows(mkdir('03_collapsed_fastas.dir'))
@collate(renameFastas,
         regex('.+/isolate[0-9]+_(.+).fasta'),
         r'03_collapsed_fastas.dir/\1.fasta')
def collapseSampleFastas(infiles, outfile):
    '''Collapse all the files belonging to a single isolate'''

    infiles = ' '.join(infiles)
    statement = "cat %(infiles)s > %(outfile)s"
    P.run()

   
@follows(mkdir('04_dereplicated_fastas.dir'))
@transform(collapseSampleFastas,
           regex('.+/(.+)\.fasta'),
           r'04_dereplicated_fastas.dir/\1_dereplicated.fasta.gz')
def dereplicateSampleFastaFiles(infile, outfile):
    '''Use usearch to dereplicate the input fasta files'''

    outf_tmp = P.snip(outfile, '.gz')
    cluster_options = '-l walltime=00:20:00'
    statement = (" usearch"
                 "  -derep_fulllength %(infile)s"
                 "  -fastaout %(outf_tmp)s"
                 "  -sizeout;"
                 " gzip %(outf_tmp)s")
    P.run()


@transform(dereplicateSampleFastaFiles, suffix('.fasta.gz'), '_sorted.fasta.gz')
def sortSampleFastaFiles(infile, outfile):
    '''Use usearch to sort the input fasta files'''

    outf_tmp = P.snip(outfile, '.gz')
    tmpf = P.getTempFilename('.')
    cluster_options = '-l walltime=00:20:00'
    statement = ("gunzip -c %(infile)s > %(tmpf)s;"
                 " usearch"
                 "  -sortbysize %(tmpf)s"
                 "  -fastaout %(outf_tmp)s"
                 "  -minsize 1;" # keep all the sequences
                 " gzip %(outf_tmp)s")
    P.run()
    os.unlink(tmpf)


###############################################################################
## Section II: Align reads to most abundant sequence
###############################################################################
@follows(mkdir('05_mapped_fastas.dir'))
@transform(sortSampleFastaFiles,
           regex('.+/(.+)_dereplicated_sorted.fasta.gz'),
           r'05_mapped_fastas.dir/\1_template.fasta')
def fetchMostAbundantSampleSequence(infile, outfile):
    '''Select the unique sequence with the most number of replicates'''
    # the first sequence is the most abundant
    outf = open(outfile, 'w')
    out_id = P.snip(outfile, '_template.fasta', strip_path=True)
    for fasta in SeqIO.parse(gzip.open(infile), 'fasta'):
        outf.write('>' + str(out_id) + '\n' + str(fasta.seq) + '\n')
        break
    outf.close()


@merge(fetchMostAbundantSampleSequence,
       '05_representative_sequences.fasta.gz')
def collateMostAbundantSampleSequences(infiles, outfile):
    '''Combine single most abundant sequence for each sample into fasta
    '''
    infiles = " ".join(infiles)
    statement = ("cat %(infiles)s | gzip > %(outfile)s")
    to_cluster = False
    P.run()


@follows(fetchMostAbundantSampleSequence)
@transform(fetchMostAbundantSampleSequence,
           regex('(.+)/(.+)_template.fasta'),
           add_inputs(r'03_collapsed_fastas.dir/\2.fasta'),
           r'\1/\2_cross_match.txt')
def mapReadsToMostAbundantSampleSequence(infiles, outfile):
    '''Map all reads to most abundant sequence using CrossMatch'''
    template_fasta, query_fasta = infiles
    cluster_options = '-l walltime=02:00:00'
    statement = ("cross_match"
                 "  -discrep_lists"
                 "  -tags"
                 "  -masklevel 0"
                 "  %(cross_match_options)s"
                 " %(query_fasta)s"
                 " %(template_fasta)s"
                 " > %(outfile)s")
    P.run()


@jobs_limit(1)
@transform(mapReadsToMostAbundantSampleSequence,
           regex('(.+)/(.+)_cross_match.txt'),
           add_inputs(r'\1/\2_template.fasta'),
           r'\1/\2_cross_match_summary.txt')
def summarizeMapping(infiles, outfile):
    '''Create a summary dataframe out of the cross match results'''
    infile, ref_db = infiles


    def _parseCrossMatchDiscrepList(infile, max_len, outfile):
        '''Takes a discrepency list output by cross match (infile), the
        maximum length of the sequences in the reference database to which
        reads were mapped, and the name of the outfile.
        Writes a summary of the mismatches in each aligned read to outfile, 
        and returns a nested dictionary of mistmatches that can be converted
        into a multi-indexed dataframe.
        '''

        # set up out dict
        out_dict = {}
        outf = open(outfile, 'w')
        outf.write('Score\tSNP\tINS\tDEL\tQuerySeqID\tQueryStart\tQueryEnd\tQueryTail'
                   '\tReverseComp\tTargetSeqID\tTargetStart\tTargetEnd\tTargetTail\n')


        for line in open(infile):
            if line.startswith('ALIGNMENT'):
                line = line.split()
                length = len(line)
                # assign alginment values to separate IDs
                if length == 13:
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, TargetSeqID, TargetStart, TargetEnd, TargetTail = line[1:]
                    ReverseComp = 'F'
                elif length == 14 and line[-1:] == '*':
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, TargetSeqID, TargetStart, TargetEnd, TargetTail = line[1:-1]
                    ReverseComp = 'F'
                elif len(line) == 14:
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, ReverseComp, TargetSeqID, TargetTail, TargetEnd, TargetStart = line[1:]
                elif len(line) == 15 and line[-1:] == '*':
                    Score, SNP, INS, DEL, QuerySeqID, QueryStart, QueryEnd, \
                        QueryTail, ReverseComp, TargetSeqID, TargetTail, TargetEnd, TargetStart = line[1:-1]
                else:
                    raise ValueError('Expecting either a 13 or 14 field line'
                                     ' output from cross_match ALIGNMENT')

                TargetTail = re.sub('\(|\)', '', TargetTail)
                QueryTail = re.sub('\(|\)', '', QueryTail)

                outf.write('\t'.join([Score, SNP, INS, DEL, QuerySeqID, QueryStart, 
                                      QueryEnd, QueryTail, ReverseComp, TargetSeqID,
                                      TargetStart, TargetEnd, TargetTail]) + '\n')

                # sanity check
                if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
                    pass
                else:
                    assert (TargetSeqID, QuerySeqID) not in out_dict.keys(), \
                        ("There are multiple matches for the same sequence to the same"
                         " target in cross_match output: %s %s " % (TargetSeqID, QuerySeqID))
                    
                # reset the SNP/INDEL record and set dict entry
                errors = ['']*max_len
                out_dict[(TargetSeqID, QuerySeqID)] = errors

            elif line.startswith('DISCREPANCY'):
                # $2 == type of discrepency, $5 == discrepency location (1-based)
                # WARNING: All indels are marked as I/D... the number of bases is
                # ignored
                disc_loc = int(line.split()[4]) - 1 # table is 0-based
                expr = re.compile('-[0-9]+')
                disc_type = line.split()[1]
                disc_type = re.sub(expr, '', disc_type)


                # update the SNP/INDEL record
                out_dict[(TargetSeqID, QuerySeqID)][disc_loc] = disc_type

            else:
                continue

        return out_dict


    # WARNING!
    if os.path.basename(outfile).endswith('_cross_match_summary.txt'):
        L.warn('There is no sanity check for duplicate matches in'
               ' default cross_match run') 

    # get the maximum length of the sequences in the reference database
    # (informs output table size)
    max_len = 0
    for fasta in SeqIO.parse(open(ref_db), 'fasta'):
        if len(str(fasta.seq)) > max_len:
            max_len = len(str(fasta.seq))

    # parse cross match output
    out_dict = _parseCrossMatchDiscrepList(infile, max_len, outfile)

    # pickle and dump out_dict
    out_pickle = P.snip(outfile, '.txt') + '.p'
    pickle.dump(out_dict, open(out_pickle, 'wb'))


@follows(mkdir('06_technical_replicates.dir'))
@subdivide(summarizeMapping,
           regex('.+/(.+)_cross_match_summary.txt'),
           r'06_technical_replicates.dir/isolate*_\1_cross_match_summary.tsv')
def splitTechnicalReplicates(infile, outfiles):
    '''Take the crossmatch summary data and split out the technical
    replicates for each sample.'''

    out_dir = '06_technical_replicates.dir'
    infile = P.snip(infile, '.txt') + '.p'
    sample_id = P.snip(os.path.basename(infile), '_cross_match_summary.p')

    df = pd.DataFrame(pickle.load(open(infile, 'rb')))

    if df.empty:
        L.warn('No reads mapped successfully for sample: %s' % sample_id)
        # Creat empty output file for each isolate, to be filtered in next step
        isolates = open(PARAMS['database_sample_names']).readlines()[1:]
        sample_isolates = []
        for i in isolates:
            isolate = i.split()[0]
            s = i.split().pop()
            s = re.sub('-|\.', '_', s)
            if s == sample_id:
                sample_isolates.append(isolate)
        for isolate in sample_isolates:
            out_file = isolate + '_' + sample_id + '_cross_match_summary.tsv'
            outfile = os.path.join(out_dir, out_file)
            pd.DataFrame().to_csv(outfile, sep='\t')
    else:
        # Output a mapping summary file for each isolate
        df = df[sample_id]

        isolates = set([x.split('_')[0] for x in df.columns])
        for isolate in isolates:
            df_sub = df.filter(regex=isolate + '_*', axis=1)
            out_file = isolate + '_' + sample_id + '_cross_match_summary.tsv'
            outfile = os.path.join(out_dir, out_file)
            df_sub.to_csv(outfile, sep='\t')

   
###############################################################################
## Section III: Subdivide and filter each technical replicate
###############################################################################
@jobs_limit(1)
@follows(mkdir('07_filtered_technical_replicates.dir'))
@subdivide(splitTechnicalReplicates,
           regex('.+/(.+)_cross_match_summary.tsv'),
           r'07_filtered_technical_replicates.dir/*/\1_true_snps.tsv')
def filterTechnicalReplicates(infile, outfile):
    '''Parse the crossmatch output to select only samples that meet 
    minimum sequencing depth requirements, and have fewer than n% of 
    their bases divergent.
    '''
    passed = os.path.join('07_filtered_technical_replicates.dir', 'passed')
    if not os.path.exists(passed): os.mkdir(passed)
    failed = os.path.join('07_filtered_technical_replicates.dir', 'failed')
    if not os.path.exists(failed): os.mkdir(failed)
    
    def _getOutliers(df, x):
        '''Dataframe containing Percent, and scaling factor for outliers'''
        df1 = df[df['Percent']!=0]
        if df.empty:
            return 0
        variants = df1['Percent']

        upper = np.percentile(variants, 75) + x*np.subtract(*np.percentile(variants, [75, 25]))
        return upper

    isolate = P.snip(infile, '_cross_match_summary.tsv', strip_path=True)
    
    df = pd.read_table(infile, index_col=0, sep='\t')

    if df.empty:
        outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
        outf.write(isolate + '\tCoverage too low: 0\n')
        outf.close()
    else:
        df = df.applymap(lambda x: np.nan if x == '' else x)

        # filter on coverage
        sequence_depth = len(df.columns)
        if sequence_depth < int(PARAMS['filtering_min_coverage']):
            outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
            outf.write(isolate + '\tCoverage too low: %i\n' % sequence_depth )
            outf.close()
        else:
            # Drop the errors that are not substitutions
            df_sub = df.applymap(lambda x: np.nan if x in ['D', 'I'] else x)
            df_sub = df_sub.count(axis=1).to_frame().reset_index()
            df_sub.columns = ['Locus', 'Count']

            # Calculate the errors as a percentage
            df_sub['Percent'] = df_sub['Count']/float(sequence_depth) * 100 

            # Calculate the threshold below which to consider SNPs as noise
            # Either extreme outliers...
            threshold = _getOutliers(df_sub,
                                     PARAMS['filtering_outlier_scaling_factor'])
            # ...or some informed threshold 
            threshold = max([threshold, int(float(PARAMS['filtering_noise_threshold']))])
        
            # Calculate the proportion of bases that are above frequency threshold
            snps = sum(df_sub['Percent'] > threshold) / float(len(df_sub.index)) * 100
            # Threshold at which to consider alignment multiple sequences...
            if snps >= float(PARAMS['filtering_identity_threshold']):
                outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
                outf.write(isolate + '\tSequence divergence too high: %f\n' % snps)
            else:
                # Fetch true snp loci
                # Caculate extreme outliers...
                # ...discard outliers that are below some informed threshold.
                true_snp_loci = tuple(df_sub[df_sub['Percent']>threshold]['Locus'])
                true_snp_ct =  tuple(df_sub[df_sub['Percent']>threshold]['Count'])
                true_snps = [','.join(x) for x in zip(map(str, true_snp_loci),
                                                      map(str, [x for x in true_snp_ct]))]
        
                # Write true snps to flat file
                outf = open(os.path.join(passed, isolate + '_true_snps.tsv'), 'w')
                outf.write('\t'.join(map(str, true_snps)))
                outf.close()


@transform(filterTechnicalReplicates,
           regex('(.+)/passed/(isolate[0-9]{4})_(.+)_true_snps.tsv'),
           add_inputs(r'05_mapped_fastas.dir/\3_template.fasta',
                       r'05_mapped_fastas.dir/\3_cross_match.txt'),
           r'\1/\2_\3_denoised.fasta')
def fetchFilteredReplicateFasta(infiles, outfile):
    '''Filter the cross_match discrepancies based on whether they are likely 
    to be genuine, or noise. Output the filtered sequence for each alignment.
    '''

    def _cross_match_iterator(infile):
        query, snps = None, []
        
        for line in open(infile):
            if not (line.startswith("ALIGNMENT") or line.startswith('DISCREPANCY')):
                continue
            
            if line.startswith('ALIGNMENT'):
                if query != None:
                    yield query, snps
                query = line.split()[5]
                snps = []

            if line.startswith('DISCREPANCY'):
                d = line.split()[1]
                if not d.startswith('S'):
                    continue
                
                assert len(d) == 1, 'Unexpected discrepancy" %s' % d
                base = line.split()[3][0]
                position = int(line.split()[4]) - 1
                snps.append((base, position))
        else:
            yield query, snps
                
    true_snps, template, alignment_file = infiles

    # Fetch the isolate ID
    isolate = os.path.basename(true_snps).split('_')[0]
    
    # Parse the true snps into a tuple
    true_snp_loci = open(true_snps).readline().split()
    true_snp_loci = tuple([int(x.split(',')[0]) for x in true_snp_loci])

    # Fetch the template fasta seqeuence
    template = str(next(SeqIO.parse(open(template), 'fasta')).seq.upper())
        
    # Parse the alignment file, outputting corrected sequences
    with open(outfile, 'w') as outf:
        for read, alignment in _cross_match_iterator(alignment_file):
            # Reads from all technical replicates appear in the alignment file
            if not read.startswith(isolate):
                continue

            seq_out = [i for i in template]
            
            for base, position in alignment:
                if position in true_snp_loci:
                    seq_out[position] = base

            outf.write('>' + read + '\n' + ''.join(seq_out) + '\n')
            


@follows(mkdir('08_filtered_16S_gene_copies.dir'))
@transform(fetchFilteredReplicateFasta,
           regex('.+/(.+)_denoised.fasta'),
           r'08_filtered_16S_gene_copies.dir/\1_unique.fasta')
def dereplicateFilteredFastaFiles(infile, outfile):
    
    cluster_options = '-l walltime=00:20:00'
    statement = (" usearch"
                 "  -derep_fulllength %(infile)s"
                 "  -fastaout %(outfile)s"
                 "  -sizeout"
                 " &> %(outfile)s.log")

    P.run()


###############################################################################
## Section IV: Summarize technical variability
###############################################################################
@follows(mkdir('09_technical_variability'))
@collate(filterTechnicalReplicates,
         regex('.+/passed/(isolate[0-9]{4})_(.+)_true_snps.tsv'),
         r'09_technical_variability/\2.tsv')
def collateTechnicalReplicates(infiles, outfile):
    '''Generate a dataframe containing the relative frequency of SNPs at
    different loci in each technical replicate.
    '''

    snp_dict = collections.defaultdict(dict)

    for infile in infiles:

        # Fetch the dataframe summarizing read mapping.
        sample_id = P.snip(infile, '_true_snps.tsv', strip_path=True)
        map_file = os.path.join('06_technical_replicates.dir',
                                sample_id + '_cross_match_summary.tsv')
        assert os.path.exists(map_file), map_file
        
        # Number of reads corresponds to the number of columns in map file
        n_reads = len(open(map_file).readline().split('\t')[1:])

        # read in the list of snp loci
        snps = open(infile).readline().split()
    
        # combine the isolate Id and the read depth to generate column headers
        new_sample_id = sample_id.split('_')[0] + '_' + str(n_reads)
        
        # add each snp frequency to dictionary
        for snp in snps:
            locus, frequency = snp.split(',')
            frequency = float(frequency) / float(n_reads)
            
            snp_dict[new_sample_id][locus] = frequency
            
    df = pd.DataFrame(snp_dict)
    df.fillna(0, inplace=True)
    df.index = [int(x) for x in df.index]
    df.sort_index(axis=0, inplace=True)
    df.index.name = 'locus'
    
    df.to_csv(outfile, sep='\t')


@follows(mkdir('10_pairwise_measurement_error'))
@jobs_limit(1)
@transform(collateTechnicalReplicates,
           regex('.+/(.+).tsv'),           
           r'10_pairwise_measurement_error/\1.tsv')
def calculateMeasurementError(infile, outfile):
    '''For each pair of technical replicates, calculate the measurement
    error across variabile loci.
    '''

    R('''rm(list=ls())''')

    zeta_w = R('''
               zeta_w <- function(df, y="value", x="locus"){
               res = anova(lm(value ~ locus, data=df))
               return(sqrt(res[["Mean Sq"]][2]))
               }
    ''')


    # Hack... I forgot about samples with only a single locus
    zeta_w2 = R('''
                zeta_w2 <- function(df){
                v = apply(df, 1, var)
                m = mean(v)
                return(sqrt(m))
                }
    ''')
    
    # Open the dataframe and check that there is more than one measurement...
    df = pd.read_table(infile, sep='\t', index_col=0)
    if df.empty:
        L.warn('Sample %s has empty dataframe' % os.path.basename(infile))
    elif len(df.columns) < 2:
        L.warn('Sample %s has only one sequencing measurement' % os.path.basename(infile))
    else:
        sample_id =  P.snip(infile, '.tsv', strip_path=True)
        outf = open(outfile, 'w')
        outf.write('SampleID\tIsolateComparison\tMeanDepth\tDepthDiff\tMeasurementError\n')

        for pair in itertools.combinations(df.columns, 2):

            # subset dataframe
            pair = list(pair)
            df1 = df[pair]

            if len(df.index) == 1:
                L.warn('Sample %s has only one variable lcous' % \
                       os.path.basename(infile))
                m_error = zeta_w2(df)[0]
            else:
                # melt the dataframe
                df1['locus'] = [str(x) for x in df1.index]
                df1 = df1.melt(id_vars='locus')
                print(df)
                
                # calculate the measurement error
                m_error = zeta_w(df1)[0]

            # Fetch the contrast, mean read depth and difference in read depth
            p1 = pair[0]
            v1 = float(p1.split('_')[1])
            p1 = p1.split('_')[0]

            p2 = pair[1]
            v2 = float(p2.split('_')[1])
            p2 = p2.split('_')[0]
            contrast = p1 + '_' + p2

            md = np.mean([v1, v2])
            dd = abs(v1 - v2)
            
            outf.write('\t'.join(map(str, [sample_id,
                                           contrast,
                                           md,
                                           dd,
                                           m_error])) + '\n')
            
        outf.close()
    

@follows(collateTechnicalReplicates)
@merge('09_technical_variability/*tsv',
       '10_pairwise_measurement_error/pairwise_measurement_error.summary')
def summarizeMeasurementErrorCalculations(infiles, outfile):
    '''Summarize the for which pairwise measurement error can be calculated'''

    with open(outfile, 'w') as outf:
        for infile in infiles:
            df = pd.read_table(infile, sep='\t', index_col=0)
            if df.empty:
                outf.write('%s\tNo reads aligned\n' % os.path.basename(infile))
            elif len(df.columns) < 2:
                outf.write('%s\tNo replicates\n' % os.path.basename(infile))
            else:
                outf.write('%s\tHas replicates\t%i\n' % (os.path.basename(infile),
                                                         len(df.columns)))

        
@follows(mkdir('11_per_sample_measurement_error'))
@jobs_limit(1)
@transform(collateTechnicalReplicates,
           regex('.+/(.+).tsv'),           
           r'11_per_sample_measurement_error/\1.tsv')
def calculatePerSampleMeasurementError(infile, outfile):
    '''Sample, calculate the measurement error across variabile
    loci for all technical replicates in one
    '''

    R('''rm(list=ls())''')

    zeta_w = R('''
               zeta_w <- function(df, y="value", x="locus"){
               res = anova(lm(value ~ locus, data=df))
               return(sqrt(res[["Mean Sq"]][2]))
               }
    ''')

    # Hack... I forgot about samples with only a single locus
    zeta_w2 = R('''
                zeta_w2 <- function(df){
                v = apply(df, 1, var)
                m = mean(v)
                return(sqrt(m))
                }
    ''')
    
    # Open the dataframe and check that there is more than one measurement...
    df = pd.read_table(infile, sep='\t', index_col=0)
    if len(df.columns) < 2:
        pass
    else:

        # Fetch the sequencing depths
        depths = [float(x.split('_')[1]) for x in df.columns]
        mean_depth = np.mean(depths)

        sample_id =  P.snip(infile, '.tsv', strip_path=True)
        outf = open(outfile, 'w')
        outf.write('SampleID\tMeanDepth\tMeasurementError\n')

        if len(df.index) == 1:
            L.warn('Sample %s has only one variable lcous' % \
                   os.path.basename(infile))
            m_error = zeta_w2(df)[0]
        else:
            # melt the dataframe
            df['locus'] = [str(x) for x in df.index]
            df = df.melt(id_vars='locus')
            
            # calculate the measurement error
            m_error = zeta_w(df)[0]

        outf.write('\t'.join(map(str, [sample_id,
                                       mean_depth,
                                       m_error])) + '\n')
        outf.close()


@follows(calculateMeasurementError)
@merge('10_pairwise_measurement_error/*tsv',
       'pairwise_measurement_error.tsv')
def collapsePairwiseMeasurementError(infiles, outfile):

    with open(outfile, 'w') as outf:
        outf.write('SampleID\tMeanDepth\tDepthDiff\tMeasurementError\n')

    for infile in infiles:
        statement = "cat %(infile)s | awk 'NR != 1' >> %(outfile)s"
        to_cluster = False
        P.run()


#@follows(calculatePerSampleMeasurementError)
@merge('11_per_sample_measurement_error/*tsv',
       'persample_measurement_error.tsv')
def collapsePerSampleMeasurementError(infiles, outfile):

    with open(outfile, 'w') as outf:
        outf.write('SampleID\tMeanDepth\tMeasurementError\n')

    for infile in infiles:
        statement = "cat %(infile)s | awk 'NR != 1' >> %(outfile)s"
        to_cluster = False
        P.run()


@follows(collapsePairwiseMeasurementError,
         collapsePerSampleMeasurementError)
def collapseMeasurementError():
    pass


###############################################################################
## Section V: Filter the pooled samples
###############################################################################
@follows(mkdir('12_filtered_samples.dir'))
@subdivide(summarizeMapping,
           regex('.+/(.+)_cross_match_summary.txt'),
           r'12_filtered_samples.dir/*/\1_true_snps.tsv')
def filterSamples(infile, outfile):
    '''Apply exactly the same filtering used to clean up technical replicates, 
    but using the pooled data'''
    
    passed = os.path.join('12_filtered_samples.dir', 'passed')
    if not os.path.exists(passed): os.mkdir(passed)
    failed = os.path.join('12_filtered_samples.dir', 'failed')
    if not os.path.exists(failed): os.mkdir(failed)
    
    def _getOutliers(df, x):
        '''Dataframe containing Percent, and scaling factor for outliers'''
        df1 = df[df['Percent']!=0]
        if df.empty:
            return 0
        variants = df1['Percent']

        upper = np.percentile(variants, 75) + x*np.subtract(*np.percentile(variants, [75, 25]))
        return upper

    infile = P.snip(infile, '.txt') + '.p'
    assert os.path.exists(infile)
    isolate = P.snip(infile, '_cross_match_summary.p', strip_path=True)
    df = pd.DataFrame(pickle.load(open(infile, 'rb')))

    to_remove = PARAMS['filtering_samples_to_remove'].split(',')
    
    if df.empty:
        outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
        outf.write(isolate + '\tCoverage too low:\t0\n')
        outf.close()
    elif isolate in to_remove:
        outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
        outf.write(isolate + '\Sample explicitly removed\n')
        outf.close()        
    else:
        df = df[isolate]
        df = df.applymap(lambda x: np.nan if x == '' else x)

        # filter on coverage
        sequence_depth = len(df.columns)
        if sequence_depth < int(PARAMS['filtering_min_coverage']):
            outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
            outf.write(isolate + '\tCoverage too low: %i\n' % sequence_depth )
            outf.close()
        else:
            # Drop the errors that are not substitutions
            df_sub = df.applymap(lambda x: np.nan if x in ['D', 'I'] else x)
            df_sub = df_sub.count(axis=1).to_frame().reset_index()
            df_sub.columns = ['Locus', 'Count']

            # Calculate the errors as a percentage
            df_sub['Percent'] = df_sub['Count']/float(sequence_depth) * 100 

            # Calculate the threshold below which to consider SNPs as noise
            # Either extreme outliers...
            threshold = _getOutliers(df_sub,
                                     PARAMS['filtering_outlier_scaling_factor'])
            # ...or some informed threshold 
            threshold = max([threshold, int(float(PARAMS['filtering_noise_threshold']))])
        
            # Calculate the proportion of bases that are above frequency threshold
            snps = sum(df_sub['Percent'] > threshold) / float(len(df_sub.index)) * 100
            # Threshold at which to consider alignment multiple sequences...
            if snps >= float(PARAMS['filtering_identity_threshold']):
                outf = open(os.path.join(failed, isolate + '_true_snps.tsv'), 'w')
                outf.write(isolate + '\tSequence divergence too high:\t%f\n' % snps)
            else:
                # Fetch true snp loci
                # Caculate extreme outliers...
                # ...discard outliers that are below some informed threshold.
                true_snp_loci = tuple(df_sub[df_sub['Percent']>threshold]['Locus'])
                true_snp_ct =  tuple(df_sub[df_sub['Percent']>threshold]['Percent'])
                true_snps = [','.join(x) for x in zip(map(str, true_snp_loci),
                                                      map(str, [x for x in true_snp_ct]))]
        
                # Write true snps to flat file
                outf = open(os.path.join(passed, isolate + '_true_snps.tsv'), 'w')
                outf.write('\t'.join(map(str, true_snps)))
                outf.close()    


@transform(filterSamples,
           regex('(.+)/passed/(.+)_true_snps.tsv'),
           add_inputs(r'05_mapped_fastas.dir/\2_template.fasta',
                       r'05_mapped_fastas.dir/\2_cross_match.txt'),
           r'\1/\2_denoised.fasta')
def fetchFilteredSampleFasta(infiles, outfile):
    '''Filter the cross_match discrepancies based on whether they are likely 
    to be genuine, or noise. Output the filtered sequence for each alignment.
    '''

    def _cross_match_iterator(infile):
        query, snps = None, []
        
        for line in open(infile):
            if not (line.startswith("ALIGNMENT") or line.startswith('DISCREPANCY')):
                continue
            
            if line.startswith('ALIGNMENT'):
                if query != None:
                    yield query, snps
                query = line.split()[5]
                snps = []

            if line.startswith('DISCREPANCY'):
                d = line.split()[1]
                if not d.startswith('S'):
                    continue
                
                assert len(d) == 1, 'Unexpected discrepancy" %s' % d
                base = line.split()[3][0]
                position = int(line.split()[4]) - 1
                snps.append((base, position))
        else:
            yield query, snps
                
    true_snps, template, alignment_file = infiles
    
    # Parse the true snps into a tuple
    true_snp_loci = open(true_snps).readline().split()
    true_snp_loci = tuple([int(x.split(',')[0]) for x in true_snp_loci])

    # Fetch the template fasta seqeuence
    template = str(next(SeqIO.parse(open(template), 'fasta')).seq.upper())
        
    # Parse the alignment file, outputting corrected sequences
    with open(outfile, 'w') as outf:
        for read, alignment in _cross_match_iterator(alignment_file):

            seq_out = [i for i in template]
            
            for base, position in alignment:
                if position in true_snp_loci:
                    seq_out[position] = base

            outf.write('>' + read + '\n' + ''.join(seq_out) + '\n')
            

@transform(fetchFilteredSampleFasta,
           regex('(.+)/(.+)_denoised.fasta'),
           r'\1/\2_unique.fasta')
def dereplicateFilteredSampleFastas(infile, outfile):
    
    cluster_options = '-l walltime=00:20:00'
    statement = (" usearch"
                 "  -derep_fulllength %(infile)s"
                 "  -fastaout %(outfile)s"
                 "  -sizeout"
                 " &> %(outfile)s.log")
    to_cluster=False
    P.run()

    
###############################################################################
## Generate 99% OTU clusters for the in silico reads that passed filtering
###############################################################################
@follows(mkdir('13_filtered_sample_otus.dir'))
@collate(filterSamples,
         regex('12_filtered_samples.dir/passed/(.+)_true_snps.tsv'),
         r'13_filtered_sample_otus.dir/isolate_otus.fasta')
def createOTUs(infiles, outfile):
    '''Make OTUs from the template sequences of all samples that passed
    filtering'''

    fasta_files = []
    for infile in infiles:
        fasta_file = P.snip(infile, '_true_snps.tsv', strip_path=True)\
                     + '_template.fasta'
        fasta_file = os.path.join('05_mapped_fastas.dir', fasta_file)
        assert os.path.exists(fasta_file)
        fasta_files.append(fasta_file)

    tmpf = P.getTempFilename('.')
    with IOTools.openFile(tmpf, 'w') as outf:
        for fasta_file in fasta_files:
            for fasta in FI.FastaIterator(IOTools.openFile(fasta_file)):
                outf.write('>' + fasta.title + ';size=1\n' + fasta.sequence + '\n')

    out_up = P.snip(outfile, '.fasta') + '_up.txt'
    statement = ("usearch"
                 "  -cluster_otus %(tmpf)s"
                 "  -otus %(outfile)s"
                 "  -uparseout %(out_up)s"
                 "  -otu_radius_pct 1.0"
                 "  -uparse_break -999"
                 "  -relabel OTU_"
                 " &> %(outfile)s.log")
    to_cluster = False
    P.run()


@subdivide(filterSamples,
           regex('(.+)/passed/(.+)_true_snps.tsv'),
           add_inputs(createOTUs),
           r'\1/OTU_*_\2.tsv')
def assignFilteredSamplesToOTUs(infiles, outfiles):
    '''Prepend OTU to sample name for filtered samples'''

    snp_file, otu_file = infiles

    # fetch sample to OTU mapping
    otu_file = P.snip(otu_file, '.fasta') + '_up.txt'
    otu_dict = {}
    for row in open(otu_file):
        sample_id = row.split()[0].split(';')[0]
        otu_id = row.split().pop()
        otu_dict[sample_id] = otu_id

    sample_id = P.snip(snp_file, '_true_snps.tsv', strip_path=True)
    otu = otu_dict[sample_id]
    out_dir = os.path.dirname(os.path.dirname(snp_file))
    out_file = os.path.join(out_dir, otu + '_' + sample_id + '.tsv')

    shutil.copyfile(snp_file, out_file)

    
@follows(mkdir('14_filter_sample_error_profiles.dir'))
@jobs_limit(1)
@subdivide(filterSamples,
           regex('.+/passed/(.+)_true_snps.tsv'),
           add_inputs(createOTUs),
           r'14_filter_sample_error_profiles.dir/*_\1.pdf')
def plotFilteredSamples(infiles, outfiles):
    '''Create a plot of the SNP profiles for each filtered sample'''

    error_profile, otu_assignment = infiles

    otu_assignment = P.snip(otu_assignment, '.fasta') + '_up.txt'
    otu_dict = {}
    for row in open(otu_assignment):
        sample_id = row.split()[0].split(';')[0]
        otu_id = row.split().pop()
        otu_dict[sample_id] = otu_id
    
    def _fetch_loci(infile):
        # Some samples have no snps...
        if not open(infile).readline():
            L.warn('Sample %s has no SNPs' % infile)
            idx = [i for i in range(1, 1501)]
            snp = [0,]*1500
            df = pd.DataFrame([idx, snp]).transpose()
        else:
            df = pd.DataFrame([x.split(',') for x in open(infile).readline().split('\t')])
            
        df.columns = ['Locus', 'Frequency']
        df = df.applymap(float)
        
        return df

    sample_id = P.snip(error_profile, '_true_snps.tsv', strip_path=True)
    otu_id = otu_dict[sample_id]
    outfile = os.path.join('14_filter_sample_error_profiles.dir',
                           otu_id + '_' + \
                           sample_id + '.pdf')

    R('''rm(list=ls())''')
    R('''require('ggplot2')''')
    df = _fetch_loci(error_profile)
    R.assign('df', df)

    R('''require('ggplot2')
         pl <- ggplot(df, aes(x=Locus, xend=Locus, y=0, yend=Frequency)) + geom_segment()
         pl <- pl + theme_bw() + theme(panel.grid=element_blank())
         pl <- pl + xlim(0, 1500) + scale_y_continuous(expand=c(0,0), limits=c(0, 100))
         pl <- pl + xlab('Position Along 16S Gene') + ylab('Frequency (%%)')
         pl <- pl + ggtitle('%s\n%s')
         pdf('%s', height=3, width=5)
         plot(pl)
         dev.off()
      ''' % (otu_id, sample_id, outfile))

    R('''rm(list=ls())''')


@follows(mkdir('15_otu_error_profiles.dir'))
@collate(assignFilteredSamplesToOTUs,
         regex('.+/(OTU_[0-9]+)_.+.tsv'),
         r'15_otu_error_profiles.dir/\1_error_profiles.tsv')
def collateOTUErrorProfiles(infiles, outfile):
    '''For each of the isolate snp profiles belonging to a single OTU...
    merge them together in a single table.
    '''

    def _fetch_loci(infile):
        # Some samples have no snps...
        if not open(infile).readline():
            L.warn('Sample %s has no SNPs' % infile)
            df = pd.DataFrame({'Locus': [], 'Frequency': []})
        else:
            df = pd.DataFrame([x.split(',') for x in open(infile).readline().split('\t')])
            df.columns = ['Locus', 'Frequency']
            df = df.applymap(float)

        df.set_index('Locus', inplace=True)
        df.columns = [P.snip(infile, '.tsv', strip_path=True),]
            
        return df

    # Fetch dataframes containing snp loci for each isolate in an OTU
    otu_error_profiles = []
    for infile in infiles:
        inf_df = _fetch_loci(infile)
        otu_error_profiles.append(inf_df)

        
    # Merge the loci into a single table
    df_out = otu_error_profiles.pop()
    if otu_error_profiles:
        for error_profile in otu_error_profiles:
            df_out = df_out.join(error_profile, how='outer')

    df_out.index = [int(x) for x in df_out.index]
    df_out.to_csv(outfile, sep='\t')

    
@jobs_limit(1)
@transform(collateOTUErrorProfiles,
           suffix('_error_profiles.tsv'),
           '_unique_profiles.tsv')
def fetchUniqueOTUErrorProfiles(infile, outfile):
    '''Collapse error profiles based on understanding of measurement error'''

    def _fetch_unique_isolates(df,
                               log_file,
                               summary_file,
                               threshold=6.58123249260359,):
        '''Receives dataframe in which columns are isolates, 
        rows are loci at which snps exist, returns a dataframe
        of unique isolates. Where unique are those which differ
        by > threshold at one or more locus.
        '''
        logf = open(log_file, 'w')
        sumf = open(summary_file, 'w')
        
        samples = df.columns.tolist()
        to_keep = []
        to_keep.append(samples.pop())
        
        logf.write('1:\tQuery %s is unique\n' % to_keep[0])
        
        n = 1
        for query in samples:
            n += 1
            matches = []
            for reference in to_keep:
                q = df[query]
                r = df[reference]
                if max(abs(r-q)) <= threshold:
                    logf.write('%i:\tQuery %s matches %s\n' % (n, query, reference))
                    matches.append(True)
                    break
                else:
                    matches.append(False)
            if not sum(matches):
                logf.write('%i:\tQuery %s is unique\n' % (n, query))
                to_keep.append(query)
                        
        sumf.write('%i\t%i\n' % (len(to_keep), n))

        logf.close()
        sumf.close()
        
        df_out = df[to_keep]
        return df_out
                    
    logfile = P.snip(infile, '.tsv') + '.log'
    sumfile = P.snip(infile, 's.tsv') + '_summary.tsv'

    df = pd.read_table(infile, sep='\t', index_col=0)
    df = df.fillna(value=0)

    # The possibility that none of the isolates have any snps
    if df.empty:
        df_out = pd.DataFrame(df.iloc[:,0])
        sumf = open(sumfile, 'w')
        sumf.write('1\t%i\n' % len(df.columns))
        sumf.close()
    else:
        df_out = _fetch_unique_isolates(df, logfile, sumfile)
    df_out.to_csv(outfile, sep='\t')


@merge(fetchUniqueOTUErrorProfiles, 'isolate_error_profiles.tsv')
def collateUniqueErrorProfiles(infiles, outfile):
    '''Create a single dataframe containing all the unique isolate profiles
    '''

    # Create a dummy dataframe with all 1500 bases present
    df = pd.DataFrame(columns=['Dummy',], index=range(1, 1501))

    for infile in infiles:
        df_tmp = pd.read_table(infile, index_col=0, sep='\t')
        df = df.join(df_tmp, how='outer')
        
    df.drop('Dummy', axis=1, inplace=True)
    df.fillna(value=0, inplace=True)
    df.index.name = 'Locus'

    df.to_csv(outfile, sep='\t')


@merge(fetchUniqueOTUErrorProfiles, 'isolate_collectors_curve.tsv')
def collateUniqueErrorProfileSummary(infiles, outfile):
    '''Collate the relationship between the number of isolates sequenced
    and the number of unique isolates for each OTU.
    '''

    infiles = [P.snip(x, '_unique_profiles.tsv') + '_error_profile_summary.tsv' \
               for x in infiles]

    infiles = ' '.join(infiles)

    statement = ("printf 'UniqueIsolates\tTotalIsolates\n' > %(outfile)s &&"
                 " cat %(infiles)s >> %(outfile)s")
    to_cluster = False
    P.run()
    
# ###############################################################################
# # Looking at technical variability. 
# ###############################################################################
# @follows(mkdir('renamed_fastas.dir'))
# @subdivide(fetchFilteredFasta,
#            regex('.+/(.+)_denoised.fasta'),
#            add_inputs(PARAMS['database_sample_names']),
#            r'renamed_fastas.dir/\1_*.fasta')
# def renameFastas(infiles, outfile):
#     '''Rename isolate fasta files to the original sample names'''

#     fasta, sample_names = infiles

#     # create a dictionary of sample names
#     d = %()s
#     header = True
#     for line in open(sample_names):
#         if header:
#             header = False
#             continue
#         line = line.split()
#         isolate = line[0]
#         sample_name = line[2]
#         sample_name = re.sub('-|\.', '_', sample_name)
        
#         assert isolate not in d.keys()
#         d[isolate] = sample_name
        
#     isolate = P.snip(fasta, '_denoised.fasta', strip_path=True)
#     sample_name = d[isolate]

#     out_file = os.path.join('renamed_fastas.dir',
#                             isolate + '_' + sample_name + '.fasta')

#     os.symlink(os.path.abspath(fasta), out_file)
    
    
# @follows(mkdir('sample_snp_profiles.dir'))
# @collate(renameFastas,
#          regex('.+/isolate[0-9]+_(.+).fasta'),
#          r'sample_snp_profiles.dir/\1.tsv')
# def collateSampleSNPProfiles(infiles, outfile):
#     '''For each sample, collate information on the true SNPs found 
#     in each technical replicate and output in a single dataframe
#     '''
#     sample_id = os.path.basename(outfile)[:-len('.tsv')]
    
#     snp_dict = %()s
#     for infile in infiles:
#         replicate_id = os.path.basename(infile).split('_')[0]

#         # Fetch the record of snp location,frequency
#         inf = os.path.join('filtered_fastas.dir', 'passed',
#                            replicate_id + '_true_snps.tsv')
        
#         d = %()s
#         snps = open(inf).readline().split()
#         for snp in snps:
#             locus, frequency = snp.split(',')
#             d[locus] = float(frequency)

#         snp_dict[replicate_id] = d
                        
#     df = pd.DataFrame(snp_dict)
#     df.fillna(value=0, inplace=True)

#     df.to_csv(outfile, sep='\t')
    
#     # While we're here... output the within locus stddev: (zeta)w...

#     # Ignoring samples with no replicates, and no snps
#     if len(infiles) > 1 and len(df.index) > 0:
#         outf = outfile[:-len('.tsv')] + '_measurement_error.tsv'
#         outf = open(outf, 'w')
        
#         # The sample variance at each locus
#         var = df.apply(lambda x: np.var(x, ddof=1), axis=1)
#         # The mean of these variance estimates
#         var_mean = np.mean(var)
#         # The square rood of the mean of the per-locus variance estimates
#         zeta_w = np.sqrt(var_mean)
    
#         # fetch the coverage for each sequenceing run
#         coverage = [len(open(x).readlines()) for x in infiles]

#         outf.write('\t'.join([sample_id,
#                               str(zeta_w),
#                               str(max(coverage)),
#                               str(min(coverage)),
#                               str(np.mean(coverage))]) + '\n')
#         outf.close()

        
#     # ...and the per locus std dev and snp frequency
#         outf = outfile[:-len('.tsv')] + '_per_locus_stdev.tsv'
#         outf = open(outf, 'w')
    
#         # The standard deviation at each locus
#         stdev = map(str, df.apply(lambda x: np.std(x, ddof=1), axis=1))
#         # The mean snp frequency at each locus
#         snp_freq = map(str, df.apply(np.mean, axis=1))

#         outf.write('\n'.join(['\t'.join(x) for x in zip(snp_freq, stdev)]))
#         outf.close()

    
if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
