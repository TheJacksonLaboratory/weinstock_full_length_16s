
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
import ftplib
import MySQLdb

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects import r as R

import CGATPipelines.Pipeline as P

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


def fetchBioSQLTaxonomy(ncbi_tax_id,
                        host='weinstock_dbs',
                        user='guest',
                        passwd='',
                        db='biosql'):
    '''Fetch the full parent taxonomy (to Kingdom) for
       a given NCBI tax id number'''

    phylogeny = collections.OrderedDict()

    # connect to BioSQL database
    db = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
    cur = db.cursor()

    # break if more than 50 iterations...
    safety = 50
    first = True
    node_rank = None
    while node_rank != 'superkingdom':
        safety -= 1
        if safety < 0:
            break

        # internal BioSQL taxon_ids don't correspond to ncbi taxon IDs
        # WARNING: There is overlap between the two. 
        field = 'taxon_id'
        if first:
            tax_id = ncbi_tax_id
            field = 'ncbi_taxon_id'
            first = False
            
        # fetch tax_id, parent_tax_id node_rank, scientific_name
        statement = ("SELECT "
                     "  taxon.taxon_id,"
                     "  taxon.parent_taxon_id,"
                     "  taxon.node_rank,"
                     "  taxon_name.name"
                     " FROM taxon JOIN taxon_name"
                     "  ON taxon.taxon_id = taxon_name.taxon_id"
                     " WHERE taxon.{} = {} "
                     "AND name_class = 'scientific name'".format(field, tax_id))
        cur.execute(statement)
        taxonomy = cur.fetchall()
        if len(taxonomy) == 1:
            old_tax_id, tax_id, node_rank, scientific_name = taxonomy[0]
            phylogeny[node_rank] = scientific_name
        elif len(taxonomy) == 0:
            E.warn('No taxonomy information for taxon ID %s beyond %s' \
                   % (ncbi_tax_id, node_rank))
            phylogeny['genus'] = 'unavailable'
            phylogeny['species'] = 'unavailable'
            
            return phylogeny

        else:
            raise ValueError('multiple entries for tax id %s: %s' \
                             % (ncbi_tax_id, str(taxonomy)))
            
    return phylogeny


###############################################################################
# Fetch Bacteroides 16S genes... 
###############################################################################
@follows(mkdir('taxon_info.dir'))
@files(PARAMS['rtg_taxonomy_file'], 'taxon_info.dir/all_taxa.tsv')
def fetchAllTaxa(infile, outfile):
    '''Fetch all taxa associated with input taxon id'''


    taxid2rank = {}
    taxid2parent = {}
    taxid2name = {}
    parent2taxid = collections.defaultdict(list)

    # create dictionaries for mapping info
    for line in IOTools.openFile(infile):
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        parent2taxid[line[1]].append(line[0])
        taxid2parent[line[0]] = line[1]
        taxid2rank[line[0]] = line[2]
        taxid2name[line[0]] = line[3]
        
    def _fetchTaxonomy(taxid, keep_parents=True, print_output=False):
        '''Receives a list of taxids. Outputs the taxid, the rank, 
        the name for this tax id and all child tax ids '''

        tax_list = parent2taxid[str(taxid)]
        tax_ids = set()
        
        i=0
        while i < len(tax_list):
            parent = tax_list[i]
            children = parent2taxid[parent]
            if children:
                tax_list.extend(children)
                if keep_parents:
                    tax_ids.add(parent)
            else:
                tax_ids.add(parent)
            i += 1
                
        if print_output:
            print 'TaxID', 'Rank', 'Name'
            for t in sorted([int(x) for x in tax_ids]):
                print str(t), taxid2rank[str(t)], taxid2name[str(t)]
        else:
            return map(str, sorted([int(x) for x in tax_ids]))

        
    with IOTools.openFile(outfile, 'w') as outf:
        outf.write('taxid\trank\tname\n')
        for tax_id in _fetchTaxonomy(PARAMS['taxon_id']):
            outf.write('\t'.join([tax_id,
                                  taxid2rank[tax_id],
                                  taxid2name[tax_id]]) + '\n')
        
                
@transform(fetchAllTaxa, suffix('.tsv'), '_genbank_ids.tsv')
def fetchTaxonGenBankAccessions(infile, outfile):
    '''Fetch a single genbank accession ID for each taxon id'''

    gb_dict = collections.defaultdict(list)
    # parse the entire taxonomy lookup table into a dictionary
    for line in IOTools.openFile(PARAMS['rtg_taxonomy_lookup_file']):
        t_id, accession  = line.strip().split('\t')
        gb_id = accession.split('|')[3]

        gb_dict[t_id].append(gb_id)

    pickle.dump(gb_dict, open('gb_dict.p', 'wb'))

    outf_discard = IOTools.openFile(outfile[:-len('.tsv')] + '_discarded.tsv', 'w')
    
    with IOTools.openFile(outfile, 'w') as outf:
        header = True
        for line in IOTools.openFile(infile):
            if header:
                header = False
                outf.write('gb_accession\ttaxid\trank\tname\n')
                continue
            taxid, rank, name = line.strip().split('\t')


            genbank_accession = gb_dict.get(taxid, None)
            if not genbank_accession:
                L.warn('No genbank accession number for TAXID: %s' % taxid)
                outf_discard.write('\t'.join([taxid, rank, name]) + '\n')
            else:
                outf.write('\t'.join([random.choice(genbank_accession),
                                      taxid,
                                      rank,
                                      name]) + '\n')
    outf_discard.close()
    

@follows(mkdir('taxon_info.dir/genbank_records'))
@split(fetchTaxonGenBankAccessions,'taxon_info.dir/genbank_records/*xml')
def fetchTaxonProjectIDs(infile, outfiles):
    '''For each taxon entry, use eutils to fetch the project accession
    information'''

    header = True
    for row in IOTools.openFile(infile):
        if header:
            header = False
            continue

        row = row.strip().split('\t')
        accession = row[0]
        taxid = row[1] 
        name = row[3]
        outfile = '_'.join(name.split()[:2]) + '_' + taxid + '.xml'
        outfile = os.path.join('taxon_info.dir/genbank_records', outfile)
        
        to_cluster = False
        statement = ("esearch -db nucleotide -query \"{accession}[ACCN]\" |"
                     " efetch -format gb "
                     " > {outfile}")
        P.run()
           

@merge(fetchTaxonProjectIDs, 'taxon_info.dir/taxon_project_ids.tsv')
def collateTaxonBioProjectNumbers(infiles, outfile):
    '''Parse the XML files, and output the project IDs'''

    outf = IOTools.openFile(outfile, 'w')
    outf.write('bioproject\ttaxon_id\ttaxon_name\n')
    outf_failed = IOTools.openFile(P.snip(outfile, '.tsv') + '_discarded.tsv', 'w')
    
    for infile in infiles:
        inf = os.path.basename(infile)[:-len('.xml')]
        tax_id = inf.split('_').pop()

        inf = IOTools.openFile(infile).readlines()
        if len(inf) == 0:
            L.warn('No genbank accession for taxon: %s' % infile)
            outf_failed.write(infile + '\n')
            continue

        name = None
        project = None
        for line in inf:
            line = line.lstrip().strip()
            if line.startswith('DBLINK'):
                assert line.split()[1] == 'BioProject:'
                project = line.split()[2]
            if line.startswith('ORGANISM'):
                name = ' '.join(line.split()[1:])

        if not project:
            L.warn('No genbank accession for taxon: %s' % infile)
            outf_failed.write(infile + '\n')
            continue

        outf.write('\t'.join([project, tax_id, name]) + '\n')

    outf.close()
    outf_failed.close()


@transform(collateTaxonBioProjectNumbers, regex('.+/(.+).tsv'), r'\1.load')
def loadTaxonSummary(infile, outfile):
    '''Load current assembly summary'''

    # fetch assembly summary table and edit
    df = pd.read_table(infile, index_col=0)
    
    # load... will throw error if table already exists
    table_name = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table_name, con=connect())
    
    open(outfile, 'w').close()

###############################################################################
## Section: Fetch NCBI RefSeq/Genbank Genome Assemblies For Taxa of Interest
###############################################################################
@originate('ncbi_assembly_summary.txt')
def fetchAssemblySummary(outfile):
    '''Fetch the NCBI assembly summary file... and remove UTF8 chars'''

    # wget assumes filename extension
    tmpf = os.path.basename(PARAMS['location_assembly_summary'])
    statement = ("wget {location_assembly_summary}")
    to_cluster = False
    P.run()
    
    statement = ("iconv -f UTF8 -t ascii//translit {tmpf} > {outfile}")
    to_cluster = False
    P.run()

    os.unlink(tmpf)


@transform(fetchAssemblySummary, suffix('.txt'), '.load')
def loadAssemblySummary(infile, outfile):
    '''Load current assembly summary'''
    
    # fetch assembly summary table and edit
    df = pd.read_table(infile, skiprows=1, index_col=0)
    df.index.name = 'assembly_accession'

    # load... will throw error if table already exists
    table_name = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table_name, con=connect(), if_exists='replace')

    open(outfile, 'w').close()


@follows(loadTaxonSummary)
@files(loadAssemblySummary,
       'ncbi_accessions.txt')
def fetchTaxonAccessions(infile, outfile):
    '''Find the accession details of all genomes belonging to taxon
    of interest'''
    
    infile = P.snip(infile[0], '.load')
    outf = IOTools.openFile(outfile, 'w')

    statement = ("SELECT"
                 "  a.taxon_id,a.taxon_name,a.bioproject,"
                 "  b.assembly_accession,b.taxid,b.species_taxid,b.ftp_path"
                 " FROM taxon_project_ids AS a"
                 " LEFT JOIN ncbi_assembly_summary AS b"
                 " WHERE a.bioproject == b.bioproject")

    df = pd.read_sql(sql=statement, con=connect(), index_col='taxon_id')

    assert df['taxid'].isnull().any() == False, 'Missing BioProjects in NCBI'

    # Fix ftp path...
    new_ftps = []
    for i in df['ftp_path']:
        genome_id = os.path.basename(i)
        genome_id = genome_id + '_genomic.fna.gz'
        new_ftp_path = i + '/' + genome_id
        new_ftps.append(new_ftp_path)

    df['ftp_path'] = new_ftps
    df.to_csv(outfile, sep='\t')
    
    outf.close()


@follows(mkdir('ncbi_assemblies.dir'))
@subdivide(fetchTaxonAccessions,
           regex('(.+).txt'),
           r'ncbi_assemblies.dir/\1_*.txt')
def splitTaxonAccessions(infile, outfiles):
    '''Split the download list into groups of 10 for download'''

    outf_stub = os.path.join('ncbi_assemblies.dir',
                             P.snip(os.path.basename(infile), '.txt'))
    outf = IOTools.openFile(outf_stub + '_00000.txt', 'w')
    for n, line in enumerate(IOTools.openFile(infile)):
        if n == 0:
            continue
        if n % 10 == 0:
            outf.close()
            outf = IOTools.openFile(outf_stub + '_' + str(n).zfill(5) + '.txt',
                                    'w') 
        outf.write(line)
    outf.close()


@subdivide(splitTaxonAccessions,
           regex('(.+)/(.+).txt'),
           r'\1/\2.dir/*gz')
def fetchGenomeAssemblyAndAnnotations(infile, outfiles):
    '''Download all assemblies using rsync, as per ncbi recommendation'''
    
    outdir = P.snip(infile, '.txt') + '.dir'
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    for genome in IOTools.openFile(infile):
        genome_path = genome.split().pop()
        genome_path = re.sub('ftp', 'rsync', genome_path)
        gff_path = P.snip(genome_path, '.fna.gz') + '.gff.gz'

        statement = (" rsync --copy-links --quiet"
                     "  %(genome_path)s"
                     "  %(outdir)s &&"
                     " rsync --copy-links --quiet"
                     "  %(gff_path)s"
                     "  %(outdir)s")
        to_cluster = False
        P.run()
    

@jobs_limit(1)
@subdivide(fetchTaxonAccessions,
           regex('(.+)_accessions.txt'),
           r'\1_16s_annotations.dir/*gz')
def fetchAssemblyRNAAnnotations(infile, outfiles):
    '''Download all annotation fasta using rsync, as per ncbi recommendation'''
    
    outdir = P.snip(infile, '_accessions.txt', strip_path=True) + '_16s_annotations.dir' 
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    header = True
    for genome in IOTools.openFile(infile):
        if header:
            header = False
            continue
        genome_path = genome.split().pop()
        genome_path = P.snip(genome_path, '_genomic.fna.gz') + '_rna_from_genomic.fna.gz'
        
        # not all entries have an '*_rna_from_genomic.fna.gz' file
        subdir = os.path.dirname(genome_path)[len('ftp://ftp.ncbi.nlm.nih.gov/'):]
        
        server = 'ftp.ncbi.nlm.nih.gov'
        user = 'anonymous'
        password = 'jethro.johnson@jax.org'

        ftp = ftplib.FTP(server)
        ftp.login(user,password)
        l = [os.path.basename(x) for x in ftp.nlst(subdir)]

        if os.path.basename(genome_path) in l:
            rsync_path = re.sub('ftp', 'rsync', genome_path)
            statement = ("rsync --copy-links --quiet"
                         " %(rsync_path)s"
                         " %(outdir)s")
            to_cluster = False
            P.run()


@jobs_limit(1)
@transform(fetchAssemblyRNAAnnotations,
           suffix('_rna_from_genomic.fna.gz'),
           '_16s_sequences.fasta.gz')
def extract16SSequences(infile, outfile):
    '''Extract any RNA fasta entry with the term 16S in the title'''

    outf = IOTools.openFile(outfile, 'w')

    for fasta in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        if re.search('16S', fasta.title, re.IGNORECASE):
            outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')

    outf.close()




###############################################################################
## Section: Fetch full taxonomy for downloaded assemblies
###############################################################################
@transform(fetchTaxonAccessions,
           suffix('_accessions.txt'),
           '_taxonomy.txt')
def fetchAssemblyTaxonomy(infile, outfile):
    '''Output the full taxonomy for the downloaded assemblies, WARNING:
    empty dictionaries are replaced with genus "unknown", species "unknown"'''

    output_dict = {}
    
    header = True
    for line in IOTools.openFile(infile):
        line = line.split('\t')
        if header:
            assert line[5] == 'species_taxid'
            assert line[3] == 'assembly_accession'
            header = False
            continue
        
        ncbi_taxon_id = line[5]

        # key = accession, value = full taxonomy
        output_dict[line[3]] = fetchBioSQLTaxonomy(ncbi_taxon_id)

    output_df = pd.DataFrame.from_dict(output_dict, orient='index')
    output_df.to_csv(IOTools.openFile(outfile, 'w'),
                     sep='\t',
                     index_label='accession_id')


@transform(fetchAssemblyTaxonomy, suffix('.txt'), '.load')
def loadAssemblyTaxonomy(infile, outfile):
    df = pd.read_table(infile)
    table = P.snip(outfile, '.load', strip_path=True)
    df.to_sql(name=table, con=connect(), if_exists='replace')
    open(outfile, 'w').close()




###############################################################################
## Section: Collate 16S gene sequences
###############################################################################
# convert the fasta headers to the genus_species_accession
@jobs_limit(1)
@transform(extract16SSequences,
           suffix('.fasta.gz'),
           add_inputs(loadAssemblyTaxonomy),
           '_renamed.fasta.gz')
def rename16SFastaSequences(infiles, outfile):
    '''Rename the 16S fasta sequences so headers are genus_species_accession'''

    fasta_file, taxonomy_file = infiles
    taxonomy_file = P.snip(taxonomy_file, '.load', strip_path=True)


    ##TO DO: Add the asm_name so it's accession_id_asm_name
    # fetch a dictionary of genus_species for each accession
    statement = ("SELECT accession_id,genus,species FROM %s" % taxonomy_file)
    cc = connect().execute(statement)

    tax_dict = {}
    for accession in cc.fetchall():
        accession_id, genus, species = accession
        assert accession_id not in tax_dict.keys(), 'Duplicate accession %s' \
            % accession_id

        # these samples have taxonomy information, but no genus level entry .
        if not genus:
            E.warn('Sample %s, which is species %s has no genus classification' \
                   % (accession_id, species))
            
        species = '_'.join(species.split())
        acc = accession_id.split('.')[0]
        
        tax_dict[accession_id] = '_'.join([species, acc])

    # fetch accession_id from file name
    accession_id = P.snip(os.path.basename(fasta_file), '_16s_sequences.fasta.gz')
    accession_id = '_'.join(accession_id.split('_')[:2])
    
    outfile = IOTools.openFile(outfile, 'w')
    n = 0
    for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
        n += 1
        out_header = tax_dict[accession_id] + '.' + str(n)
        outfile.write('>' + out_header + '\n' + fasta.sequence + '\n')


#@follows(mkdir('taxon_16s_genes.dir'))
@collate(rename16SFastaSequences,
         regex('(.+)_annotations.dir/.+_sequences_renamed.fasta.gz'),
         r'\1_genes.fasta.gz')
def collate16SFastaSequences(infiles, outfile):

    in_dir = 'ncbi_16s_annotations.dir'
    infiles = in_dir + '/' + '*_sequences_renamed.fasta.gz'

    statement = (" rm -f %(outfile)s &&"
                 " for FILE in %(infiles)s;"
                 "  do cat $FILE >> %(outfile)s;"
                 " done")
    
    to_cluster = False
    P.run()

    
@transform(collate16SFastaSequences,
           suffix('.fasta.gz'),
           '_stats.tsv.gz')
def fetch16SGeneStats(infile, outfile):

    statement = ("zcat %(infile)s |"
                 " python %(scriptsdir)s/fasta2table.py"
                 "  --log %(outfile)s.log"
                 "  --section=length,na,gaps |"
                 " gzip > %(outfile)s")
    to_cluster=False
    P.run()


@follows(fetch16SGeneStats)
@transform(collate16SFastaSequences,
           regex('(.+)_16s_genes.fasta.gz'),
           add_inputs(r'\1_16s_genes_stats.tsv.gz'),
           r'\1_16s_genes_filtered.fasta.gz')
def filter16SGenes(infiles, outfile):
    '''Drop anything with gaps and drop anything with sequence length
    below 1400bp'''

    fasta_file, stat_file = infiles

    # fetch a list of those entries to be filtered
    to_drop = set()
    df = pd.read_table(stat_file, index_col=0)
    to_drop.update(df[df['ngap_regions'] > 0].index.tolist())
    to_drop.update(df[df['length'] < 1400].index.tolist())

    # drop from fasta file
    with IOTools.openFile(outfile, 'w') as outf:
        for fasta in FastaIterator.FastaIterator(IOTools.openFile(fasta_file)):
            if fasta.title in to_drop:
                continue
            else:
                outf.write('\n'.join(['>' + fasta.title, fasta.sequence]) + '\n')   


###############################################################################
## Section: Align 16S genes
###############################################################################
@transform(filter16SGenes,
           suffix('_filtered.fasta.gz'),
           '_aligned.fasta.gz')
def align16SGenes(infile, outfile):
    '''Use muscle to create a global alignment of all 16S genes'''

    outfile = P.snip(outfile, '.gz')
    statement = ("muscle"
                 " -in <(gunzip -c %(infile)s)"
                 " -out %(outfile)s"
                 " -maxhours 24 &&"
                 " gzip -f %(outfile)s")

    to_cluster=False
    P.run()

###############################################################################
## Section: Build Tree
###############################################################################
@transform(align16SGenes,
           suffix("_aligned.fasta.gz"),
           "_tree.tre")
def createOTUAlignmentTree(infile, outfile):
    """
    Use fasttree to create phylogeny
    """

    tmpf = P.getTempFilename('.', suffix='.fasta')
    
    statement = ("zcat %(infile)s > %(tmpf)s &&"
                 " make_phylogeny.py"
                 "  -i %(tmpf)s"
                 "  -o %(outfile)s"
                 "  -r tree_method_default"
                 "  -t fasttree")
    to_cluster = False
    P.run()

    os.unlink(tmpf)





if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
