#! /home/johnsj/devel/labscripts/analysis/ruffus/ruffus2.7/bin/python

import os
import argparse
import logging as L

# pass the input and output directories as positional args
parser = argparse.ArgumentParser()
parser.add_argument('indir', type=str,
                    help='directory of the files to be analysed')
parser.add_argument('outdir', type=str,
                    help='directory in which files are to be linked')
args = parser.parse_args()


# create out_directory
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)

# create symlinks for all the file in the input directory
sample_ids = []
out_dir = os.path.abspath(args.outdir)
for inf in os.listdir(args.indir):
    inf = os.path.join(args.indir, inf)
    if not os.path.isdir(inf):
        print 'file %s is not a directory' % inf
        continue

    # get the JAX ID for the sample
    jax_id = inf.split('.').pop()
    assert jax_id not in sample_ids, 'Multiple samples with the same ID'

    src = os.path.join(inf, 'ccs_filter_trim_clean.fasta')
    link_name = os.path.join(os.path.abspath(args.outdir), jax_id + '.fasta')
    os.symlink(src, link_name)
    L.info('Created symlink for %s' % link_name)
