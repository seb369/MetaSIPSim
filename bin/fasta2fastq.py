#!/usr/bin/env python

# Function for converting fasta files into fastq files
## Here we take the fasta files generated with MetaSIPSim and convert them into fastq files
## complete with sequencing errors and quality scores. Error and quality based on error models.
## Error models as well as some code and functions are taken from InSilicoSeq:
### Hadrien Gourlé, Oskar Karlsson-Lindsjö, Juliette Hayer, Erik Bongcam-Rudloff (2019)
### Simulating Illumina metagenomic data with InSilicoSeq. Bioinformatics, 35,3:521–522,
### https://doi.org/10.1093/bioinformatics/bty630

import os
import sys
import pysam
import numpy as np
import argparse
import time
import multiprocessing
import gzip
from scipy import stats
from random import random
from Bio import SeqIO


from iss import util
from iss import modeller
from iss.error_models import kde

import argparse
parser = argparse.ArgumentParser(description='Converting a fasta file into a fastq file. Sequence qualities and errors are based on the supplied error model.')
parser.add_argument("i", type=str, help="The input fasta file")
parser.add_argument("d", type=str, help="Is this a [forward] or [reverse] read set")
parser.add_argument("e", type=str, help="Error model file. For more info on this, see the InSilicoSeq documentation.")
parser.add_argument("l", type=int, help="Read length")
parser.add_argument("t", type=str, help="Path for the temporary directory")
parser.add_argument("p", type=int, help="number of processors")
parser.add_argument('--version', '-v', action='version',
                    version='version 0.0.0')

args = parser.parse_args()

## Function for breaking up input fasta into batches
def batch_iterator(iterator, batch_size):
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

## Function for error modification and writing fastq
def error_modification(fasta, ErrorModel, direction):
    if len(fasta.seq) >= args.l:
        if len(fasta.seq) > args.l:
            fasta.seq = fasta.seq[0:args.l]
        fasta = ErrorModel.introduce_error_scores(fasta, direction)
        fasta.seq = ErrorModel.mut_sequence(fasta, direction)
        mutseq = fasta.seq.tomutable()
        pos = 0
        for nucl in mutseq:
            if nucl in ['R','Y','K','M','S','W','B','D','H','V']:
                mutseq[pos] = 'N'
            pos += 1
        fasta.seq = mutseq.toseq()
        return fasta.format('fastq')
    else:
        print "Skipping sequence " +  fasta
        return "ERROR_len0"
    
def genfasta_MP(fastafile):
    fastq_out = fastafile.replace('fasta', "fastq")
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    fastqentry = str()
    for fasta in fasta_sequences:
        entry = error_modification(fasta, ErrorModel, args.d)
        if entry != 'ERROR_len0':
            fastqentry+=entry
    print "writing fastq file: " + fastq_out
    with gzip.open(fastq_out, 'wb') as output_fq:
        output_fq.write(fastqentry)

# Main function
start_time = time.time()
       
fasta_sequences = SeqIO.parse(open(args.i),'fasta')
ErrorModel = kde.KDErrorModel(args.e)

if not os.path.isdir(args.t):
    os.makedirs(args.t)

chunk_size = int(10000000/(args.p*2))
print "Chunk size: " + str(chunk_size)

for i, batch in enumerate(batch_iterator(fasta_sequences, chunk_size)):
    filename = os.path.join(args.t, ("subfasta_%i.fasta" % (i + 1)))
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
    
files = [os.path.join(args.t, f) for f in os.listdir(args.t) if 'subfasta_' in f] 

pool = multiprocessing.Pool(processes=args.p)
pool.map(genfasta_MP, files)

fastq_out = args.i.replace('fasta', "fastq.gz")
print("writing to " + fastq_out)

files = [os.path.join(args.t, f) for f in os.listdir(args.t) if 'subfastq_' in f] 
catcomm = 'cat '
catcomm = catcomm + ' '.join(files) + ' > ' + fastq_out
os.system(catcomm)

files = [os.path.join(args.t, f) for f in os.listdir(args.t) if 'subfast' in f] 
rmcomm = 'rm '
rmcomm = rmcomm + ' '.join(files)
os.system(rmcomm)

outtext = 'It took ' + str(round(time.time() - start_time, 3)) + ' seconds to convert from fasta to fastq. This program is now done!'

print outtext