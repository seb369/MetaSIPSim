#!/usr/bin/env python

# Function for generating simulated Metagenomic datasets without SIP
## Author: Samuel Barnett, 2019
## Read generation approximately follows the same pipeline as in SIPSim_metagenome. 
## Main difference is that fragment/read abudnances are derived solely from the species abundance table

import sys
import os
import numpy as np
import pandas as pd
import time
from datetime import datetime
import uuid
import ConfigParser
import multiprocessing
import csv
from Bio import SeqIO
import pyfasta
import StringIO
import argparse
import gzip
import pickle
from make_fragments import *
from make_reads import *
from write_fasta import *

###############################################################################################################################
############################################ Writing fragments into text file #################################################
###############################################################################################################################

### This section generates fragments and caluclates abundances then writes the tables of fragments to files. 
### Each file only contails a portion of the fragments and therefore will be concatinated together in the main function. 
### This further calls the function "generate_fragments" found in the file "make_fragments.py".

#### This function generates the inital fragments
def generate_fragments(genome_idx_chunk):
    workername = str(multiprocessing.Process())
    workername = str(uuid.uuid3(uuid.NAMESPACE_DNS, workername))
    temp_frag_outfile = tempdir+'/'+workername+'.initfragments.temp.txt'
    for cover in range(1, coverage_of_fragments+1):
        all_frag_entries = '\n'.join(genome_idx_chunk.apply(lambda row: get_frag_entries(row[0], os.path.join(genomeDir, row[1]), frag_length_distribution), axis=1))
        if cover == 1: 
            with open(temp_frag_outfile, 'w') as outfile:
                outfile.write(all_frag_entries)
                outfile.write('\n')
        else: 
            with open(temp_frag_outfile, 'a') as outfile:
                outfile.write(all_frag_entries)
                outfile.write('\n')

#### This function calculates fragment abundances from inital fragment map               
def write_fragments(frag_map_input):
    frag_map_input = frag_map_input[['taxon_name','scaffoldID','library','rel_abund_perc','fragStart','fragLength']].rename(index=str, columns={"rel_abund_perc":"abundance","fragStart":"fragment_start","fragStart":"fragment_length"})

    workername = str(multiprocessing.Process())
    workername = str(uuid.uuid3(uuid.NAMESPACE_DNS, workername))
    temp_frag_outfile = tempdir+'/'+workername+'.fragments.temp.txt.gz'
    frag_map_input.to_csv(temp_frag_outfile, sep='\t', index=False, header=False, compression='gzip')
    
###############################################################################################################################
###################################################### Generate/Select Reads ##################################################
###############################################################################################################################

### This section generates random reads from fragments generated earlier.
### This calls the function "GENreads" from file "make_reads.py".
            
def select_reads(frag_map_input):
    columns=['taxon_name', 'scaffoldID', 'library', 'OriBD', 'percent_incorp',
             'fragment_start', 'fragment_length',
             'read', 'forward_start', 'reverse_start', 'read_length']
    i = 0
    read_df_sub = []
    frag_map_input['OriBD'] = 0
    frag_map_input['percent_incorp'] = 0
    frag_map_input.apply(GENreads, axis=1, 
                         read_list_sub=read_df_sub, 
                         insert_avg=insert_avg, 
                         insert_stddev=insert_stddev,
                         max_read_length=max_read_length)
    read_df_sub = pd.DataFrame(read_df_sub, columns=columns).drop(columns=['OriBD', 'percent_incorp'])
    return read_df_sub

###############################################################################################################################
############################################ Writing read list into text file #################################################
###############################################################################################################################

### This section writes the list of reads out to text files. 
### All files only contain a subset of the reads and therefore are combined later on in the main function.
### This function is used if a list of reads is the desired output.

def write_readlist(readlist_chunk):
    workername = str(multiprocessing.Process())
    workername = str(uuid.uuid3(uuid.NAMESPACE_DNS, workername))
    temp_read_outfile = tempdir+'/'+workername+'.reads.temp.txt.gz'
    readlist_chunk.to_csv(temp_read_outfile, sep='\t', index=False, header=False, compression='gzip')
    
###############################################################################################################################
################################################ Writing sequence into fasta ##################################################
###############################################################################################################################

### This section pulls the sequences of the sequenced reads from the original genome files 
### and generates both forward and reverse fasta files containing read sequences.
### All files only contain a subset of the reads and therefore are concatinated later on in the main function.

def build_fasta_idxs(row, genomeDir, fasta_idxs):
    taxon_name = row['taxon_name']
    fasta_file = os.path.join(genomeDir, row['fasta_file'])
    taxon_fasta_idx = pyfasta.Fasta(fasta_file)
    fasta_idxs[taxon_name] = taxon_fasta_idx
    
def sequence_entry(read_map_subset):
    workername = str(multiprocessing.Process())
    workername = str(uuid.uuid3(uuid.NAMESPACE_DNS, workername))
    sub_F_fasta = tempdir+'/'+workername+'.forward.temp.fasta'
    sub_R_fasta = tempdir+'/'+workername+'.reverse.temp.fasta'
    coord = range(1, read_map_subset.shape[0]+1)
    read_map_subset['coordinates'] = coord
    forward_list = []
    reverse_list = []
    read_map_subset.apply(write_entry, axis=1, 
                          F_entry=forward_list, 
                          R_entry=reverse_list, 
                          fasta_idxs=fasta_idxs,
                          workuuid=workername)
    forward_entry = '\n'.join(forward_list)
    forward_list = None
    forward_entry = forward_entry + '\n'
    reverse_entry = '\n'.join(reverse_list)
    reverse_list = None
    reverse_entry = reverse_entry + '\n'
    with gzip.open(sub_F_fasta, 'wb') as forward:
        forward.write(forward_entry)
        forward_entry = None
    with gzip.open(sub_R_fasta, 'wb') as reverse:
        reverse.write(reverse_entry)
        reverse_entry = None

###############################################################################################################################
#################################################### Input variables ##########################################################
###############################################################################################################################

### Set start time for simulation to allow tracking of how long this takes
simulation_start_time = time.time()

### This section imports the many variables supplied in the configuration file. 
### For help generating the configuration file see the function (generate_config_file.py).
### For more info describing each variable see the manual.

parser = argparse.ArgumentParser(description='Simulating high throughput metagenome sequencing reads from a non-SIP experiment. Read generation approximately follows the same pipeline as in SIPSim_metagenome. The Main difference is that fragment/read abudnances are derived solely from the species abundance table, not from GC content and gradient position.')
parser.add_argument('cfg', type=str, help='The configuration file with all other variables used for this simulation.')
parser.add_argument('--version', '-v', action='version',
                    version='version 0.0.0')

args = parser.parse_args()

config_file = args.cfg

config = ConfigParser.ConfigParser()
config.read(config_file)

#### Temporary directory
tempdir = str(config.get('Other', 'temp_directory'))
if not os.path.exists(tempdir):
    os.makedirs(tempdir)

#### Threads for multiprocessing
threads = int(config.get('Other', 'threads'))

#### Log file
logfile = str(config.get('Other', 'logfile'))
outtext = 'Running nonSIP_metagenome\nThis program was writen by Samuel Barnett (seb369@cornell.edu)\n\nThis run was started on ' + datetime.now().strftime('%d/%m/%y at %H:%M:%S')
print outtext + '\n\n'
with open(logfile, 'w') as log:
    log.write(outtext + '\n\n\n')
    
#### Simulation endpoint ("read_sequence", "read_list", or "fragment_list")
endpoint = str(config.get('Other', 'endpoint'))
if endpoint == 'read_sequences':
    outtext = 'You have chosen to get sequences of simulated, sequenced SIP metagenome reads.'
elif endpoint == 'read_list':
    outtext = 'You have chosen to get an list of simulated, sequenced SIP metagenome reads.'
elif endpoint == 'fragment_list':
    outtext = 'You have chosen to get an list of simulated SIP metagenome fragments with adjusted buoyant densities.'
else:
    sys.exit('You must select "read_sequences", "read_list", or "fragment_list" in the configuration file.')
print outtext
with open(logfile, 'a') as log:
    log.write(outtext + '\n')
    
#### Libraries/gradients to be simulated (list of libraries separated by commas)
lib_string = str(config.get('Library', 'library_list'))
lib_list = [int(s) for s in lib_string.split(',')]

#### Fragment generation variables
genome_index_file = str(config.get('Fragment', 'genome_index_file'))
genomeDir = str(config.get('Fragment', 'genomeDir'))
frag_length_distribution = str(config.get('Fragment', 'frag_length_distribution'))
coverage_of_fragments = int(config.get('Fragment', 'coverage_of_fragments'))
temp_frag_file = str(config.get('Fragment', 'temp_fragment_file'))

#### Community composition file and incorporator stats
comm_file = str(config.get('Community', 'community_file'))

#### Sequencing read stats (only used if generating reads)
if endpoint == 'read_sequences' or endpoint == 'read_list':
    max_read_length = int(config.get('Sequencing', 'max_read_length'))
    insert_avg = float(config.get('Sequencing', 'avg_insert_size'))
    insert_stddev = float(config.get('Sequencing', 'stddev_insert_size'))
    finalN = int(config.get('Sequencing', 'final_number_of_sequences'))
    iterN = int(config.get('Sequencing', 'number_of_sequences_per_iteration'))
    num_iterations = finalN/iterN

## Iterations for fragments (only used if you are generating a list of fragments)
elif endpoint == 'fragment_list':
    num_iterations = int(config.get('Fragment', 'number_of_iterations'))
    
###############################################################################################################################
#################################################### Main function ############################################################
###############################################################################################################################

### This section includes the main function to generate metagenomic reads.
#### Build your sequence dictionary. 
#### This takes up a lot of memory and stays for all the iterations but it makes the process much faster. 
#### This is unnecessary unless you are actually making a fasta file in the end.
if endpoint == 'read_sequences':
    genome_index = pd.read_csv(genome_index_file, header=None, sep='\t')
    genome_index.columns = ['taxon_name', 'fasta_file']
    fasta_idxs = {}
    genome_index.apply(build_fasta_idxs, axis=1,
                       genomeDir = genomeDir,
                       fasta_idxs = fasta_idxs)
    del genome_index

#### Get pools started to allow for multiprocessing. 
#### This is done now as all the important variables used by later functions are already called.
pool = multiprocessing.Pool(processes=threads)

#### Generate fragments from input genome files
#### This will write these initial fragments to a temprary file that will be deleted at the end of the program.
#### Also, a separate fragment file will be made for each iteration if multiple iterations are used for memory saving.
outtext = 'Building fragments'
print '\n' + outtext
with open(logfile, 'a') as log:
    log.write('\n' + outtext + '\n')

frag_time = time.time()

if not os.path.isdir('./' + temp_frag_file):
    os.makedirs('./' + temp_frag_file)
else:
    outtext = 'Temporary fragment directory already exists. Your current files may be overwriten.'
    print '\n' + outtext + '\n'
    with open(logfile, 'a') as log:
        log.write('\n' + outtext + '\n')

genome_index = pd.read_csv(genome_index_file, header=None, sep='\t')

for n in range(1, num_iterations+1):
    iterfrag_file = temp_frag_file + '/frags_Iter' + str(n) + '.txt'
    if __name__ == '__main__':
        # calculate the chunk size as an integer
        chunk_size = int(genome_index.shape[0]/threads)
        genome_index = [genome_index.ix[genome_index.index[i:i + chunk_size]] for i in range(0, genome_index.shape[0], chunk_size)]
        pool.map(generate_fragments, genome_index)
        del genome_index
    else:
        outtext = 'Error with multiprocessing'
        print outtext
        with open(logfile, 'a') as log:
            log.write(outtext + '\n')
        sys.exit(outtext)
    
    # Now Concatinating all the intermediate fragment files together to get final file
    if n == 1:
        with open(iterfrag_file, 'wb') as fragoutfile:
            fragoutfile.write('taxon_name\tscaffoldID\tfragStart\tfragLength\tfragGC\n')
    cat_com = 'cat'
    rm_com = 'rm'
    for file in os.listdir(tempdir):
        if file.endswith('.initfragments.temp.txt'):
            cat_com = ' '.join([cat_com, os.path.join(tempdir, file)])
            rm_com = ' '.join([rm_com, os.path.join(tempdir, file)])
    cat_com = ' '.join([cat_com, '>>', iterfrag_file])
    os.system(cat_com)
    os.system(rm_com)
    del cat_com
    del rm_com

outtext = 'It took ' + str(round(time.time() - frag_time, 3)) + ' seconds to build the fragments\n\n----------'    
print '\n' + outtext + '\n'
with open(logfile, 'a') as log:
    log.write('\n' + outtext + '\n\n')

del frag_time

#### Run the simulation for each library/gradient used.
N = 0

for lib in lib_list:
    lib_time = time.time()    
    outtext = 'Starting library ' + str(lib)
    print outtext + '\n'
    with open(logfile, 'a') as log:
        log.write(outtext + '\n\n')
        
    start_time = time.time()
    
    # Generating output file names
    frag_outfile_name = 'nonSIP_library_' + str(lib) + '_fragments.txt.gz'
    outfile_test = frag_outfile_name

    if endpoint == 'read_sequences':
        fasta_prefix = 'nonSIP_library_' + str(lib) + '_reads'
        F_fasta = fasta_prefix + '_f.fasta.gz'
        R_fasta = fasta_prefix + '_r.fasta.gz'
        outfile_test = F_fasta
        
    elif endpoint == 'read_list':
        read_outfile_name = 'nonSIP_library_' + str(lib) + '_reads.txt.gz'
        outfile_test = read_outfile_name
        
    # If the output already exists we will skip it to avoid accidently overwriting the pre-existing file.
    if os.path.exists(outfile_test):
        print outfile_test + ' already exists. Please move, rename, or delete file before retrying this library!'
    else:
        for n in range(1, num_iterations+1):
            if num_iterations > 1:
                outtext = 'Starting iteration ' + str(n)
                print '\n' + outtext
                with open(logfile, 'a') as log:
                    log.write('\n' + outtext + '\n')
            newtime = time.time()

            iteration_fragfile = temp_frag_file + '/frags_Iter' + str(n) + '.txt'

            sub_df = pd.read_csv(iteration_fragfile, sep='\t')

            comm = pd.read_csv(comm_file, sep='\t')
            abund_table = comm.loc[comm['library'] == lib]
            comm = None

            sub_df['library'] = lib

            sub_df = sub_df.join(abund_table.set_index(['taxon_name','library']), on=['taxon_name','library'])
            sub_df = sub_df.fillna(0)
            abund_table = None
            incorp_table = None

            # Generating fragment files with adjusted fragment abundances based on BD.
            ## First generating the fragment file if first iteration and/or setting the dataframe headers.
            outtext = 'Writing fragments to file'
            print outtext
            with open(logfile, 'a') as log:
                log.write(outtext + '\n')
            header = ['taxon_name', 'scaffoldID', 'library', 'abundance', 'fragment_start', 'fragment_length']
            if n == 1:
                with gzip.open(frag_outfile_name, 'wb') as fragoutfile:
                    fragoutfile.write('\t'.join(header) + '\n')
            ## Adjusting abundances and writing fragment files.
            if __name__ == '__main__':
                # calculate the chunk size as an integer
                chunk_size = int(sub_df.shape[0]/threads)
                sub_df = [sub_df.ix[sub_df.index[i:i + chunk_size]] for i in range(0, sub_df.shape[0], chunk_size)]
                pool.map(write_fragments, sub_df)
                sub_df = None
            else:
                outtext = 'Error with multiprocessing'
                print outtext
                with open(logfile, 'a') as log:
                    log.write(outtext + '\n')
                sys.exit(outtext)
            ## Concatinating fragment files
            cat_com = 'cat'
            rm_com = 'rm'
            for file in os.listdir(tempdir):
                if file.endswith('.fragments.temp.txt.gz'):
                    cat_com = ' '.join([cat_com, os.path.join(tempdir, file)])
                    rm_com = ' '.join([rm_com, os.path.join(tempdir, file)])
            cat_com = ' '.join([cat_com, '>>', frag_outfile_name])
            os.system(cat_com)
            os.system(rm_com)
            del cat_com
            del rm_com
            
            outtext = 'It took ' + str(round(time.time() - newtime, 3)) + ' seconds to write fragment file.'
            print outtext
            with open(logfile, 'a') as log:
                log.write(outtext + '\n')
            newtime = time.time()
            
            # If you have chosen to get reads, then this continues
            if endpoint == 'read_list' or endpoint == 'read_sequences':
                # Number of reads for this iteration
                number_reads = iterN
                # Import fragment map
                sub_df = pd.read_csv(frag_outfile_name, sep='\t')
                # Get relative abundance of all potential reads from all fractions
                sub_df['abundance'] = (sub_df['abundance']*sub_df['fragment_length'])/(2*max_read_length)
                total_abundance = sum(sub_df['abundance'])
                sub_df['abundance'] = sub_df['abundance']/total_abundance
                # Selecting reads based on relative abundance estimate
                sub_df = sub_df.sample(n=number_reads, replace=True, weights='abundance', axis=0)
                sub_df = sub_df.groupby(sub_df.columns.tolist()).size().reset_index().rename(columns={0:'readcount'})
                # Variables for final readmap
                columns=['taxon_name', 'scaffoldID', 'library', 'abundance', 'fragment_start', 'fragment_length',
                         'read', 'forward_start', 'reverse_start', 'read_length']
                start_time = time.time()

                # Generating reads from the fragments
                if __name__ == '__main__':
                    # calculate the chunk size as an integer
                    chunk_size = int(sub_df.shape[0]/threads)
                    sub_df = [sub_df.ix[sub_df.index[i:i + chunk_size]] for i in range(0, sub_df.shape[0], chunk_size)]
                    sub_df = pool.map(select_reads, sub_df)

                else:
                    outtext = 'Error with selecting reads multiprocessing'
                    print outtext
                    with open(logfile, 'a') as log:
                        log.write(outtext + '\n')
                    sys.exit(outtext)

                outtext = 'It took ' + str(round(time.time() - newtime, 3)) + ' seconds to generate reads. Now building map.'
                print outtext
                with open(logfile, 'a') as log:
                    log.write(outtext + '\n')
                newtime = time.time()
                
                # If you have selected to get just a list of reads, writing them to a file here.
                if endpoint == 'read_list':
                    outtext = 'Writing read list to file'
                    print outtext
                    with open(logfile, 'a') as log:
                        log.write(outtext + '\n')
                    columns=['taxon_name', 'scaffoldID', 'library','fragment_start', 'fragment_length',
                             'read', 'forward_start', 'reverse_start', 'read_length']
                    header = '\t'.join(columns)
                    if n == 1:
                        with gzip.open(read_outfile_name, 'wb') as readoutfile:
                            readoutfile.write(header + '\n')
                    if __name__ == '__main__':
                        pool.map(write_readlist, sub_df)
                        sub_df = None
                    else:
                        outtext = 'Error with multiprocessing'
                        print outtext
                        with open(logfile, 'a') as log:
                            log.write(outtext + '\n')
                        sys.exit(outtext)
                    # Combining all the partial read lists.
                    cat_com = 'cat'
                    rm_com = 'rm'
                    for file in os.listdir(tempdir):
                        if file.endswith('.reads.temp.txt.gz'):
                            cat_com = ' '.join([cat_com, os.path.join(tempdir, file)])
                            rm_com = ' '.join([rm_com, os.path.join(tempdir, file)])
                    cat_com = ' '.join([cat_com, '>>', read_outfile_name])
                    os.system(cat_com)
                    os.system(rm_com)
                    del cat_com
                    del rm_com

                    outtext = 'It took ' + str(round(time.time() - newtime, 3)) + ' seconds to write read list file'
                    print outtext
                    with open(logfile, 'a') as log:
                        log.write(outtext + '\n')
                    newtime = time.time()
                    
                # Else, if you have selected to get the read sequences, generating the fasta files here.
                else:
                    outtext = 'Writing read sequences'
                    print outtext
                    with open(logfile, 'a') as log:
                        log.write(outtext + '\n')
                    if __name__ == '__main__':
                        pool.map(sequence_entry, sub_df)
                        sub_df = None
                    else:
                        outtext = 'Error with multiprocessing'
                        print outtext
                        with open(logfile, 'a') as log:
                            log.write(outtext + '\n')
                        sys.exit(outtext)
                    # Combining all the partial fasta files.
                    cat_com = 'cat'
                    rm_com = 'rm'
                    for file in os.listdir(tempdir):
                        if file.endswith('.forward.temp.fasta'):
                            cat_com = ' '.join([cat_com, os.path.join(tempdir, file)])
                            rm_com = ' '.join([rm_com, os.path.join(tempdir, file)])
                    cat_com = ' '.join([cat_com, '>>', F_fasta])
                    os.system(cat_com)
                    os.system(rm_com)
                    cat_com = 'cat'
                    rm_com = 'rm'
                    for file in os.listdir(tempdir):
                        if file.endswith('.reverse.temp.fasta'):
                            cat_com = ' '.join([cat_com, os.path.join(tempdir, file)])
                            rm_com = ' '.join([rm_com, os.path.join(tempdir, file)])
                    cat_com = ' '.join([cat_com, '>>', R_fasta])
                    os.system(cat_com)
                    os.system(rm_com)
                    del cat_com
                    del rm_com

                    outtext = 'It took ' + str(round(time.time() - newtime, 3)) + ' seconds to make this fasta'
                    print outtext
                    with open(logfile, 'a') as log:
                        log.write(outtext + '\n')
                        
        outtext = 'It took ' + str(round(time.time() - start_time, 3)) + ' seconds to run library ' + str(lib) + ' iteration ' + str(n)
        print outtext + '\n'
        with open(logfile, 'a') as log:
            log.write(outtext + '\n\n')
                            
        
    outtext = 'It took ' + str(round(time.time() - lib_time, 3)) + ' seconds to run the whole library ' + str(lib) + '\n\n----------'
    print outtext + '\n'
    with open(logfile, 'a') as log:
        log.write(outtext + '\n\n')
        
# Cleanup temporary fragment files for each iteration
rm_com = 'rm'
for file in os.listdir(temp_frag_file):
    if file.startswith('frags_Iter'):
        rm_com = ' '.join([rm_com, os.path.join(temp_frag_file, file)])
os.system(rm_com)
del rm_com

outtext = 'It took ' + str(round(time.time() - simulation_start_time, 3)) + ' seconds to run the entire simulation. This program is now done!'
print outtext
with open(logfile, 'a') as log:
    log.write(outtext + '\n\n\n')
        
pool.close()
pool.terminate()
pool.join()