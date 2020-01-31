# make_fragments.py
## Writen by Samuel Barnett
## This file contains the function for finding BD of DNA fragments to be used in conjunction with SIPSim_metagenome.

from scipy.stats import norm as normdist
from scipy.stats import skewnorm as skewnormdist
from scipy.stats import uniform as uniformdist
from scipy.stats import truncnorm as truncnormdist
import random
from Bio import SeqIO
from supplemental_functions import *

## Functions for generating fragments from input genome files
def get_frag_GC(start, end, scaffold):
    frag_seq = scaffold.seq[start:end]
    GC_count = frag_seq.count('G') + frag_seq.count('C') + frag_seq.count('c') + frag_seq.count('g')
    N_count = frag_seq.count('N') + frag_seq.count('n')
    fragGC = float(GC_count)/(len(frag_seq)-N_count)*100
    return(fragGC)

def rand_frag_size(dist_string):
    dist_list = dist_string.split(',')
    if dist_list[0] == 'skewed-normal':
        rand_size = abs(int(round(skewnormdist.rvs(float(dist_list[3]), loc=float(dist_list[1]), scale=float(dist_list[2])))))
    elif dist_list[0] == 'normal':
        rand_size = abs(int(round(normdist.rvs(loc=float(dist_list[1]), scale=float(dist_list[2])))))
    elif dist_list[0] == 'uniform':
        rand_size = abs(int(round(uniformdist.rvs(loc=float(dist_list[1]), scale=float(dist_list[2])))))
    elif dist_list[0] == 'truncated-normal':
        rand_size = abs(int(round(truncnormdist.rvs(float(dist_list[3]), float(dist_list[4]), loc=float(dist_list[1]), scale=float(dist_list[2])))))
    return(rand_size)

def get_scaff_entries(genome_name, scaffold, dist_string):
    scaffold_length = len(scaffold.seq)
    frag_len_list = list()
    frag_sum = 0
    while frag_sum < scaffold_length:    
        frag_len_list.append(rand_frag_size(dist_string))
        frag_sum = sum(frag_len_list)

    if frag_sum > scaffold_length:
        del frag_len_list[-1]
        if len(frag_len_list) > 0:
            new_frag = scaffold_length - sum(frag_len_list)
            if new_frag < min(frag_len_list) :
                new_frag = min(frag_len_list) + scaffold_length - sum(frag_len_list)
                frag_len_list.remove(min(frag_len_list))                
        else:
            new_frag = scaffold_length
        frag_len_list.append(new_frag)
    random.shuffle(frag_len_list)

    df = pd.DataFrame({"init_lengths": frag_len_list})
    df = df.cumsum(axis=0)
    df['fragStart'] = df['init_lengths'] - frag_len_list
    df['fragEnd'] = df['init_lengths'] - 1
    df['fragLength'] = df['fragEnd'] - df['fragStart']
    df['fragGC'] = df.apply(lambda row: get_frag_GC(row.fragStart, row.fragEnd, scaffold), axis=1)
    df['scaffoldID'] = scaffold.description
    df['taxon_name'] = genome_name
    scaff_entries = '\n'.join(df.apply(lambda row: '\t'.join([row.taxon_name, row.scaffoldID, 
                                                              str(row.fragStart), str(row.fragLength),
                                                              str(row.fragGC)]), axis=1))
    return(scaff_entries)

def get_frag_entries(genome_name, fasta_path, dist_string):
    fasta = SeqIO.parse(open(fasta_path), 'fasta')
    scaff_entry_list = [get_scaff_entries(genome_name, scaffold, dist_string) for scaffold in fasta]
    return('\n'.join(scaff_entry_list))

## Function for caclulating abundance of fragments
def GENfragments(row, frag_list_sub, BD_models, BD_shift, diffVar, D, B, w, I, r_max, A, r, minBd, maxBd, pDBL):
    
    percent_incorp = float(row['percent_incorporation'])
    sd_incorp = float(row['sd_incorporation'])
    BDi = float(row['fragBD'])
    length = int(row['fragLength'])
    start = int(row['fragStart'])
    taxon_name = str(row['taxon_name'])
    scaffoldID = str(row['scaffoldID'])
    library = int(row['library'])
    genome_abundance = float(row['rel_abund_perc'])
    
    # BD models
    BDtubefit = BD_models[library]

    # Call all functions that will be used so that you can save time from looking them up later
    lapp = frag_list_sub.append
    RandNorm = random.gauss
    Rand = random.randint
    RandRand = random.random
    normdistcdf = normdist.cdf

    # Next add incorporation status, only if in treatment library
    if percent_incorp != 0:
        perc_incorp = RandNorm(percent_incorp, sd_incorp)/100
        if perc_incorp > 1:
            perc_incorp = 1
        if perc_incorp < 0:
            perc_incorp = 0
        BD = BDi + (perc_incorp * BD_shift)
    elif percent_incorp == 0:
        BD = BDi
        perc_incorp = 0
        
    # Now find standard Deviation for diffusion
    sig = diffusion_SD(BD, length, diffVar)
        
    # Now find max and min of DBL for fragment
    x = BD2distFromAxis(BD, D, B, w, I)
    angle_tube_pos = axisDist2angledTubePos(x, r_max, A, r)
    DBL_range = angle_tube_pos_to_BD(angle_tube_pos, BDtubefit)
        
    # Now find the abundance of the fragment normally in this BD window
    window_perc = normdistcdf(maxBd, loc=BD, scale=sig) - normdistcdf(minBd, loc=BD, scale=sig)
    lumen_abundance = (1 - pDBL) * genome_abundance * window_perc
             
    # Now find the abundance of the fragment in the DBL of this BD window
    dbl_abundance = pDBL * genome_abundance * DBL_prop(DBL_range, maxBd, minBd)
        
    # Now add togethere the lumen_abundance and the DBL_abundance for total window abundance
    abundance = lumen_abundance + dbl_abundance
    
    lapp([taxon_name, scaffoldID, library, BDi, perc_incorp, abundance, start, length])
    