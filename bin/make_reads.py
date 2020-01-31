# make_reads.py
## Writen by Samuel Barnett
## This file contains the function for generating sequencing reads from DNA fragments to be used in conjunction with SIPSim_metagenome.py and nonSIP_metagenome.py.

from scipy.stats import norm as normdist
import random
from supplemental_functions import *


def GENreads(row, read_list_sub, insert_avg, insert_stddev, max_read_length):
    
    perc_incorp = float(row['percent_incorp'])
    BDi = float(row['OriBD'])
    length = int(row['fragment_length'])
    start = int(row['fragment_start'])
    taxon_name = str(row['taxon_name'])
    scaffoldID = str(row['scaffoldID'])
    library = int(row['library'])
    readcount = int(row['readcount'])

    # Call all functions that will be used so that you can save time from looking them up later
    lapp = read_list_sub.append
    RandNorm = random.gauss
    Rand = random.randint
    RandRand = random.random
    normdistcdf = normdist.cdf
    
    for read in range(1, readcount+1):
        cur_insert = int(RandNorm(insert_avg, insert_stddev))
        limit = length + start
        start1 = int(Rand(start, limit))
        check = RandRand()
        if check < 0.5:
            start2 = start1 + 2*max_read_length + cur_insert
            if start2 > limit:
                start2 = start1 - 2*max_read_length - cur_insert
        else:
            start2 = start1 - 2*max_read_length - cur_insert
            if start2 < start:
                start2 = start1 + 2*max_read_length + cur_insert
        lapp([taxon_name, scaffoldID, library, BDi, perc_incorp, start, length, 
                    read, start1, start2, max_read_length])
