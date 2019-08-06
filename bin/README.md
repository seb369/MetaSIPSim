# Brief explination of each file in MetaSIPSim/bin/

## Executable scripts
SIPSim_metagenome.py: This file contains the main script for simulating metagenomic-SIP datasets. 
    To run this script you need a configuration file. For more information on running this script 
    please refer to the main page or the examples. The command to run this script is:
        python SIPSim_metagenome.py [Config file]

nonSIP_metagenome.py: This file contains the main script for simulating conventional shotgun metagenome 
    datasets. To run this script you need a configuration file which can be the same one used for 
    simulating the metagenomic-SIP dataset. For more information on running this script please refer to 
    the main page or the examples. The command to run this script is:
        python nonSIP_metagenome.py [Config file]
        
fasta2fastq.py: This file contains the script for converting a fasta file, like the output from 
    SIPSim_metagenome.py and nonSIP_metagenome.py, into a fastq file format like you would get from a 
    standard Illumina sequencing run. Sequencing error and quality uses an error model developed
    from an actual sequencing run. This script uses functions derived or directly from InSilicoSeq
    (https://doi.org/10.1093/bioinformatics/bty630). Error models must be generated before hand. 
    See the manual for InSilicoSeq to do this or simply use those provided by that software.

## Functions used by the scripts
make_fragments.py: This file contains functions for adjusting a fragments abundance based on its position
    in a buoyant density gradient. These functions are called by SIPSim_metagenome.py.

make_reads.py: This file contains functions for generating simulated next-generation sequencing reads
    from simulated DNA fragments. These functions are called by both 
    SIPSim_metagenome.py and nonSIP_metagenome.py.

write_fasta.py: This file contains functions for writing simulated next-generation sequencing reads
    into fasta format. These functions are called by both SIPSim_metagenome.py and nonSIP_metagenome.py.

supplimental_functions.py: This file contains functions used by both SIPSim_metagenome.py and 
    make_fragments.py for calculating fragment buoyant densities and building gradient models. Most 
    functions here are derived from SIPSim (https://github.com/nick-youngblut/SIPSim).
