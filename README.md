# MetaSIPSim

Welcome to the project page for MetaSIPSim. This is a tool for simulating metagenomic stable isotope probing (metagenomic-SIP) datasets, including next-generation read libraries. Metagenomic-SIP allows us to link microbial genomes to function within an environmental sample. Metagenomic-SIP improves our ability to recover, assemble, and bin genomes from target, isotopically labeled microbes, compared to conventional shotgun metagenomics.

MetaSIPSim can be used to simulate read libraries for use in testing experimental parameters and in the development of analytical tools for metagenomic-SIP methodologies. This can save on costs and time where preliminary experiments or *in vitro* mock communities would otherwise be used.

## Requirements 

MetaSIPSim can be run using linux or unix based operating systems and has been successfully and fully tested using an Ubuntu version 16.04.4 system as well as macOSX version 10.12.6. MetaSIPSim is currently incompatible with windows.

MetaSIPSim runs with python 2.7 and is currently incompatible with python 3.6. For those familiar with it, we highly recommend using anaconda environments for version control and installing the dependencies. 

## Installation

### Dependencies

To run MetaSIPSim scripts you must install the following python modules. All are available through conda.

* numpy (≥ 1.16.3)
* pandas (≥ 0.24.2)
* Bioptyhon (≥ 1.73)
* scipy (≥ 1.2.1)
* pyfasta (≥ 0.5.2)

If you want to use the included script for converting fasta files to fastq format, you also need to install 

* InSilicoSeq (≥ 1.3.6)

### Installing MetaSIPSim

To install MetaSIPSim, simply clone this repo. You can then call the MetaSIPSim scripts directly from the cloned directory.

## Scripts

* **SIPSim_metagenome.py:** This script simulates metagenomic-SIP datasets based on input experimental parameters and a simulated community. Output from this script can be a table of DNA fragments with their abundances in a given gradient window, a table of next-gen sequencing reads recovered in a gradient window, or these recovered reads in fasta format.

* **nonSIP_metagenome.py:** This script simulates conventional shotgun metagenomic datasets. This can be used to compare metagenomic-SIP methodologies to conventional metagenomics. Input can the identical to that of SIPSim_metagenome.py. Output from this script can be a table of DNA fragments with their abundances in a given gradient window, a table of next-gen sequencing reads recovered in a gradient window, or these recovered reads in fasta format.

* **fasta2fastq.py:** This script converts fasta formatted read libraries generated by SIPSim_metagenome.py and nonSIP_metagenome.py into fastq format with sequencing errors and quality scores. This script relies heavily on InSilicoSeq functions.

## Simulating datasets

### Metagenome Simulations
For both SIPSim_metagenome.py and nonSIP_metagenome.py, input is a configuration file containing paths to input files or directories, experimental parameters including gradient variables, and sequencing variables. For a more thorough description of the configuration file, please see MetaSIPSim_config.pdf. This same configuration file can be used for both simulations however the gradient parameters are not used for nonSIP_metagenome.py

Once the configuration file is generated, you can run these scripts by running:

`python PATH_TO_MetaSIPSim/bin/SIPSim_metagenome.py CONFIG_FILE`

`python PATH_TO_MetaSIPSim/bin/nonSIP_metagenome.py CONFIG_FILE`

### Converting fasta to fastq

To convert the simulated read library in fasta to fastq format, simply run:

`python PATH_TO_MetaSIPSim/bin/fasta2fastq.py i d e l t p`

This script requires the following inputs:
* **i:** The input fasta file
* **d:** Indicates whether the input file contains [forward] or [reverse] reads
* **e:** The error model file, which can be generated with InSilicoSeq or one of the ones provided
* **l:** The length of a read in bp
* **t:** The path to a temporary directory to store intermediate files in
* **p:** The number of processors to use

## Tutorials

example_simulation.ipynb contains a simple example for simulating metagenomic-SIP libraries as well as conventional shotgun libraries using the same parameters. For more in depth examples you can also check out the validation notebooks or the case study scripts.

## Citation
TBD

## Versions

* **0.1.0:** Initial commit of development version.

## Acknowledgements:

Many thanks to Nicholas Youngblut, the author of the original SIPSim toolkit, for advice and mentorship during development of this tool. Some of the functions, particularly those involved in dealing with fragment buoyant densities, are derived or taken from SIPSim.

For more information on SIPSim, please refer to:

Youngblut ND, Barnett SE, Buckley DH. SIPSim: A Modeling Toolkit to Predict Accuracy and Aid Design of DNA-SIP Experiments. Frontiers in Microbiology 2018;9:570. doi: [10.3389/fmicb.2018.00570](https://doi.org/10.3389/fmicb.2018.00570)

or the project site https://github.com/nick-youngblut/SIPSim.

The script for converting read library files from fasta to fastq format relies heavily on code from InSilicoSeq. For more information on this program please refer to:

Gourlé H, Karlsson-Lindsjö O, Hayer J, Bongcam-Rudloff E. Simulating Illumina data with InSilicoSeq. Bioinformatics 2018;35(3):521-522. doi: [10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)

## License
TBD

