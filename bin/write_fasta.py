# write_fasta.py
## Writen by Samuel Barnett
## This file contains functions for finding the reverse compliment of a DNA sequence and for writing a fasta formated string in conjunction with the SIPSim_metagenome program.

def reverse_compliment(sequence):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 
            'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 
            'B':'V', 'V':'B', 'D':'H', 'H':'D', 
            'N':'N', '.':'.', '-':'-'}
    # Note sequences should be a string object
    new_sequence = str()
    sequence = sequence[::-1]
    for base in sequence:
        new_sequence += comp[base]
    sequence = None
    return new_sequence

def write_entry(row, F_entry, R_entry, fasta_idxs, workuuid):
    Fapp = F_entry.append
    Rapp = R_entry.append
    
    taxon_name = row['taxon_name']
    scaffold = row['scaffoldID']
    library = str(row['library'])
    read_length = int(row['read_length'])
    forward_start = int(row['forward_start'])
    reverse_start = int(row['reverse_start'])
    reverse_start = int(row['reverse_start'])
    coord = str(row['coordinates'])
    
    read = str(fasta_idxs[taxon_name][scaffold][min(forward_start,reverse_start):max(forward_start,reverse_start)])
    if forward_start < reverse_start:
        forward_seq = read[0:read_length]
        reverse_seq = reverse_compliment(read[len(read)-read_length:len(read)])
    elif forward_start > reverse_start:
        reverse_seq = read[0:read_length]
        forward_seq = reverse_compliment(read[len(read)-read_length:len(read)])
    read = None
    Fapp('>' + workuuid.replace('-','_') + ':1:NCC1701:1:1:' + coord + ':' + coord + ' 1:N:0:' + library + '\n' + forward_seq)
    Rapp('>' + workuuid.replace('-','_') + ':1:NCC1701:1:1:' + coord + ':' + coord + ' 2:N:0:' + library + '\n' + reverse_seq)
    #Fapp('>'+scaffold+'____'+OriBD+'____'+percent_incorp+'____'+fragment_start+'____'+fragment_length+'\n'+forward_seq)
    #Rapp('>'+scaffold+'____'+OriBD+'____'+percent_incorp+'____'+fragment_start+'____'+fragment_length+'\n'+reverse_seq)