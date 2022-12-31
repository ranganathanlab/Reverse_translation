from Bio.Seq import Seq
from Bio import SeqIO
import os
import argparse

def get_seq(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    records_seq = [i.seq for i in records]
    headers = [i.description for i in records]
    return records_seq, headers

if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help='Input Sequence Alignment')
    parser.add_argument("outname", help='Output Sequence Alignment')
    options = parser.parse_args()
    
    seq, head = get_seq(options.alignment)
    seq = [str(i) for i in seq]
    
    savepath = 'Inputs/'
    output_name = savepath + options.outname + ".fasta"
    
    with open(output_name, 'w') as f:
        for i,item in enumerate(seq):
            f.write(">%s\n" %head[i])
            f.write("%s\n" %seq[i].replace('-',''))