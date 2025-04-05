#! /usr/bin/env python3
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--padding", '-p', type=int, 
        help="How much padding was added to the beginning and end of \
all seqs (in bases)?", default=1000)
    parser.add_argument("--alignment", "-a", 
        help="Alignment in FASTA format", required=True)
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    seqs = {}
    seqid = None
    seq = ""
    f = open(options.alignment, 'r')
    for line in f:
        line = line.rstrip()
        if line[0] == ">":
            if seqid is not None:
                seqs[seqid] = seq
            seq = ""
            seqid = line[1:].split(' ')[0]
        else:
            seq += line.upper()
    if seqid is not None:
        seqs[seqid] = seq

    seq1idx = 0
    seq2idx = 0
    seq3idx = 0
        
    consensus = []
    
    s1len = 0
    for base in seqs['seq1']:
        if base != '-':
            s1len += 1
    s1len -= options.padding

    for idx in range(0, len(seqs['seq1'])):
        seq1base = seqs['seq1'][idx]
        seq2base = seqs['seq2'][idx]
        seq3base = seqs['seq3'][idx]
        
        if seq1idx > s1len-1:
            break
        else:
            if seq1idx >= options.padding:
                consbase = '-'
                if seq1base == seq2base:
                    consbase = seq1base
                else:
                    consbase = seq3base
                if consbase != '-':
                    consensus.append(consbase)

            if seq1base != '-':
                seq1idx += 1
            if seq2base != '-':
                seq2idx += 1
            if seq3base != '-':
                seq3idx += 1
    
    print(">chrM")
    print("".join(consensus))


if __name__ == '__main__':
    sys.exit(main(sys.argv))
