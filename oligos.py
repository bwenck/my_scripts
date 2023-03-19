#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:25:02 2019

This is a program designed to analyze a sequence to assist in designing 
oligonucleotides.

@author: breewenck
"""

import string


def main(): 
    while True:
        rna = input('Please enter a sequence: ')
        if rna == '':
            break
        else:
            print('Your cleaned up sequence is:', clean_up(rna))
            print('The reverse complement of your sequence is:', rev_comp(clean_up(rna)))
            print('Length:', len(clean_up(rna)), 'nts')
            print('GC content:', gc_content(clean_up(rna)), '%')
            print('Tm:', tm(clean_up(rna)), '℃')
"""
This function is designed to request the user for a sequence and then pass the 
input to other functions to analyze the sequence.

"""          
def clean_up(seq):
    non_DNA = ''
    new_seq = ''
    nts = seq
    num_trans = str.maketrans('ATGCU', '12345') 
    num = seq.translate(num_trans) 
    for nt in seq:
        if nt == 'U':
            non_DNA += 'T' 
        elif nt == 'A':
            non_DNA += 'A'
        elif nt == 'C':
            non_DNA += 'C'
        elif nt == 'G':
            non_DNA += 'G'
        elif nt in string.whitespace:
            continue
            new_seq += nt
        else:
            if nts != num:
                break
            
    return non_DNA
"""
This function is designed to analyze an inputed sequence and clean the sequence
to result in DNA nucleotides.

"""
def rev_comp(seq):
    trans_table = str.maketrans('ATGC','TACG')
    rc = seq.translate(trans_table)
           
    return rc
"""
This function is designed to return the reverse complement of the cleaned up 
inputed sequence.

"""  
def gc_content(seq):
    count = 0
    total = len(seq)
    gc_total = seq.count('G') + seq.count('C')
    gc_content = float(gc_total) / total * 100 
    for i in seq:
        if i == 'C' or i == 'G':
            count += 1
    
    return round(gc_content, 2)
"""
This function is designed to determine the G/C content of the inputed sequence.

"""

def tm(seq):
    count = 0
    a = seq.count('A')
    t = seq.count('T')
    c = seq.count('C')
    g = seq.count('G')
    at = [a, t]
    gc = [c, g]
    tm = (2*sum(at)) + (4*sum(gc))
    for i in seq:
        if i == 'C' and i == 'G' and i == 'A' and i == 'T':
            count += 1
            
    return tm
"""
This function is designed to return the melting temperature of the cleaned up 
inputed sequence.

"""
if __name__ == '__main__':
    main()
    
    
