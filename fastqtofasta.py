#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:00:44 2019

This program is designed to either convert a fastq file into a fasta file 
or to trim the fastq file

@author: breewenck
"""

def main():
    """
    This function is designed to determine what the user requires and then pass 
    the input to other functions for analysis. 
    """
    
    while True:
        file = input('If you want a fasta file please type "fasta" otherwise enter "fastq" if you want a trimmed fastq file: ')
        if file == 'fasta':
            in_file = input('Please enter a file: ')
            out_file = input('Please name your output file: ')
            print(fastq_to_fasta(in_file,out_file))    
        elif file == 'fastq':
            input_file = input('Please enter a file: ')
            output_file = input('Please name your output file: ')
            trim_5p = int(input("Please enter the number of nts to remove from the 5' end: "))
            trim_3p = int(input("Please enter the number of nts to remove from the 3' end: "))
            print(fastq_trimmer(input_file,output_file,trim_5p,trim_3p))
                
        return 'Function complete'
    

def fastq_to_fasta(in_file, out_file):
    """
    This function is designed to convert a fastq file to a fasta file.
    """
    try:
        input_handle = open(in_file)
        output_handle = open(out_file, 'w')
    except:
        return 'File not found'
    
    with input_handle, output_handle:
        ln = 0
        for line in input_handle:
            new_line = line.rstrip()
            ln += 1
            if (ln + 3) % 4 == 0:
                output_handle.write(f'>{new_line[1:]}\n')
            elif (ln + 3) % 4 == 1:
                output_handle.write(f'{new_line}\n')
            else:
                continue
        return ln/4  
    
             
def fastq_trimmer(input_file, output_file, trim_5p, trim_3p):
    """
    This function is designed to trim the 5' and 3' ends from the sequences and
    the quality score lines as designated by the user.
    """
    try:
        input_handle = open(input_file)
        output_handle = open(output_file, 'w')
    except:
        return 'File not found'
    
    with input_handle, output_handle:
        ln = 0
        for line in input_handle:
            new_line = line.rstrip()
            ln += 1
            if (ln + 3) % 4 == 1:
                output_handle.write(f'{new_line[trim_5p:-trim_3p]}\n')
            elif (ln + 3) % 4 == 3:
                output_handle.write(f'{new_line[trim_5p:-trim_3p]}\n')
            else:
                output_handle.write(f'{new_line}\n')
            
        return ln/4
 

if __name__ == '__main__':
    main()          