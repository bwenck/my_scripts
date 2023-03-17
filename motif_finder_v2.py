#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 09:41:11 2021

@author: breewenck
"""

import argparse
import re

def main():
    """
    Call arg_parse, motif_finder, and line_finder functions
    """
    input_file, motif = arg_parse()
    print(new_file(input_file))
    print(motif_finder(input_file, motif))
    print(line_finder(input_file, motif))

    
def arg_parse():
    """
    This function is designed to add input_file and motif arguments and parse 
    to the function from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True, help="fasta file") 
    parser.add_argument('-m', '--motif', required=True, help="desired motif; N is a valid nt")
    
    args = parser.parse_args()
    
    return args.input_file, args.motif 
    
def new_file(input_file):
	"""
    This function is designed to create a new file to be able to count the number of 
    characters to determine which line(s) the motif is located on
    """
	try:
		input_handle = open(input_file)
	except:
		return 'Please enter a valid fasta file'
	
	output_file = open('file_for_counting.txt', 'w')
	with input_handle, output_file:
		for line in input_handle:
			if not line.startswith('>'):
				line = line.strip()
				output_file.write(line)
				
	return 'file_for_counting.txt made in working directory'

def motif_finder(input_file, motif):
    """
    This function is designed to assess the input_file for the desired motif and 
    return how many times the motif is found within the input_file, even if 
    there is overlap. 
    """
    try:
    	input_handle = open(input_file)
    except:
        return 'Please enter a valid fasta file'
    
    with input_handle:
    	seqs = ''
    	count = 0
    	pos1 = 0
    	lnct = 0
    	d1 = {}
    	u_motif = motif.strip().upper()
    	n_motif = u_motif.replace('N', '.')
    	for line in input_handle:
    		if not line.startswith('>'):
    			line = line.strip()
    			seqs += line
    			lnct += 1
    			lol = len(line)
    	while True:
    		search_motif = re.search(n_motif, seqs[pos1:])
    		if search_motif is None:
    			break
    		else:
    			pos1 += search_motif.start() + 1
    			count += 1
    			
    	return f'{count} instance(s) of motif identified'
                
def line_finder(input_file, motif):
	"""
    This function is designed to assess the input_file for the desired motif and 
    return the location(s) of the motif within the input_file, even if 
    there is overlap.
    """
    
	try:
		input_handle = open(input_file)
	except:
		return 'Please enter a valid fasta file'
		
	input_file = open('file_for_counting.txt')	
	with input_handle, input_file:
		infile = input_file.read()
		nofc = len(infile)
		seqs = ''
		pos1 = 0
		lines = 0
		u_motif = motif.strip().upper()
		n_motif = u_motif.replace('N', '.')
		for line in input_handle:
			if not line.startswith('>'):
				line = line.strip()
				seqs += line
				lines += 1
		while True:
			search_motif = re.search(n_motif, seqs[pos1:])
			if search_motif is None:
				break
			else:
				pos1 += search_motif.start() + 1
				num_char = nofc/(lines)
				lpos = round((pos1/num_char) + 2)	
				
				
		return f'Your motif is located on line(s): {lpos} in your file'
    	
if __name__ == '__main__':
    main()