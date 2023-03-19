#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Feb 11 09:25:02 2023

This is a module designed to take a pdb file, restrict the residues supplied 
by the user,run the pose through mutational analysis, and determine the 
delta delta G of each mutation. 

@author: breewenck
"""

#Python
from pyrosetta import *
import pyrosetta; pyrosetta.init()
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import cleanATOM
import sys
import argparse
from rosetta import *


#Core Includes
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.core.select.movemap import *

#Protocols
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn
from pyrosetta.rosetta.protocols import minimization_packing as pack_min


from pyrosetta.rosetta.core import * 
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.core.scoring import * 
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.protocols.minimization_packing import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.conformation import *
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.toolbox.mutants import mutate_residue


def main():
	
	pdb_file, restrict, minimize, non_protein = arg_parse()
	if restrict == 0:
		print(clean_up(pdb_file, non_protein))	
	else:
		print(clean_up(pdb_file, non_protein))
		rpose = pose_from_pdb(pdb_file + '.relax.pdb')
		print(ptp_file(rpose))
		print(restrict_pose_residues(rpose, restrict, pdb_file, minimize))
		
				
def arg_parse():
    """
    This function is designed to add pdb_file and restrict arguments and parse 
    to the function from the command line
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('pdb_file', 
    					help="PDB file name without .png, must be in working directory (ex. 1hpv and NOT 1hpv.png)") 
    parser.add_argument('-r', '--restrict', nargs='+', 
    					help="Enter 0 for no restrictions or enter the residues to restrict from mutation (ex. 32A, 24B - an integer and letter)")
    parser.add_argument('-m', '--minimize',
    					help='perform minimization after each mutation? yes or no')
    parser.add_argument('-rnp', '--non_protein', help='remove non-protein? yes or no')
    
    args = parser.parse_args()
    
    return args.pdb_file, args.restrict, args.minimize, args.non_protein 

def clean_up(pdb_file, non_protein):
	"""
    This function is designed to clean up the pdb file, relax it, and get it ready
    for restrictions and/or mutations.
    """
	cleanATOM(pdb_file + '.pdb')
	pose = pose_from_pdb(pdb_file + '.clean.pdb')
	if non_protein == 'yes':
		pyrosetta.rosetta.core.pose.remove_nonprotein_residues(pose)
	sfxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.constrain_relax_to_start_coords(True)
	relax.set_scorefxn(sfxn)
	relax.apply(pose)
	pose.dump_pdb(pdb_file + '.relax.pdb')


def cat(in1, in2, in3 = '', in4 = ''):
	"""
    This function is designed to concatenate a minimum of two up to four values.
    """
	return f'{in1}{in2}{in3}{in4}'
	
def pdb_to_pose_dict(k, val): 
	"""
    This function is designed to build a dictionary for PDB:pose reference
    """
	d = {k: val}
	return d

def ptp_file(rpose):
	"""
	This function is designed to create a text file to reference PDB to pose values; 
	determine which pose residues (integers) correspond to the PDB residues in the file 
	to restrict from mutations.
	"""
	count = 0
	output_file = open('pdb_to_pose_values.txt', 'w')	
	with output_file:
		for i in range (1, rpose.total_residue() + 1):	
			res = rpose.pdb_info().number(i)
			c = rpose.pdb_info().chain(i) 
			count += 1
			k = cat(res,c)
			d = pdb_to_pose_dict(k, count)
			output_file.write(f'{d}\n')
			
if __name__ == '__main__':
	main()