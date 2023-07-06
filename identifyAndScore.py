# Author: Raj Chari
# Date: October 5th, 2016

# Main wrapper script for sgRNA Scorer 2.0 (RUN THIS ONLY, DO NOT RUN THE OTHER THREE SCRIPTS)
from __future__ import division

import sys
import subprocess
import os
import argparse
import platform
from collections import defaultdict
from Bio import SeqIO

def validPAM(pamSequence):
	goodPAM = True
	validCharacters = ['A','C','T','G','K','M','R','Y','S','W','B','V','H','D','N']
	for character in pamSequence:
		if character not in validCharacters:
			goodPAM = False
			break
	return goodPAM


def runPipeline(inputFile,outputFile,spacerLength,pamSequence,pamOrientation):
	# call the identify function
	print('Identifying sgRNAs in input file')
	outputFile1 = inputFile.name.replace('.fasta','.putative.fasta')
	command1 = 'python identifyPutativegRNASites.V2.py -i ' + inputFile.name + ' -p ' + pamSequence + ' -q ' + pamOrientation + ' -s ' + spacerLength + ' -o ' + outputFile1
	p = subprocess.Popen(command1,shell=True)
	p.communicate()
	# next call the SVM function
	print('Classifying identified sgRNA sequences')
	outputFile2 = inputFile.name.replace('.fasta','.SVMOutput.tab')
	command2 = 'python generateSVMFile.V2.py -g Cas9.High.tab -b Cas9.Low.tab -i ' + outputFile1 + ' -s ' + spacerLength + ' -p ' + pamOrientation + ' -l ' + str(len(pamSequence)) + ' -o ' + outputFile2
	p = subprocess.Popen(command2,shell=True)
	p.communicate()
	# finally call the make table function to put it into a table
	print('Making final output file')
	command3 = 'python makeFinalTable.V2.py -g ' + outputFile1 + ' -s ' + outputFile2 + ' -o ' + outputFile + ' -p ' + pamOrientation
	p = subprocess.Popen(command3,shell=True)
	p.communicate()	
	# delete the temporary files
	if platform.system()=='Windows':
		delCommand = 'del '
	else:
		delCommand = 'rm '
	delOutput = delCommand + outputFile1 + ' ' + outputFile2
	p = subprocess.Popen(delOutput,shell=True)
	p.communicate()
	return


def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-s','--spacerLength',required=True)
	parser.add_argument('-p','--pamOrientation',required=True)
	parser.add_argument('-l','--pamSequence',required=True)
	opts = parser.parse_args(argv)
	# initialize variables
	validPAMSequence = ''
	spacerLength = '20'
	pamOrientation = '3'
	fastaDictionary = defaultdict(str)

	# make sure the input is proper, first check PAM
	if validPAM(opts.pamSequence)==False:
		sys.exit('PAM sequence has invalid characters. PAM sequence must only contain A,C,T,G,K,M,R,Y,S,W,B,V,H,D,N')
	else:
		validPAMSequence = opts.pamSequence
	# next check that the spacer is only numbers
	if opts.spacerLength.isdigit()==False:
		sys.exit('Spacer length must be an integer with a minimum size of 14')
	else:
		spacerLength = opts.spacerLength
	# check PAM orientation
	if opts.pamOrientation!='5' and opts.pamOrientation!='3':
		sys.exit('Valid PAM orientations are 5 or 3')
	else:
		pamOrientation = opts.pamOrientation
	# go through the input file
	for record in SeqIO.parse(opts.input,'fasta'):
		fastaDictionary[str(record.id)] = str(record.seq)
	if any(fastaDictionary):
		runPipeline(opts.input,opts.output,spacerLength,validPAMSequence,pamOrientation)
	else:
		sys.exit('Invalid sequence file. Please make sure file is in FASTA format')

if __name__ == '__main__':
    main(sys.argv[1:])