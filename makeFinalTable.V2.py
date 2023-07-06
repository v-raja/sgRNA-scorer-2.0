from __future__ import division

import sys
import argparse
import os
from collections import defaultdict
import subprocess
import numpy as np
from smtplib import SMTP
from email.mime.text import MIMEText
from Bio import SeqIO
from Bio.Seq import Seq

# major function call that creates output file
def makeOutputFile(gRNAFile,pamOrientation,svmOutput,outputFile):
	# now get the scores of the sequences
	scoreArray = defaultdict(float)
	for line in svmOutput:
		line = line.rstrip('\r\n')
		parts = line.split('\t')
		scoreArray[parts[0]] = parts[1] 
	svmOutput.close()

	# write the output file
	outputFile.write('SeqID\tSequence\tScore\n')

	for record in SeqIO.parse(gRNAFile,'fasta'):
		if pamOrientation=='5':
			stw = str(record.seq.reverse_complement())
		else:
			stw = str(record.seq)
		outputFile.write(str(record.id) + '\t' + stw + '\t' + str(scoreArray[str(record.id)]) + '\n')

	# close file handles
	outputFile.close()
	gRNAFile.close()
	return

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g','--gRNAFile', type=argparse.FileType('r'),required=True)
	parser.add_argument('-s','--svmOutput', type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--outputFile',type=argparse.FileType('w'),required=True)
	parser.add_argument('-p','--pamOrientation',required=True)
	opts = parser.parse_args(argv)
	makeOutputFile(opts.gRNAFile,opts.pamOrientation,opts.svmOutput,opts.outputFile)

if __name__ == '__main__':
    main(sys.argv[1:])
