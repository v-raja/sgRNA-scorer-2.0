# Date: July 20th, 2016
# Author: Raj Chari
# Script function: generate SVM distances for test sequences based on the model from good and bad
# Input arguments:
# -g (--good) -> sequences of high activity sgRNAs for Sp and St1
# -b (--bad) -> sequences of low activity sgRNAs for Sp and St1
# -i (--input) -> input sequences of the putative sgRNAs identified
# -o (--svmOutput) -> output file with SVM distances for each sequence
# -s (--spacerLength) -> legnth of spacer sequence to work from
# -p (--pamOrientation) -> depending on if the PAM is on the 3' or 5' end will make a difference on how the model is applied

from __future__ import division

import sys
import argparse
import os

from collections import defaultdict
from sklearn.svm import SVC
from Bio import SeqIO

# binary encoding
encoding = defaultdict(str)
encoding['A'] = '0001'
encoding['C'] = '0010'
encoding['T'] = '0100'
encoding['G'] = '1000'

# add encoding for ambiguous bases
encoding['K'] = '1100'
encoding['M'] = '0011'
encoding['R'] = '1001'
encoding['Y'] = '0110'
encoding['S'] = '1010'
encoding['W'] = '0101'
encoding['B'] = '1110'
encoding['V'] = '1011'
encoding['H'] = '0111'
encoding['D'] = '1101'
encoding['N'] = '1111'

def generateSVMOut(goodFile,badFile,inputFile,spacerLength,pamOrientation,pamLength,svmOutputFile):
	# make a giant x list and y list
	xList = []
	yList = []
	offset = 0
	# if the spacer length provided is > 20, only take up to 20 bases
	if int(spacerLength) >= 20:
		spacerLengthInt = 20
		offSetGuide = int(spacerLength) - 20
	else:
		spacerLengthInt = int(spacerLength)
		offSetGuide = 0

	# calculate offSet for SVM
	if int(spacerLength) < 20:
		offSetModel = 20 - spacerLengthInt
	else:
		offSetModel = 0

	# go through each list
	for sequence in goodFile:
		sequence = sequence.rstrip('\r\n')
		entryList = []
		x = offSetModel
		while x < spacerLengthInt + offSetModel:
			y = 0
			while y < 4:
				entryList.append(int(encoding[sequence[x]][y]))
				y += 1
			x += 1
		xList.append(entryList)
		yList.append(1)

	# go through the bad
	for sequence in badFile:
		sequence = sequence.rstrip('\r\n')
		entryList = []
		x = offSetModel
		while x < spacerLengthInt + offSetModel:
			y = 0
			while y < 4:
				entryList.append(int(encoding[sequence[x]][y]))
				y += 1
			x += 1
		xList.append(entryList)
		yList.append(-1)

	# close files
	goodFile.close()
	badFile.close()

	# calculate all SVMs
	clfLinear = SVC(kernel='linear')
	clfLinear.fit(xList,yList)

	## debugging
	maxScore = 0
	minScore = 0

	for record in SeqIO.parse(inputFile,'fasta'):
		sequence = str(record.seq).upper()
		seqRev = sequence[::-1]
		entryList = []
		entryListR = []
		# do this if the PAM is at the 3' end
		x = offSetGuide
		while x < offSetGuide + spacerLengthInt:
			y = 0
			while y < 4:
				entryList.append(int(encoding[sequence[x]][y]))
				entryListR.append(int(encoding[seqRev[x]][y]))
				y += 1
			x += 1
		# predict based on the entry
		myEntry = []
		if pamOrientation=='3':
			myEntry.append(entryList)
		else:
			myEntry.append(entryListR)
		prediction = clfLinear.predict(myEntry)
		score = clfLinear.decision_function(myEntry)

		# print prediction
		svmOutputFile.write(str(record.id) + '\t' + str(score[0]) + '\n')

	# close test file
	inputFile.close()
	svmOutputFile.close()

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g','--good', type=argparse.FileType('r'),required=True)
	parser.add_argument('-b','--bad', type=argparse.FileType('r'),required=True)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True)
	parser.add_argument('-s','--spacerLength',required=True)
	parser.add_argument('-p','--pamOrientation',required=True)
	parser.add_argument('-l','--pamLength',required=True)
	parser.add_argument('-o','--svmOutput',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	generateSVMOut(opts.good,opts.bad,opts.input,opts.spacerLength,opts.pamOrientation,opts.pamLength,opts.svmOutput)

if __name__ == '__main__':
    main(sys.argv[1:])
