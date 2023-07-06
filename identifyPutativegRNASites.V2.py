# script to identify gRNA sites in a list of FASTA files
# Date: October 5th, 2016

from __future__ import division

import sys
import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

# list of valid characters
validCharacters = ['A','C','T','G','K','M','R','Y','S','W','B','V','H','D','N']

# reverse complement table
rcTable = defaultdict(str)
rcTable['A'] = 'T'
rcTable['C'] = 'G'
rcTable['G'] = 'C'
rcTable['T'] = 'A'

def generateRecognitionSites(spacerLength,pamSequence,pamOrientation):
	allSites = []
	sites = []
	sitesRC = []
	# build spacer sequence
	x = 0
	spacer = ''
	while x < int(spacerLength):
		spacer = spacer + '\S'
		x += 1
	# next, make PAM sequences
	for character in pamSequence:
		returnedSites = generateSiteLists(character)
		myList = []
		myListRC = []
		if len(sites)==0:
			sites = returnedSites[0]
			sitesRC = returnedSites[1]
		else:
			for site in sites:
				for site2 in returnedSites[0]:
					newSite = site + site2
					myList.append(newSite)
			for site in sitesRC:
				for site2 in returnedSites[1]:
					newSite = site2 + site
					myListRC.append(newSite)
			sites = myList
			sitesRC = myListRC
	# add to all sites
	targetSites = []
	for site in sites:
		if pamOrientation=='5':
			targetSite = site + spacer
		else:
			targetSite = spacer + site
		targetSites.append(targetSite)
	# do the reverse complement
	targetSitesRC = []
	for site in sitesRC:
		if pamOrientation=='5':
			targetSite = spacer + site
		else:
			targetSite = site + spacer
		targetSitesRC.append(targetSite)	

	# return them
	allSites.append(targetSites)
	allSites.append(targetSitesRC)
	return allSites

# function to generate site lists
def generateSiteLists(character):
	sites = []
	sitesRC = []
	if character=='N':
		site = '\S'
		sites.append(site)
		sitesRC.append(site)
	elif character=='K':
		site1 = 'G'
		site2 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='M':
		site1 = 'A'
		site2 = 'C'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='R':
		site1 = 'A'
		site2 = 'G'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='Y':
		site1 = 'C'
		site2 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='S':
		site1 = 'G'
		site2 = 'C'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='W':
		site1 = 'A'
		site2 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
	elif character=='B':
		site1 = 'C'
		site2 = 'G'
		site3 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
		sites.append(site3)
		sitesRC.append(rcTable[site3])				
	elif character=='V':
		site1 = 'A'
		site2 = 'C'
		site3 = 'G'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
		sites.append(site3)
		sitesRC.append(rcTable[site3])		
	elif character=='H':
		site1 = 'A'
		site2 = 'C'
		site3 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
		sites.append(site3)
		sitesRC.append(rcTable[site3])		
	elif character=='D':
		site1 = 'A'
		site2 = 'G'
		site3 = 'T'
		sites.append(site1)
		sitesRC.append(rcTable[site1])
		sites.append(site2)
		sitesRC.append(rcTable[site2])
		sites.append(site3)
		sitesRC.append(rcTable[site3])		
	else:
		site = character
		siteRC = rcTable[character]
		sites.append(site)
		sitesRC.append(siteRC)
	listReturn = []
	listReturn.append(sites)
	listReturn.append(sitesRC)
	return listReturn

# function to identify sgRNA sites
def identifysgRNASites(inputFASTA,recognitionSites,targetSize,outputFile):
	# initialized variables
	count = 1
	totalgRNAsFound = 0
	gRNASeen = defaultdict(str)
	# target size = spacer length + length of PAM
	sites = recognitionSites[0]
	sitesRC = recognitionSites[1]
	# go through each record and identify target sites
	for record in SeqIO.parse(inputFASTA,'fasta'):
		seqID = record.id
		strSeq = str(record.seq)
		end = len(strSeq) - targetSize
		index = 0
		while index <= end:
			for toFind in sites:
				myIter = re.finditer(toFind,strSeq[index:])
				for gRNA in myIter:
					if 'N' not in str(gRNA.group(0)) and str(gRNA.group(0)) not in gRNASeen:
						outputFile.write('>' + str(seqID) + '_Plus_' + str(gRNA.start() + index) + '\n' + str(gRNA.group(0)) + '\n')
						count = count + 1
						totalgRNAsFound = totalgRNAsFound + 1
						gRNASeen[str(gRNA.group(0))] = 'Y'
			for toFindRC in sitesRC:
				myIterRC = re.finditer(toFindRC,strSeq[index:])
				for gRNA in myIterRC:
					if 'N' not in str(gRNA.group(0)) and str(gRNA.group(0)) not in gRNASeen:
						seq = Seq(gRNA.group(0))
						seqRC = seq.reverse_complement()
						outputFile.write('>' + str(seqID) + '_Minus_' + str(gRNA.start() + index) + '\n' + str(seqRC) + '\n')
						count = count + 1
						totalgRNAsFound = totalgRNAsFound + 1
						gRNASeen[str(gRNA.group(0))] = 'Y'
			index += 1
	outputFile.close()
	return

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i','--inputFASTA',type=argparse.FileType('r'),required=True)
	parser.add_argument('-p','--pamSequence',required=True)
	parser.add_argument('-q','--pamOrientation',required=True)
	parser.add_argument('-s','--spacerLength',required=True)
	parser.add_argument('-o','--outputFile',type=argparse.FileType('w'),required=True)
	opts = parser.parse_args(argv)
	# generate the recognition sites
	print('Generating recognition sites')
	recognitionSites = generateRecognitionSites(opts.spacerLength,opts.pamSequence,opts.pamOrientation)
	targetLength = int(opts.spacerLength) + len(opts.pamSequence)
	print('Using recognition sites to find target sites')
	identifysgRNASites(opts.inputFASTA,recognitionSites,targetLength,opts.outputFile)

if __name__ == '__main__':
    main(sys.argv[1:])