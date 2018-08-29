import os
import sys
import re
import time
import argparse
import numpy as np


"""*****************************************
********            MAIN            ********
*****************************************"""


def main():

	#	Parse input arguments
	#
	parser = argparse.ArgumentParser(description='track.dat --> track.bed')
	parser.add_argument('-i', type=str, required=True, metavar='<str>', help="input_track.dat")
	parser.add_argument('-r', type=str, required=True, metavar='<str>', help="ref.fa")
	parser.add_argument('-b', type=str, required=True, metavar='<str>', help="comma-separated list of bins for bed output values", default=None)
	parser.add_argument('-o', type=str, required=True, metavar='<str>', help="output.bed")
	args = parser.parse_args()

	OUT_FILE   = args.o
	BED_BINS   = sorted([int(n) for n in args.b.split(',')])

	REF_SEQ    = args.r
	splt = args.i.split('/')
	if len(splt) == 1:
		IND = './'
		INF = args.i
	else:
		IND = '/'.join(splt[:-1])+'/'
		INF = splt[-1]
	PAD_LEN   = 100

	inputTrack = np.fromfile(args.i,'<u4')


	#	Read in ref and notate position of each sequence
	#
	sys.stdout.write('reading in ref...\n')
	sys.stdout.flush()
	f = open(REF_SEQ,'r')
	ref = ''
	refNames = []
	PADDING  = ''.join(['N']*PAD_LEN)
	nSeq = 0
	for line in f:
		if line[0] == '>':
			refNames.append((len(ref),'padding'))
			ref += PADDING
			refNames.append((len(ref),line[1:-1]))
			nSeq += 1
		else:
			ref += line[:-1].upper();
	f.close()
	refLens = {}
	for i in xrange(len(refNames)-1):
		if refNames[i][1] != 'padding':
			refLens[refNames[i][1]] = refNames[i+1][0]-refNames[i][0]
	refLens[refNames[-1][1]] = len(ref)-refNames[-1][0]
	REF_LEN = len(ref)
	ref = ''


	print 'generating bed output...'
	tt = time.time()
	f = open(OUT_FILE,'w')
	chrPos  = [n[0] for n in refNames if n[1] != 'padding']
	chrName = [n[1] for n in refNames if n[1] != 'padding']
	currentChr = -1
	nextChrAt  = chrPos[0]
	regStart   = 0
	prevVal    = 0
	PRINT_EVERY = 1000000
	for i in xrange(len(inputTrack)):
		if (i+1)%PRINT_EVERY == 0:
			print i+1,'/',len(inputTrack),'[{0:.2%}]'.format(i/float(len(inputTrack)))
			#break
		if i == nextChrAt:
			currentChr += 1
			if currentChr < len(chrName)-1:
				nextChrAt = chrPos[currentChr+1]
			regStart = 0
			prevVal  = 0
			#print i, chrName[currentChr]
		if inputTrack[i] >= BED_BINS[0]:
			myBedVal = None
			for bb in BED_BINS:
				if inputTrack[i] >= bb:
					myBedVal = bb
				else:
					break
			if myBedVal != prevVal:
				if prevVal > 0:
					#bedRegion = (chrName[currentChr],regStart,i,prevVal)
					#print bedRegion
					f.write(chrName[currentChr] + '\t' + str(regStart-chrPos[currentChr]) + '\t' + str(i-chrPos[currentChr]) + '\t' + str(prevVal) + '\n')
				regStart = i
				prevVal = myBedVal
		else:
			if prevVal > 0:
				#bedRegion = (chrName[currentChr],regStart,i,prevVal)
				#print bedRegion
				f.write(chrName[currentChr] + '\t' + str(regStart-chrPos[currentChr]) + '\t' + str(i-chrPos[currentChr]) + '\t' + str(prevVal) + '\n')
				prevVal = 0
	f.close()
	print OUT_FILE

	print time.time()-tt,'(sec)'


if __name__ == '__main__':
	main()

