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
	args = parser.parse_args()


	MULTIPLICITY = False
	MAPPABILITY  = False
	NON_UNIQUE   = False
	if '_multiplicityTrack.dat' in args.i:
		MULTIPLICITY = True
	elif '_mappabilityTrack.dat' in args.i or '_mapabilityTrack.dat' in args.i:	# early versions of neat had this typo. ugh...
		MAPPABILITY  = True
	elif '_nonUniqueTrack.dat' in args.i:
		NON_UNIQUE   = True
	else:
		print '\nError: Invalid track type specified.\n'
		exit(1)

	BED_BINS   = sorted([int(n) for n in args.b.split(',')])

	REF_SEQ    = args.r
	splt = args.i.split('/')
	if len(splt) == 1:
		IND = './'
		INF = args.i
	else:
		IND = '/'.join(splt[:-1])+'/'
		INF = splt[-1]
	SEED_KMER = int(re.findall(r"_k\d+",INF)[0][2:])
	PAD_LEN   = int(re.findall(r"_p\d+",INF)[0][2:])
	EDIT_DIST = int(re.findall(r"_e\d+",INF)[0][2:])

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
	fn = args.i[:-4] + '.bed'
	f = open(fn,'w')
	chrPos  = [n[0] for n in refNames if n[1] != 'padding']
	chrName = [n[1] for n in refNames if n[1] != 'padding']
	currentChr = -1
	nextChrAt  = chrPos[0]
	regStart   = 0
	prevVal    = 0
	PRINT_EVERY = 1000000
	for i in xrange(len(inputTrack)):
		if i%PRINT_EVERY == 0:
			print i,'/',len(inputTrack),'[{0:.2%}]'.format(i/float(len(inputTrack)))
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
	print fn

	print time.time()-tt,'(sec)'


if __name__ == '__main__':
	main()

