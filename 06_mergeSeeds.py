import os
import sys
import re
import time
import argparse
import numpy as np

from scipy.sparse import coo_matrix
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx

UINT4_MAX = 2**32 - 1
REV_COMP = {'A':'T','C':'G','G':'C','T':'A','N':'N'}

def matToVec(i,j,N):
	if (i <= j):
		return i * N - (i - 1) * i / 2 + j - i
	else:
		return j * N - (j - 1) * j / 2 + i - j

def vecSize(N):
	return N*(N+1)/2

def vecToMat2(x,N):
	vs = N*(N+1)/2
	triRoot  = int((sqrt(8*(vs-1-x)+1)-1)/2)
	i = N - triRoot - 1
	j = N - vs + x + (triRoot*(triRoot+1))/2
	return (i,j)

COLOR_BINS = [-1,75,100,500,1000,5000,10000]

def getColor(i,N,colormap='jet'):
	cm = mpl.get_cmap(colormap) 
	cNorm  = colors.Normalize(vmin=0, vmax=N+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	colorVal = scalarMap.to_rgba(i)
	return colorVal

def grouper(iterable, diff):
	prev = None
	group = []
	for item in iterable:
		if not prev or item - prev <= diff:
			group.append(item)
		else:
			yield group
			group = [item]
		prev = item
	if group:
		yield group

def plot_coo_matrix(m):
	if not isinstance(m, coo_matrix):
		m = coo_matrix(m)
	fig = mpl.figure()
	ax = fig.add_subplot(111, axisbg='white')
	ax.plot(m.col, m.row, '.', color='black', ms=0.01)
	ax.set_xlim(0, m.shape[1])
	ax.set_ylim(0, m.shape[0])
	ax.set_aspect('equal')
	for spine in ax.spines.values():
		spine.set_visible(False)
	ax.invert_yaxis()
	ax.set_aspect('equal')
	#ax.set_xticks([])
	#ax.set_yticks([])
	return ax

def plot_PosDict(d,nBin,PNG_MAP):
	# Using contourf to provide my colorbar info, then clearing the figure
	Z = [[0,0],[0,0]]
	levels = range(len(COLOR_BINS)+1)
	CS3 = mpl.contourf(Z, levels, cmap=mpl.get_cmap('jet'))
	mpl.clf()

	if PNG_MAP:
		fig = mpl.figure(figsize=(100,80))
		mpl.rcParams.update({'font.size': 100, 'font.weight':'bold', 'lines.linewidth': 16})
	else:
		fig = mpl.figure(figsize=(20,16))
	ax  = fig.add_subplot(111, axisbg='white')

	for cbin in xrange(len(COLOR_BINS)-1):
		rowDat   = []
		colDat   = []
		cchr     = chr(cbin)
		colorVal = getColor(cbin,len(COLOR_BINS))
		for k in d.keys():
			if d[k] == cchr:
				(myRow,myCol) = vecToMat2(k,nBin)
				rowDat.append(myRow)
				colDat.append(myCol)
				del d[k]
		#print cbin, len(rowDat)
		if len(rowDat):
			if PNG_MAP:
				ax.plot(colDat, rowDat, '.', ms=1.00 + 1.00*cbin, markerfacecolor=colorVal, markeredgecolor=colorVal)
			else:
				ax.plot(colDat, rowDat, '.', ms=0.10 + 0.10*cbin, markerfacecolor=colorVal, markeredgecolor=colorVal)

	rowDat   = []
	colDat   = []
	colorVal = getColor(len(COLOR_BINS)-1,len(COLOR_BINS))
	for k in d.keys():
		(myRow,myCol) = vecToMat2(k,nBin)
		rowDat.append(myRow)
		colDat.append(myCol)
	if len(rowDat):
		if PNG_MAP:
			ax.plot(colDat, rowDat, '.', ms=1.00 + 1.00*(len(COLOR_BINS)-1), markerfacecolor=colorVal, markeredgecolor=colorVal)
		else:
			ax.plot(colDat, rowDat, '.', ms=0.10 + 0.10*(len(COLOR_BINS)-1), markerfacecolor=colorVal, markeredgecolor=colorVal)
	
	ax.set_xlim(0, nBin)
	ax.set_ylim(0, nBin)
	ax.set_aspect('equal')
	for spine in ax.spines.values():
		spine.set_visible(False)
	ax.invert_yaxis()
	ax.set_aspect('equal')

	cb = mpl.colorbar(CS3,fraction=0.046, pad=0.04)
	cb.set_ticks(np.arange(0,len(COLOR_BINS),1)+0.5)
	ctick = ['< '+str(COLOR_BINS[1]),] + ['> '+str(COLOR_BINS[i]) for i in xrange(1,len(COLOR_BINS))]
	cb.set_ticklabels(ctick)
	cb.set_label('repeat size (nt)')

	return (fig,ax)


SORT_AFTER_EVERY = 10000000
BP_BIN           = 50000

def condenseListOfRegions(l):

	delList = [False]
	prevR   = l[0]
	for i in xrange(1,len(l)):
		if l[i][1] <= prevR[1]:
			delList.append(True)
		else:
			delList.append(False)
			prevR = l[i]
	condensedInput = [l[i] for i in xrange(len(l)) if delList[i] == False]
	prevR = condensedInput[-1]
	delList = [False]
	for i in xrange(len(condensedInput)-2,-1,-1):
		if condensedInput[i][0] == prevR[0] and condensedInput[i][1] <= prevR[1]:
			delList.append(True)
		else:
			delList.append(False)
			prevR = condensedInput[i]
	delList.reverse()

	return [condensedInput[i] for i in xrange(len(condensedInput)) if delList[i] == False]

def r_intersect(a,b):
	return (a[0] <= b[1] and b[0] <= a[1])

"""*****************************************
********            MAIN            ********
*****************************************"""


def main():

	#	Parse input arguments
	#
	parser = argparse.ArgumentParser(description='Extended Seeds -> Plots + Tracks + Etc')
	parser.add_argument('-i', type=str, required=True,  metavar='<str>', help="input_directory/")
	parser.add_argument('-r', type=str, required=True,  metavar='<str>', help="ref.fa")
	parser.add_argument('-R', type=str, required=True,  metavar='<str>', help="ref name")
	parser.add_argument('-o', type=str, required=True,  metavar='<str>', help="output_prefix")
	parser.add_argument('-k', type=int, required=True,  metavar='<int>', help="K: seed kmer length")
	parser.add_argument('-p', type=int, required=True,  metavar='<int>', help="P: ref padding length")
	#parser.add_argument('-e', type=int, required=True,  metavar='<int>', help="E: edit distance used")
	parser.add_argument('-b', type=str, required=False, metavar='<str>', help="comma-separated list of bins for bed output values", default=None)
	#parser.add_argument('-l', type=str, required=False, metavar='<str>', help="lowComplexitySeeds.txt", default=None)
	parser.add_argument('--map',        required=False, action='store_true', help='produce repeat map', default=False)
	parser.add_argument('--merge',      required=False, action='store_true', help='produce merged repeat track', default=False)
	parser.add_argument('--mult', type=str,      required=False, metavar='<str>', help='produce multiplicity track with this repeat track', default=None)
	parser.add_argument('--mappability',required=False, action='store_true', help='compute mappability instead of non-uniqueness', default=False)
	parser.add_argument('--png',        required=False, action='store_true', help='output png map instead of pdf', default=False)
	parser.add_argument('--no-par', type=str,      required=False, metavar='<str>', help='exclude X/Y PAR regions ("hg19" or "hg38")', default='')
	args = parser.parse_args()

	MAKE_MAP     = args.map
	PNG_MAP      = args.png
	OUTPUT_MERGE = args.merge
	MAPPABILITY  = args.mappability
	if args.mult == None:
		MULTIPLICITY = False
	else:
		MULTIPLICITY = True
		inputRepeatTrack = np.fromfile(args.mult,'<u4')
		indexPairs = {}
	if not(MAKE_MAP) and not(OUTPUT_MERGE) and not(MULTIPLICITY):
		print "\nYou didn't give me anything to do. Choose --map, --merge, or --mult.\n"
		exit(1)
	if (MAKE_MAP and OUTPUT_MERGE) or (OUTPUT_MERGE and MULTIPLICITY) or (MULTIPLICITY and MAKE_MAP):
		print "\nCan only do one of the following at a time: --map, --merge, or --mult.\n"
		exit(1)

	MAKE_BED = False
	if args.b != None:
		MAKE_BED = True
		BED_BINS = sorted([int(n) for n in args.b.split(',')])

	no_par = args.no_par
	if no_par not in ['','hg19','hg38']:
		print '\nError: --no-par accepts only hg19 or hg38.\n'
		exit(1)

	REF_SEQ    = args.r
	REF_NAME   = args.R
	INF        = args.i
	if INF[-1] != '/':
		INF += '/'
	SEED_KMER  = args.k
	PAD_LEN    = args.p
	#EDIT_DIST  = args.e
	OUT_PREFIX = args.o
	LOW_COMPLEXITY_FILE = None

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

	chrX_start = [n[0] for n in refNames if n[1] == 'chrX'][0]
	chrY_start = [n[0] for n in refNames if n[1] == 'chrY'][0]
	if no_par == 'hg38':
		PAR_XY = [(chrX_start+10000,chrX_start+2781480), (chrY_start+10000,chrY_start+2781480)]
		PAR_XY+= [(chrX_start+155701382,chrX_start+156030896), (chrY_start+56887902,chrY_start+57217416)]
	elif no_par == 'hg19':
		PAR_XY = [(chrX_start+60000,chrX_start+2699521), (chrY_start+10000,chrY_start+2649521)]
		PAR_XY+= [(chrX_start+154931043,chrX_start+155260561), (chrY_start+59034049,chrY_start+59363567)]

	#	Make sure all jobs are present
	#
	listing = [n for n in os.listdir(INF) if ('_job_' in n and 'extendedSeeds' in n)]
	jobList = [int(n.split('_')[-3]) for n in listing]
	jobTot  = int(listing[0].split('_')[-1][:-4])
	anyMissing = False
	for i in range(1,jobTot+1):
		if i not in jobList:
			print '\nError: Missing data for job #'+str(i)+'\n'
			anyMissing = True
	if anyMissing:
		exit(1)
	#listing = [INF]

	if MAKE_MAP:
		N        = REF_LEN/BP_BIN + 1
		repSizes = {}

	if OUTPUT_MERGE:
		bigRegionList    = []
		newItemsToInsert = {}

	#	Read in extended seeds
	tt = time.time()
	firstTimeThrough = True
	for li in xrange(len(listing)):

		fn = listing[li]
		sys.stdout.write('reading file: '+fn+'\n')
		sys.stdout.flush()

		f = open(INF+fn,'r')
		everyOther = False
		prevRegion = None
		nPairs     = 0
		for line in f:
			####splt   = line.strip().split(' ')
			####if everyOther:
			####	region = (int(splt[0]),int(splt[1]))
			####	isRC = [(region[0]>=REF_LEN),(prevRegion[0]>=REF_LEN)]
			####	if isRC[0]:
			####		region = (region[0]-REF_LEN, region[1]-REF_LEN)
			####	if isRC[1]:
			####		prevRegion = (prevRegion[0]-REF_LEN, prevRegion[1]-REF_LEN)
			####else:
			####	prevRegion = (int(splt[0]),int(splt[1]))

			splt = line.strip().split('\t')
			prevRegion = (int(splt[0]),int(splt[1]))
			region = (int(splt[2]),int(splt[3]))
			isRC = [(region[0]>=REF_LEN),(prevRegion[0]>=REF_LEN)]
			if isRC[0]:
				region = (region[0]-REF_LEN, region[1]-REF_LEN)
			if isRC[1]:
				prevRegion = (prevRegion[0]-REF_LEN, prevRegion[1]-REF_LEN)

			# process region pairs
			if True:

				regionPairs = [(region,prevRegion)]
				####if firstTimeThrough:
				####	firstTimeThrough = False
				####	if LOW_COMPLEXITY_FILE != None:
				####		lc_regions = []
				####		f2 = open(LOW_COMPLEXITY_FILE,'r')
				####		for line2 in f2:
				####			splt2  = sorted([int(n) for n in line2.strip().split('\t')])
				####			groups = dict(enumerate(grouper(splt2,SEED_KMER), 1))
				####			for v in groups.values():
				####				if v[0] >= REF_LEN:
				####					lc_regions.append((v[0]-REF_LEN,v[-1]+SEED_KMER-REF_LEN))
				####				else:
				####					lc_regions.append((v[0],v[-1]+SEED_KMER))
				####		f2.close()
				####		for i in xrange(len(lc_regions)-1):
				####			for j in xrange(i+1,len(lc_regions)):
				####				regionPairs.append((lc_regions[i],lc_regions[j]))

				for regionPair in regionPairs:
					if regionPair[0][0] < regionPair[1][0]:
						(region,prevRegion) = regionPair
					else:
						(prevRegion,region) = regionPair

					#print '- r, pr:',region, prevRegion

					# region validation
					if region[0] < 0 or prevRegion[0] < 0 or region[1] >= REF_LEN or prevRegion[1] >= REF_LEN:
						print 'Warning: Skipping invalid regions:', prevRegion, region
						continue
					hasN = False
					#if 'N' in ref[region[0]:region[1]] or 'N' in ref[prevRegion[0]:prevRegion[1]]:
					#	hasN = True
					hasOL = False
					if prevRegion[0] < region[1]:
						hasOL = True
					# ignore regions that overlap
					if hasOL:
						continue

					#if (hasN or hasOL) and region[1]-region[0] > 10000:
					#	print (hasN, hasOL), region, prevRegion, region[1]-region[0],
					#	print ref[region[0]:region[1]].count('N'),
					#	print ref[prevRegion[0]:prevRegion[1]].count('N')
					#	print ''
					#	print ref[region[0]:region[0]+50]
					#	print ref[region[1]-50:region[1]]
					#	print ''
					#	print ref[prevRegion[0]:prevRegion[0]+50]
					#	print ref[prevRegion[1]-50:prevRegion[1]]
					#	print ''

					rLen = region[1] - region[0]
					nPairs += 1

					# exclude repeats corresponding to PAR regions in chrX/Y, if desired...
					# shortcut: assume X comes before Y
					if no_par == 'hg19' or no_par == 'hg38':
						if r_intersect(PAR_XY[0],region) and r_intersect(PAR_XY[1],prevRegion):
							print 'skipping PAR1:', prevRegion, region
							continue
						if r_intersect(PAR_XY[2],region) and r_intersect(PAR_XY[3],prevRegion):
							print 'skipping PAR2:', prevRegion, region
							continue

					if MAKE_MAP:
						myColor = -1
						for ci in COLOR_BINS:
							if rLen >= ci:
								myColor += 1
						bitRange1 = xrange(prevRegion[0]/BP_BIN, prevRegion[1]/BP_BIN + 1)
						bitRange2 = xrange(region[0]/BP_BIN, region[1]/BP_BIN + 1)
						for b1 in bitRange1:
							for b2 in bitRange2:
								if b1 > b2:
									(b1,b2) = (b2,b1)
								myInd = matToVec(b1,b2,N)
								if myInd in repSizes:
									if myColor > repSizes[myInd]:
										repSizes[myInd] = chr(myColor)
								else:
									repSizes[myInd] = chr(myColor)

					if OUTPUT_MERGE:
						newItemsToInsert[prevRegion] = True
						newItemsToInsert[region] = True
						if nPairs%SORT_AFTER_EVERY == 0:
							print nPairs,'region pairs read. list size:',len(bigRegionList),'-->',
							bigRegionList.extend(condenseListOfRegions(sorted(newItemsToInsert.keys())))
							bigRegionList = condenseListOfRegions(sorted(bigRegionList))
							print len(bigRegionList)
							newItemsToInsert = {}

					if MULTIPLICITY:
						for ii in xrange(region[1] - region[0]):
							ri = region[0] + ii
							oi = prevRegion[0] + ii
							if rLen == inputRepeatTrack[ri]:
								if ri not in indexPairs:
									indexPairs[ri] = set()
								indexPairs[ri].add(oi)

							if rLen == inputRepeatTrack[oi]:
								if oi not in indexPairs:
									indexPairs[oi] = set()
								indexPairs[oi].add(ri)

				everyOther = False
			else:
				everyOther = True
		#break

	""" ***********************************************************
				PROCESS ALL REGIONS TO PRODUCE TRACKS
	*********************************************************** """

	if MAKE_MAP:
		print 'plotting self-similarity map...'
		xt = [n[0]/BP_BIN for n in refNames if n[1] != 'padding']
		xl = [n[1] for n in refNames if n[1] != 'padding']
		for i in xrange(len(xl)):
			if not len(re.findall(r"chr\d+",xl[i])+re.findall(r"Chr\d+",xl[i])):
				xl[i] = ''
		prevChr = xt[-1]
		for i in xrange(len(xt)-2,-1,-1):
			if xt[i] == prevChr:
				xl[i] += ' '+xl[i+1]
				del xt[i+1]
				del xl[i+1]
			else:
				prevChr = xt[i]
		(fig,ax) = plot_PosDict(repSizes,N,PNG_MAP)
		mpl.xticks(xt,xl,rotation=70)
		mpl.yticks(xt,xl)
		mpl.grid(color='black',linewidth=1,linestyle='-')

		if PNG_MAP:
			fn = OUT_PREFIX + '_selfSimilarityMap.png'
		else:
			fn = OUT_PREFIX + '_selfSimilarityMap.pdf'
		print fn
		mpl.savefig(fn)

	if OUTPUT_MERGE or MULTIPLICITY:

		if MULTIPLICITY:
			print 'creating multiplicity track...'
			repCov = np.zeros(REF_LEN,dtype='<u4')
			multList_by_repSize = {}
			for ii in sorted(indexPairs.keys()):
				repCov[ii] = len(indexPairs[ii]) + 1
				myR = inputRepeatTrack[ii]
				if myR not in multList_by_repSize:
					multList_by_repSize[myR] = []
				multList_by_repSize[myR].append(len(indexPairs[ii]) + 1)

			fn = OUT_PREFIX + '_multiplicityTrack.dat'

		elif OUTPUT_MERGE:
			# One final merge
			print nPairs,'region pairs read. list size:',len(bigRegionList),'-->',
			bigRegionList.extend(condenseListOfRegions(sorted(newItemsToInsert.keys())))
			bigRegionList = condenseListOfRegions(sorted(bigRegionList))
			print len(bigRegionList)

			print 'sorting repeats by length...'
			bigRegionList = sorted([(n[1]-n[0],n) for n in bigRegionList])

			repCov = np.zeros(REF_LEN,dtype='<u4')
			if MAPPABILITY:
				print 'creating mappability track...'
				# this will numerically overflow if repeat size is > 2^32. But that will never realistically happen..
				for n in bigRegionList:
					if n[0] > 2 * SEED_KMER:
						indsToDo = n[0] - 2 * SEED_KMER
						for i in xrange(SEED_KMER,SEED_KMER+(indsToDo+1)/2):
							if repCov[n[1][0]+i] < i:
								repCov[n[1][0]+i] = i
						for i in xrange(SEED_KMER,SEED_KMER+indsToDo/2):
							if repCov[n[1][1]-1-i] < i:
								repCov[n[1][1]-1-i] = i
				fn = OUT_PREFIX + '_mappabilityTrack.dat'

			else:
				print 'creating repeat track...'
				for n in bigRegionList:
					rVal = min([n[0],UINT4_MAX])
					for i in xrange(n[1][0],n[1][1]):
						repCov[i] = rVal
				fn = OUT_PREFIX + '_nonUniqueTrack.dat'

		print fn
		repCov.tofile(fn)

		if MAKE_BED:
			print 'generating bed output...'
			if MULTIPLICITY:
				fn = OUT_PREFIX + '_multiplicityTrack.bed'
			elif MAPPABILITY:
				fn = OUT_PREFIX + '_mappabilityTrack.bed'
			else:
				fn = OUT_PREFIX + '_nonUniqueTrack.bed'
			f = open(fn,'w')
			chrPos  = [n[0] for n in refNames if n[1] != 'padding']
			chrName = [n[1] for n in refNames if n[1] != 'padding']
			currentChr = -1
			nextChrAt  = chrPos[0]
			regStart   = 0
			prevVal    = 0
			for i in xrange(len(repCov)):
				if i == nextChrAt:
					currentChr += 1
					if currentChr < len(chrName)-1:
						nextChrAt = chrPos[currentChr+1]
					regStart = 0
					prevVal  = 0
					#print i, chrName[currentChr]
				if repCov[i] >= BED_BINS[0]:
					myBedVal = None
					for bb in BED_BINS:
						if repCov[i] >= bb:
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
			print fn

		print 'generating plot data...'

		if MULTIPLICITY:

			xDat = sorted(multList_by_repSize.keys())
			yDat = []
			for k in xDat:
				#multList_by_repSize[k] = np.median(multList_by_repSize[k])
				yDat.append(str(np.median(multList_by_repSize[k]))+' '+str(np.max(multList_by_repSize[k])))

			fn = OUT_PREFIX + '_multiplicity_plotData.txt'

		else:
			countDict = {}
			for i in xrange(len(repCov)):
				if repCov[i] not in countDict:
					countDict[repCov[i]] = 0
				countDict[repCov[i]] += 1

			TOTAL_BP = REF_LEN - nSeq * PAD_LEN
			xDat  = []
			yDat  = []
			prevY = 0.0
			for k in sorted(countDict.keys(),reverse=True):
				if k == 0:
					continue
				xDat.append(k)
				percentage = countDict[k]/float(TOTAL_BP)
				yDat.append(percentage + prevY)
				prevY = yDat[-1]
			xDat = xDat[::-1]
			yDat = yDat[::-1]

			if MAPPABILITY:
				fn = OUT_PREFIX + '_mappability_plotData.txt'
			else:
				fn = OUT_PREFIX + '_nonUnique_plotData.txt'

		print fn
		f = open(fn,'w')
		for i in xrange(len(xDat)):
			f.write(str(xDat[i]) + '\t' + str(yDat[i]) + '\n')
		f.close()

	print time.time()-tt,'(sec)'


if __name__ == '__main__':
	main()

