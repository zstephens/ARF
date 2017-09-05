import sys
import re
import bisect
import copy

SEPARATE_LOWCOMPLEXITY = False
REMOVE_PERIODIC        = True

from collections import defaultdict
RC_DICT = defaultdict(lambda:'N')
RC_DICT['A'] = 'T'
RC_DICT['C'] = 'G'
RC_DICT['G'] = 'C'
RC_DICT['T'] = 'A'
def RC(s):
	return ''.join([RC_DICT[n] for n in s[::-1]])

# detect periodicty (<=6) of a kmer
def getPeriodicity(s):
	for period in range(1,7):
		isPeriodic = True
		pattern = s[:period]
		for i in xrange(period,len(s),period):
			pat2 = s[i:i+period]
			if pat2 != pattern[:len(pat2)]:
				isPeriodic = False
				break
		if isPeriodic:
			return period
	return None

if len(sys.argv) != 3:
	print '\nusage:\n\npython pruneSeeds.py ref.fa input.txt\n'
	exit(1)
else:
	ref_file  = sys.argv[1]
	in_file   = sys.argv[2]
	SEED_KMER = int(re.findall(r"_k\d+",in_file)[0][2:])
	PAD_LEN   = int(re.findall(r"_p\d+",in_file)[0][2:])

print 'reading ref...'
f = open(ref_file,'r')
ref = ''
PADDING  = ''.join(['N']*PAD_LEN)
nSeq = 0
for line in f:
	if line[0] == '>':
		ref += PADDING
		nSeq += 1
	else:
		ref += line[:-1].upper();
f.close()
REF_LEN = len(ref)
#print 'REF_LEN:',REF_LEN
#REF_LEN = 119668450

print 'reading input seed regions...'
f = open(in_file,'r')
allReps = []
for line in f:
	splt = sorted([int(n) for n in line.split('\t')])
	allReps.append(splt)
f.close()
allReps = sorted(allReps)

print 'removing redundant seeds that will extend to same repeat...'
print len(allReps),'-->',

#prevGroup = []
prunedReps = []
for i in xrange(len(allReps)):
	if i == 0 or i == len(allReps)-1:
		prunedReps.append(allReps[i])
		#prevGroup = copy.deepcopy(allReps[i])
		continue
	if allReps[i] == allReps[i-1]:
		continue

	# check for periodicity
	if REMOVE_PERIODIC:
		if len(allReps[i]) >= 5:
			stretches = [[],[],[]]	# 1, 2, 3
			for gapLen in xrange(1,len(stretches)+1):
				currentStart = 0
				currentlyIn  = (allReps[i][1] - allReps[i][0] == gapLen)
				for j in xrange(1,len(allReps[i])):
					myGap = allReps[i][j] - allReps[i][j-1]
					if myGap == gapLen:
						if currentlyIn == False:
							currentStart = j
						currentlyIn = True
					else:
						if currentlyIn == True:
							stretches[gapLen-1].append([currentStart,j])
						currentlyIn = False
			indsWeDontNeed = []
			indsWeNeed     = []
			for gapInd in xrange(len(stretches)):
				for stretchInd in xrange(len(stretches[gapInd])):
					if stretches[gapInd][stretchInd][1] - stretches[gapInd][stretchInd][0] + 1 >= 5:
						indsWeNeed.extend([stretches[gapInd][stretchInd][0],stretches[gapInd][stretchInd][0]+1,stretches[gapInd][stretchInd][1],stretches[gapInd][stretchInd][1]-1])
						indsWeDontNeed.extend(range(stretches[gapInd][stretchInd][0]+2,stretches[gapInd][stretchInd][1]-1))
			indsWeDontNeed = sorted(list(set(indsWeDontNeed)))
			indsWeNeed     = sorted(list(set(indsWeNeed)))

			delList = []
			for j in xrange(len(allReps[i])):
				if j in indsWeDontNeed and j not in indsWeNeed:
					delList.append(j)
			delList = sorted(delList,reverse=True)
			for j in delList:
				del allReps[i][j]

	hasPrevMate = [False for n in allReps[i]]
	hasNextMate = [False for n in allReps[i]]
	gDelta1 = allReps[i][0] - allReps[i-1][0]
	gDelta2 = allReps[i+1][0] - allReps[i][0]

	for j in xrange(len(allReps[i])):
		for k in xrange(len(allReps[i-1])):
			myDelta = -1
			if allReps[i][j] < REF_LEN and allReps[i-1][k] < REF_LEN:
				myDelta = allReps[i][j] - allReps[i-1][k]
			elif allReps[i][j] >= REF_LEN and allReps[i-1][k] >= REF_LEN:
				myDelta = allReps[i-1][k] - allReps[i][j]
			if myDelta == gDelta1 and myDelta > 0 and myDelta < SEED_KMER:
				hasPrevMate[j] = True

		for k in xrange(len(allReps[i+1])):
			myDelta = -1
			if allReps[i+1][k] < REF_LEN and allReps[i][j] < REF_LEN:
				myDelta = allReps[i+1][k] - allReps[i][j]
			elif allReps[i+1][k] >= REF_LEN and allReps[i][j] >= REF_LEN:
				myDelta = allReps[i][j] - allReps[i+1][k]
			if myDelta == gDelta2 and myDelta > 0 and myDelta < SEED_KMER:
				hasNextMate[j] = True

	if all(hasPrevMate) and all(hasNextMate):
		pass
	else:
		prunedReps.append(allReps[i])

prunedReps = sorted(prunedReps)

print len(prunedReps)

if SEPARATE_LOWCOMPLEXITY:
	print 'separating low complexity seeds...'
	easySeeds = []
	lowCSeeds = []
	for n in prunedReps:
		if n[0] >= REF_LEN:
			ind = n[0] - REF_LEN
			s = RC(ref[ind:ind+SEED_KMER])
		else:
			ind = n[0]
			s = ref[ind:ind+SEED_KMER]

		myPer = getPeriodicity(s)
		if myPer != None:
			lowCSeeds.append([m for m in n])
		else:
			easySeeds.append([m for m in n])
	print len(prunedReps),'-->',len(easySeeds),'+',len(lowCSeeds),'(low)'
else:
	easySeeds = prunedReps
	lowCSeeds = []

if in_file[-4:] == '.txt':
	out_file1 = in_file[:-4]+'_pruned.txt'
	if SEPARATE_LOWCOMPLEXITY:
		out_file2 = in_file[:-4]+'_lowComplexity.txt'
else:
	print 'input file should be .txt'
	exit(1)
f = open(out_file1,'w')
for i in xrange(len(easySeeds)):
	for j in xrange(len(easySeeds[i])-1):
		f.write(str(easySeeds[i][j])+'\t')
	f.write(str(easySeeds[i][-1])+'\n')
f.close()
if SEPARATE_LOWCOMPLEXITY:
	f = open(out_file2,'w')
	for i in xrange(len(lowCSeeds)):
		for j in xrange(len(lowCSeeds[i])-1):
			f.write(str(lowCSeeds[i][j])+'\t')
		f.write(str(lowCSeeds[i][-1])+'\n')
	f.close()

