import sys
import re
import bisect
import copy

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

f = open(in_file,'r')
allReps = []
for line in f:
	splt = sorted([int(n) for n in line.split('\t')])
	allReps.append(splt)
f.close()
print 'done.'

allReps = sorted(allReps)

print 'removing redundant seeds that will extend to same repeat...'
print len(allReps),'-->',

prevGroup = []
prunedReps = []
for i in xrange(len(allReps)):
	n = copy.deepcopy(allReps[i])
	if i == 0:
		prunedReps.append(n)
		prevGroup = n
		continue
	if allReps[i] == allReps[i-1]:
		continue
	gDelta = n[0]-prevGroup[0]
	if gDelta < SEED_KMER:
		temp = [n[0]]
		for j in xrange(1,len(n)):
			myMinusOneIsInLastGroup = False
			for k in xrange(1,len(prevGroup)):
				if n[j] == prevGroup[k]+gDelta:
					myMinusOneIsInLastGroup = True
			myPlusOneIsInLastGroupAndImRC = False
			if n[j] >= REF_LEN:
				for k in xrange(1,len(prevGroup)):
					if n[j] == prevGroup[k]-gDelta:
						myPlusOneIsInLastGroupAndImRC = True
			if (not myMinusOneIsInLastGroup) and (not myPlusOneIsInLastGroupAndImRC):
				temp.append(n[j])
		if len(temp) > 1:
			if any([(m<REF_LEN) for m in temp]):
				prunedReps.append(temp)
	else:
		if any([(m<REF_LEN) for m in n]):
			prunedReps.append(n)
	prevGroup = n
print len(prunedReps)

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

if in_file[-4:] == '.txt':
	out_file1 = in_file[:-4]+'_pruned.txt'
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
f = open(out_file2,'w')
for i in xrange(len(lowCSeeds)):
	for j in xrange(len(lowCSeeds[i])-1):
		f.write(str(lowCSeeds[i][j])+'\t')
	f.write(str(lowCSeeds[i][-1])+'\n')
f.close()

