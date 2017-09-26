import sys
import os

PRINT_EVERY = 100000

IN_DIR = sys.argv[1]
if IN_DIR[-1] != '/':
	IN_DIR += '/'

NORMAL_JOBS = True
SPLIT_JOBS  = -1
if len(sys.argv) > 2:
	if sys.argv[2] == 'ignore':
		NORMAL_JOBS = False
	if 'split_' in sys.argv[2]:
		SPLIT_JOBS = int(sys.argv[2].split('_')[1])

if NORMAL_JOBS:
	listing = [n for n in os.listdir(IN_DIR) if ('_job_' in n and 'exactMatches' in n)]
	jobList = [int(n.split('_')[-3]) for n in listing]
	jobTot  = int(listing[0].split('_')[-1][:-4])
	anyMissing = False
	for i in range(1,jobTot+1):
		if i not in jobList:
			print '\nError: Missing data for job #'+str(i)+'\n'
			anyMissing = True
	if anyMissing:
		exit(1)
	if SPLIT_JOBS > 0:
		listing  = [listing[n::SPLIT_JOBS] for n in xrange(SPLIT_JOBS)]
		OUT_FILE = ['_'.join(listing[0][0].split('_')[:-4])+'_split_'+str(n+1)+'.txt' for n in xrange(SPLIT_JOBS)]
	else:
		listing  = [listing]
		OUT_FILE = ['_'.join(listing[0][0].split('_')[:-4])+'.txt']
else:
	listing  = [[n for n in os.listdir(IN_DIR) if n[-4:] == '.txt']]
	OUT_FILE = ['all_raw.txt']

#
#
#
for fff_ind in xrange(len(listing)):
	fList = [open(IN_DIR+n,'r') for n in listing[fff_ind]]
	fInd  = [0 for n in fList]
	fNext = [n.readline() for n in fList]
	myMin = min([int(n.split('\t')[0]) for n in fNext])

	ouf = open(IN_DIR+OUT_FILE[fff_ind],'w')
	while True:
		outDict = {}
		delList = []
		intList = [int(n.split('\t')[0]) for n in fNext]
		myMin   = min(intList)
		for i in xrange(len(fList)):
			if intList[i] == myMin:
				outDict[fNext[i]] = True
				while True:
					fNext[i]    = fList[i].readline()
					if fNext[i] == '':
						delList.append(i)
						break
					intList[i] = int(fNext[i].split('\t')[0])
					if intList[i] == myMin:
						outDict[fNext[i]] = True
					else:
						break
		for n in delList[::-1]:
			del fList[n]
			del fNext[n]
		for k in sorted(outDict.keys()):
			ouf.write(k)
		if len(fList) == 0:
			break
