import sys
import os

IN_DIR = sys.argv[1]
if IN_DIR[-1] != '/':
	IN_DIR += '/'

NORMAL_JOBS = True
if len(sys.argv) > 2:
	if sys.argv[2] == 'ignore':
		NORMAL_JOBS = False

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
	OUT_FILE = '_'.join(listing[0].split('_')[:-4])+'.txt'
else:
	listing  = [n for n in os.listdir(IN_DIR) if n[-4:] == '.txt']
	OUT_FILE = 'all_raw.txt'

fList = [open(IN_DIR+n,'r') for n in listing]
fInd  = [0 for n in fList]

PRINT_EVERY = 100000

fNext = [n.readline() for n in fList]
myMin = min([int(n.split('\t')[0]) for n in fNext])

ouf = open(IN_DIR+OUT_FILE,'w')
while True:
	outDict = {}
	delList = []
	intList = [int(n.split('\t')[0]) for n in fNext]
	myMin   = min(intList)
	for i in xrange(len(fList)):
		if intList[i] == myMin:
			outDict[fNext[i]] = True
			while True:
				fNext[i]   = fList[i].readline()
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
	for k in outDict.keys():
		ouf.write(k)
	if len(fList) == 0:
		break
