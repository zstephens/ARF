import os
import argparse
from collections import defaultdict

# filter all junk characters in reference to Ns
REF_FILTER = defaultdict(lambda:'N')
REF_FILTER['A'] = 'A'
REF_FILTER['C'] = 'C'
REF_FILTER['G'] = 'G'
REF_FILTER['T'] = 'T'

# RC dict
# add an option to make this RNA friendly?
REV_COMP = defaultdict(lambda:'N')
REV_COMP['A'] = 'T'
REV_COMP['C'] = 'G'
REV_COMP['G'] = 'C'
REV_COMP['T'] = 'A'

def makedir(d):
	if not os.path.isdir(d):
		print 'creating:',d
		os.system('mkdir '+d)

def main():

	parser = argparse.ArgumentParser(description='NEAT-genReads V2.0')
	parser.add_argument('-r',    type=str, required=True,  metavar='<str>',  help="* ref.fa")
	parser.add_argument('-R',    type=str, required=True,  metavar='<str>',  help="* reference name (e.g. hg19)")
	parser.add_argument('-i',    type=str, required=True,  metavar='<str>',  help="* input seeds_pruned.txt")
	parser.add_argument('-o',    type=str, required=True,  metavar='<str>',  help="* outputDir/")
	parser.add_argument('-k',    type=int, required=True,  metavar='<int>',  help="* seed kmer length")
	parser.add_argument('-p',    type=int, required=True,  metavar='<int>',  help="* reference padding length")
	parser.add_argument('--job', type=int, required=False, metavar=('<int>','<int>'), default=(1,1), help='ids for parallel jobs', nargs=2)
	args = parser.parse_args()

	(REF_FILE, REF_NAME, IN_FILE, OUT_DIR, SEED_KMER, PAD_LEN) = (args.r, args.R, args.i, args.o, args.k, args.p)
	(JOB_ID, JOB_TOT) = args.job

	if OUT_DIR[-1] != '/':
		OUT_DIR += '/'
	makedir(OUT_DIR)
	infile_name = IN_FILE.split('/')[-1]
	outfile_name = '_'.join(infile_name.split('_')[:-1])+'_exactMatches_job_'+str(JOB_ID)+'_of_'+str(JOB_TOT)+'.txt'
	# make sure we can write to a valid outfile before doing all the expensive work
	ouf = open(OUT_DIR+outfile_name,'w')
	ouf.close()

	#
	#	READ IN REFERENCE SEQ AND CREATE REVERSE-COMPLEMENT
	#
	print 'reading ref...'
	ref = ''
	PADDING = ''.join(['N']*PAD_LEN)
	f = open(REF_FILE,'r')
	for line in f:
		if line[0] == '>':
			ref += PADDING
		else:
			soi = line[:-1].upper()
			for i in xrange(len(soi)):
				ref += REF_FILTER[soi[i]]
	f.close()
	ref_RC = ''
	for i in xrange(len(ref)-1,-1,-1):
		ref_RC += REV_COMP[ref[i]]
	print 'done.'

	RLEN = len(ref)
	print '\n'+REF_FILE
	print len(ref),'bp\n'

	#
	#	READ IN SEEDS
	#
	all_seeds = []
	f = open(IN_FILE,'r')
	myLine = 0
	totalDict = {}
	totalTot  = 0
	for line in f:
		myLine += 1
		if myLine%JOB_TOT != JOB_ID - 1:
			continue
		print '[',myLine,']'
		v = [int(n) for n in line.strip().split('\t')]
		subDict = {}
		subTot  = 0
		for i in xrange(len(v)-1):
			for j in xrange(i+1,len(v)):
				isRC = [v[i]>=len(ref),v[j]>=len(ref)]
				if all(isRC):	# skip if both in RC (because it's equivalent to both not in RC)
					continue
				if isRC[0]:
					myRef1 = ref_RC
					s1     = 2*len(ref)-v[i]-SEED_KMER
					e1     = 2*len(ref)-v[i]
				else:
					myRef1 = ref
					s1     = v[i]
					e1     = v[i] + SEED_KMER
				if isRC[1]:
					myRef2 = ref_RC
					s2     = 2*len(ref)-v[j]-SEED_KMER
					e2     = 2*len(ref)-v[j]
				else:
					myRef2 = ref
					s2     = v[j]
					e2     = v[j] + SEED_KMER

				# extend to the left
				while True:
					if s1 == 1 or s2 == 1:
						break
					if myRef1[s1-1] != myRef2[s2-1] or myRef1[s1-1] == 'N' or myRef2[s2-1] == 'N':
						break
					s1 -= 1
					s2 -= 1

				# extend to the right
				while True:
					if e1 >= len(ref)-1 or e2 >= len(ref)-1:
						break
					if myRef1[e1] != myRef2[e2] or myRef1[e1] == 'N' or myRef2[e2] == 'N':
						break
					e1 += 1
					e2 += 1

				# if either index was in RC, lets correct it back to len(ref)+ind
				if isRC[0]:
					e1_temp = len(ref)-s1
					s1_temp = len(ref)-e1
					e1 = e1_temp+len(ref)
					s1 = s1_temp+len(ref)
				if isRC[1]:
					e2_temp = len(ref)-s2
					s2_temp = len(ref)-e2
					e2 = e2_temp+len(ref)
					s2 = s2_temp+len(ref)

				# add to dictionary
				if s1 < s2:
					subDict[(s1,e1,s2,e2)] = True
				else:
					subDict[(s2,e2,s1,e1)] = True
				subTot += 1

		#print subTot,'-->',len(subDict)
		# 808396 --> 517538
		totalTot  += subTot
		for k in subDict.keys():
			totalDict[k] = True 

	f.close()

	print totalTot,'-->',len(totalDict)
	sorted_keys = sorted([k for k in totalDict.keys()])
	ouf = open(OUT_DIR+outfile_name,'a')
	for k in sorted_keys:
		ouf.write('\t'.join([str(n) for n in k])+'\n')
	ouf.close()


if __name__ == '__main__':
	main()

