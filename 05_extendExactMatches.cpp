#include <iostream>
#include <fstream>
#include <algorithm>

#include <unistd.h>
#include <stdlib.h>
#include <ctime>

#include <string>
#include <sstream>
#include <vector>

using namespace std;

// if we encounter more than this many non-ACGT characters in a row, stop extending
const unsigned int INVALID_SPAN = 5;
// for sanity purposes, lets not consider any edit distance threshold above this
const unsigned int MAX_EDT = 50000;

void strComp(string &x)
{
	for (unsigned int i=0; i<x.size(); i++)
	{
		if (x[i] == 'A')
			x[i] = 'T';
		else if (x[i] == 'C')
			x[i] = 'G';
		else if (x[i] == 'G')
			x[i] = 'C';
		else if (x[i] == 'T')
			x[i] = 'A';
	}
}

// string contains only non-ACGT characters
bool isInvalidString(const string& s) {
	unsigned int totalCount = 0;
	for (unsigned int i=0; i<s.size(); i++)
		if (s[i] != 'A' && s[i] != 'C' && s[i] != 'G' && s[i] != 'T')
			totalCount++;
	return totalCount == s.size();
}

// some code from stack overflow for splitting c-strings
template<typename Out>
void split(const string &s, char delim, Out result)
{
	stringstream ss;
	ss.str(s);
	string item;
	while (getline(ss, item, delim))
	{
		*(result++) = item;
	}
}

vector<string> split(const string &s, char delim)
{
	vector<string> elems;
	split(s, delim, std::back_inserter(elems));
	return elems;
}

unsigned int min_3(unsigned int x, unsigned int y, unsigned int z)
{
	return x < y ? min(x,z) : min(y,z);
}

//
// edit distance, modified from:
//
// http://codereview.stackexchange.com/questions/10130/edit-distance-between-two-strings
//
// -- this version returns the final row/column of the DP matrix. assumes strings are same size
//
vector<unsigned int> edit_distance(const string& A, const string& B)
{
	unsigned int NA = A.size();
	unsigned int NB = B.size();
	vector< vector<unsigned int> > M(NA + 1, vector<unsigned int>(NB + 1));
	for (unsigned int a = 0; a <= NA; ++a)
		M[a][0] = a;
	for (unsigned int b = 0; b <= NB; ++b)
		M[0][b] = b;
	for (unsigned int a = 1; a <= NA; ++a)
		for (unsigned int b = 1; b <= NB; ++b)
		{
			unsigned int x = M[a-1][b] + 1;
			unsigned int y = M[a][b-1] + 1;
			unsigned int z = M[a-1][b-1] + (A[a-1] == B[b-1] ? 0 : 1);
			M[a][b] = min_3(x,y,z);
		}
	//return M[A.size()][B.size()];
	vector<unsigned int> M_out(A.size()+B.size()+1);
	for (unsigned int i=0; i<A.size()+1; i++)
		M_out[i] = M[A.size()][i];
	for (unsigned int i=0; i<A.size(); i++)
	{
		//cout << A.size()+1+i << " --> [" << A.size()-1-i << "][" << A.size() << "]" << endl;
		M_out[A.size()+1+i] = M[A.size()-1-i][A.size()];
	}
	return M_out;
}

//
//	MAIN()
//
int main(int argc, char *argv[])
{
	//
	//  PARSE INPUT ARGUMENTS
	//
	int c;
	char *refFile = NULL;
	char *refName = NULL;
	char *inSeeds = NULL;
	char *outDir  = NULL;
	unsigned int kmerLen  = 0;
	unsigned int paddLen  = 0;
	unsigned int editDist = 0;
	unsigned int jobID    = 0;
	unsigned int jobTot   = 0;
	unsigned int fixedErr = 0;
	unsigned int minEMLen = 0;
	while ((c = getopt(argc, argv, "r:R:i:o:k:p:e:j:J:F:M:")) != -1)
	{
		switch (c)
		{
			case 'r':
				refFile = optarg;
				break;
			case 'R':
				refName = optarg;
				break;
			case 'i':
				inSeeds = optarg;
				break;
			case 'o':
				outDir  = optarg;
				break;
			case 'k':
				kmerLen = atoi(optarg);
				break;
			case 'p':
				paddLen = atoi(optarg);
				break;
			case 'e':
				editDist = atoi(optarg);
				break;
			case 'j':
				jobID = atoi(optarg);
				break;
			case 'J':
				jobTot = atoi(optarg);
				break;
			case 'F':
				fixedErr = atoi(optarg);
				break;
			case 'M':
				minEMLen = atoi(optarg);
				break;
			case '?':
			default:
				fprintf(stderr, "Usage: %s [-r ref.fa -R refName -i input_exactMatches.txt -o outputDir/ -k kmerLen -p paddingLen -e editDist -F fixed_error_percent -M minSeedSize -j jobID -J jobTotal]\n",argv[0]);
				exit(1);
		}
	}

	//
	//  DETERMINE INPUT MODE
	//
	bool fixed_error_mode;
	if (fixedErr > 0)
		fixed_error_mode = true;
	else
		fixed_error_mode = false;

	//
	//  READ INPUT REFERENCE FASTA
	//
	ifstream inStream(refFile);
	string inLine, ref, ref_RC, padding;
	for (unsigned int i=0; i<paddLen; i++)
		padding += "N";
	int nSeq = 0;
	unsigned long long int refLen = 0, paddedRefLen;
	vector<string> refNames;
	vector<unsigned long long int> refIndices;
	while (getline(inStream, inLine))
	{
		if (inLine[0] != '>')
		{
			transform(inLine.begin(), inLine.end(), inLine.begin(), ::toupper);
			string inLineComp(inLine);
			strComp(inLineComp);			
			ref_RC += inLineComp;
			// check for exceeding ref.maxsize() here?
			ref    += inLine;
			refLen += inLine.length();
		}
		else
		{
			ref    += padding;
			ref_RC += padding;
			nSeq++;
			refLen += paddLen;
			refNames.push_back(inLine.substr(1,inLine.length()-1));
			refIndices.push_back(refLen);
		}
	}
	inStream.close();
	reverse(ref_RC.begin(), ref_RC.end());
	paddedRefLen = ref.length();
	cout << refFile << "\n" << paddedRefLen - nSeq*padding.length() << " bp\n";


	//
	//  READ INPUT SEEDS
	//
	ifstream inStream2(inSeeds);
	string soi, *r1, *r2, rs1, rs2, s_temp1, s_temp2;
	vector<string> sop;
	unsigned long long int s1, s1_out;
	unsigned long long int e1, e1_out;
	unsigned long long int s2, s2_out;
	unsigned long long int e2, e2_out;
	unsigned long long int s1_temp, s2_temp, e1_temp, e2_temp;
	bool rc1, rc2;
	
	vector<unsigned int> dp_vals(2*editDist+1), dp_vals_next(2*editDist+1), match_mask(2*editDist+1);
	unsigned int currentDist, cEDT, cEDT_current, tvm;	// cEDT = current edit distance threshold
	unsigned int myR1_invalid, myR2_invalid;

	unsigned int slidingWin = min(editDist,INVALID_SPAN);
	long int endVal = -1;

	cout << "Writing outfile..." << endl;
	string fName(outDir);
	ostringstream oss;
	if (fName[fName.size()-1] != '/')
		fName += '/';
	if (fixed_error_mode == true)
		oss << fName << refName << "_extendedSeeds_F" << fixedErr << "_job_" << jobID << "_of_" << jobTot << ".txt";
	else
		oss << fName << refName << "_extendedSeeds_e" << editDist << "_job_" << jobID << "_of_" << jobTot << ".txt";
	fName = oss.str();
	ofstream outF;
	outF.open(fName.c_str());

	unsigned long long int linesSoFar = 0, linesProcessed = 0, printEvery = 10000;

	while (getline(inStream2, soi))
	{
		// read input coordinates from text-file
		sop = split(soi,'\t');
		linesSoFar++;
		//cout << soi << flush;
		if (sop.size() == 4 && linesSoFar%jobTot == jobID-1)
		{
			linesProcessed++;
			if (linesProcessed%printEvery == 0)
				cout << "extended seeds: " << linesProcessed << endl;
			s1 = strtoull(sop[0].c_str(),0,10);
			e1 = strtoull(sop[1].c_str(),0,10);
			s2 = strtoull(sop[2].c_str(),0,10);
			e2 = strtoull(sop[3].c_str(),0,10);
		}
		else
			continue;
		if (s1 >= paddedRefLen)
		{
			rc1 = true;
			r1  = &ref_RC;
			s1_temp = s1 - paddedRefLen;
			e1_temp = e1 - paddedRefLen;
			s1 = paddedRefLen - e1_temp;
			e1 = paddedRefLen - s1_temp;
		}
		else
		{
			rc1 = false;
			r1 = &ref;
		}
		if (s2 >= paddedRefLen)
		{
			rc2 = true;
			r2  = &ref_RC;
			s2_temp = s2 - paddedRefLen;
			e2_temp = e2 - paddedRefLen;
			s2 = paddedRefLen - e2_temp;
			e2 = paddedRefLen - s2_temp;
		}
		else
		{
			rc2 = false;
			r2 = &ref;
		}

		// skip if below minimum extended_match length
		if (minEMLen > 0 && e1-s1 < minEMLen)
			continue;

		// determine initial edit distance if in fixed error mode
		if (fixed_error_mode == true)
		{
			cEDT = ((e1-s1)*fixedErr)/100; // (s1,e1) & (s2,e2) are same size, so we can use either one
			slidingWin = min(cEDT,INVALID_SPAN);
			dp_vals.resize(2*cEDT+1);
			dp_vals_next.resize(2*cEDT+1);
			match_mask.resize(2*cEDT+1);
		}
		else
		{
			cEDT = editDist;
		}
		cEDT = min(cEDT,MAX_EDT);
		//cout << "=== " << e1-s1 << " " << cEDT << endl;
		//cout << s1 << " " << e1 << " " << s2 << " " << e2 << endl;

		// bounds check to ensure we are at least cEDT away from ref boundaries
		cEDT = min_3(cEDT,s1,s2);
		cEDT = min_3(cEDT,paddedRefLen-e1,paddedRefLen-e2);

		// for easy cases lets spit out the answer with no computation necessary
		if (cEDT == 0)
		{
			s1_out = s1;
			s2_out = s2;
			e1_out = e1;
			e2_out = e2;
		}
		else
		{

			/////////////////////////////////////////////////////////
			//
			//                  EXTEND TO THE RIGHT
			//
			/////////////////////////////////////////////////////////

			// check for invalid patch of characters
			endVal = -1;
			for (unsigned int i=0; i<=cEDT-slidingWin; i++)
			{
				if (e1+i+slidingWin >= paddedRefLen || e2+i+slidingWin >= paddedRefLen)
				{
					endVal = paddedRefLen - max(e1,e2); // check for going past end of sequence
					break;
				}
				if (isInvalidString(r1->substr(e1+i,slidingWin)) || isInvalidString(r2->substr(e2+i,slidingWin)))
				{
					endVal = i;
					break;
				}
			}
			if (endVal >= 0)
			{
				e1_out = e1 + endVal;
				e2_out = e2 + endVal;
			}
			else
			{
				// get dynamic programming matrix values to initialize seed extension algorithm
				dp_vals = edit_distance(r1->substr(e1,cEDT),r2->substr(e2,cEDT));
				e1_out  = e1+cEDT;
				e2_out  = e2+cEDT;
				myR1_invalid = 0;
				myR2_invalid = 0;
				// normal edit distance
				//currentDist = dp_vals[cEDT];

				// lowest score is the minimum of last row/col of DP matrix because
				// we are fine with ignoring trailing deletions in the alignment
				currentDist = *min_element(dp_vals.begin(),dp_vals.end());

				cEDT_current = cEDT;

				//
				// repeatedly extend the optimal alignment until we reach the specified edit distance
				//
				while (true)
				{
					// bounds check
					if (e1_out >= paddedRefLen-1 || e2_out >= paddedRefLen-1)
						break;
					// check again for string of invalid values...
					if ((*r1)[e1_out] != 'A' && (*r1)[e1_out] != 'C' && (*r1)[e1_out] != 'G' && (*r1)[e1_out] != 'T')
						myR1_invalid++;
					else
						myR1_invalid = 0;
					if ((*r2)[e2_out] != 'A' && (*r2)[e2_out] != 'C' && (*r2)[e2_out] != 'G' && (*r2)[e2_out] != 'T')
						myR2_invalid++;
					else
						myR2_invalid = 0;
					if (myR1_invalid >= INVALID_SPAN || myR2_invalid >= INVALID_SPAN)
					{
						e1_out -= INVALID_SPAN;
						e2_out -= INVALID_SPAN;
						break;
					}

					// construct vector of relevant bp match diagonals	
					for (unsigned int i=0; i<cEDT; i++)
					{
						match_mask[i] = ((*r1)[e1_out] == (*r2)[e2_out-cEDT+i] ? 0 : 1);
					}
					for (int i=cEDT; i>=0; i--)
					{
						match_mask[cEDT+i] = ((*r1)[e1_out-i] == (*r2)[e2_out] ? 0 : 1);
					}

					// construct the new row/col of DP matrix
					dp_vals_next[0] = min(dp_vals[0]+match_mask[0], dp_vals[1]+1);
					dp_vals_next[2*cEDT] = min(dp_vals[2*cEDT]+match_mask[2*cEDT], dp_vals[2*cEDT-1]+1);
					for (unsigned int i=1; i<cEDT; i++)
					{
						dp_vals_next[i] = min_3(dp_vals_next[i-1]+1, dp_vals[i]+match_mask[i], dp_vals[i+1]+1);
					}
					for (unsigned int i=2*cEDT-1; i>=cEDT+1; i--)
					{
						dp_vals_next[i] = min_3(dp_vals[i-1]+1, dp_vals[i]+match_mask[i], dp_vals_next[i+1]+1);
					}
					dp_vals_next[cEDT] = min_3(dp_vals[cEDT]+match_mask[cEDT], dp_vals_next[cEDT-1]+1, dp_vals_next[cEDT+1]+1);

					dp_vals = dp_vals_next;
					//currentDist = dp_vals[cEDT];
					currentDist = *min_element(dp_vals.begin(),dp_vals.end());

					// if in fixed error mode, do we need to update the threshold?
					if (fixed_error_mode == true && currentDist >= cEDT_current)
					{
						tvm = max(e1_out-s1,e2_out-s2);
						cEDT_current = (tvm*fixedErr)/100;
						cEDT_current = min(cEDT_current,MAX_EDT);
					}

					if (currentDist < cEDT_current)
					{
						e1_out++;
						e2_out++;
					}
					else
					{
						e1_out -= max(myR1_invalid,myR2_invalid);
						e2_out -= max(myR1_invalid,myR2_invalid);
						break;
					}
				}
			}

			/////////////////////////////////////////////////////////
			//
			//                  EXTEND TO THE LEFT
			//
			// severe amounts of code duplication, for performance.
			//
			/////////////////////////////////////////////////////////

			// check for invalid patch of characters
			endVal = -1;
			for (unsigned int i=0; i<=cEDT-slidingWin; i++)
			{
				if (s1-i-1 <= 0 || s2-i-1 <= 0)
				{
					endVal = -2;
					break;
				}
				if (isInvalidString(r1->substr(s1-i-1,slidingWin)) || isInvalidString(r2->substr(s2-i-1,slidingWin)))
				{
					endVal = i;
					break;
				}
			}
			if (endVal >= 0)
			{
				s1_out = s1 - endVal - 1 + slidingWin;
				s2_out = s2 - endVal - 1 + slidingWin;
			}
			else if (endVal == -2)	// bounds check
			{
				s1_out = s1 - min(s1,s2) + 1;
				s2_out = s2 - min(s1,s2) + 1;
			}
			else
			{
				// get dynamic programming matrix values to initialize seed extension algorithm
				s_temp1 = r1->substr(s1-cEDT,cEDT);
				s_temp2 = r2->substr(s2-cEDT,cEDT);
				reverse(s_temp1.begin(), s_temp1.end());
				reverse(s_temp2.begin(), s_temp2.end());
				dp_vals = edit_distance(s_temp1,s_temp2);
				s1_out  = s1-cEDT;
				s2_out  = s2-cEDT;
				myR1_invalid = 0;
				myR2_invalid = 0;
				currentDist = *min_element(dp_vals.begin(),dp_vals.end());

				cEDT_current = cEDT;

				//
				// repeatedly extend the optimal alignment until we reach the specified edit distance
				//
				while (true)
				{
					// bounds check
					if (s1_out <= 0 || s2_out <= 0)
						break;
					// check again for string of invalid values...
					if ((*r1)[s1_out-1] != 'A' && (*r1)[s1_out-1] != 'C' && (*r1)[s1_out-1] != 'G' && (*r1)[s1_out-1] != 'T')
						myR1_invalid++;
					else
						myR1_invalid = 0;
					if ((*r2)[s2_out-1] != 'A' && (*r2)[s2_out-1] != 'C' && (*r2)[s2_out-1] != 'G' && (*r2)[s2_out-1] != 'T')
						myR2_invalid++;
					else
						myR2_invalid = 0;
					if (myR1_invalid >= INVALID_SPAN || myR2_invalid >= INVALID_SPAN)
					{
						s1_out += INVALID_SPAN-1;
						s2_out += INVALID_SPAN-1;
						break;
					}

					// construct vector of relevant bp match diagonals	
					for (unsigned int i=0; i<cEDT; i++)
					{
						match_mask[i] = ((*r1)[s1_out-1] == (*r2)[s2_out+cEDT-1-i] ? 0 : 1);
					}
					for (int i=cEDT; i>=0; i--)
					{
						match_mask[cEDT+i] = ((*r1)[s1_out-1+i] == (*r2)[s2_out-1] ? 0 : 1);
					}

					// construct the new row/col of DP matrix
					dp_vals_next[0] = min(dp_vals[0]+match_mask[0], dp_vals[1]+1);
					dp_vals_next[2*cEDT] = min(dp_vals[2*cEDT]+match_mask[2*cEDT], dp_vals[2*cEDT-1]+1);
					for (unsigned int i=1; i<cEDT; i++)
					{
						dp_vals_next[i] = min_3(dp_vals_next[i-1]+1, dp_vals[i]+match_mask[i], dp_vals[i+1]+1);
					}
					for (unsigned int i=2*cEDT-1; i>=cEDT+1; i--)
					{
						dp_vals_next[i] = min_3(dp_vals[i-1]+1, dp_vals[i]+match_mask[i], dp_vals_next[i+1]+1);
					}
					dp_vals_next[cEDT] = min_3(dp_vals[cEDT]+match_mask[cEDT], dp_vals_next[cEDT-1]+1, dp_vals_next[cEDT+1]+1);

					dp_vals = dp_vals_next;
					currentDist = *min_element(dp_vals.begin(),dp_vals.end());

					// if in fixed error mode, do we need to update the threshold?
					if (fixed_error_mode == true && currentDist >= cEDT_current)
					{
						tvm = max(e1-s1_out,e2-s2_out);
						cEDT_current = (tvm*fixedErr)/100;
						cEDT_current = min(cEDT_current,MAX_EDT);
					}

					if (currentDist < cEDT_current)
					{
						s1_out--;
						s2_out--;
					}
					else
					{
						s1_out += max(myR1_invalid,myR2_invalid);
						s2_out += max(myR1_invalid,myR2_invalid);
						break;
					}
				}
			}
		}
		//cout << max(e1_out-s1_out,e2_out-s2_out) << " " << cEDT << " " << (max(e1_out-s1_out,e2_out-s2_out)*fixedErr)/100 << endl;

		// convert coordinates back to RC, if necessary
		if (rc1)
		{
			e1_temp = paddedRefLen - s1_out;
			s1_temp = paddedRefLen - e1_out;
			e1_out = e1_temp + paddedRefLen;
			s1_out = s1_temp + paddedRefLen;
		}
		if (rc2)
		{
			e2_temp = paddedRefLen - s2_out;
			s2_temp = paddedRefLen - e2_out;
			e2_out = e2_temp + paddedRefLen;
			s2_out = s2_temp + paddedRefLen;
		}

		// write to outfile
		outF << s1_out << "\t" << e1_out << "\t" << s2_out << "\t" << e2_out << "\n";

		//break;
		
	}

	outF.close();
	inStream2.close();

	return 0;
}

