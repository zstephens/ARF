#include <iostream>
#include <fstream>
#include <algorithm>

#include <unistd.h>
#include <stdlib.h>
#include <ctime>

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

typedef unordered_map< string, vector<unsigned long long int> > kmerInds;

const unsigned int PADDING_SIZE = 100;
const unsigned int BATCH_SIZE   = 10000000;

//
//	USAGE:	./kmerSeeds -r -i myRef.fa -o outFilePrefix -s finalKmerSize
//

string hashKmer(const string kmer)
{
	unsigned int sOut = (kmer.size()/4)+1;
	string kpad = kmer;
	for (unsigned int i=0; i<sOut-kmer.size(); i++)
		kpad += (char)0;
	string outStr;
	for (unsigned int i=0; i<kpad.size(); i += 4)
	{
		outStr += (char)( (kpad[i]<<6) + (kpad[i+1]<<4) + (kpad[i+2]<<2) + (kpad[i+3]) );
	}
	return outStr;
}

string revComp(const string x)
{
	string outStr;
	for (int i=x.size()-1; i>=0; i--)
	{
		if (x[i] == 'A')
			outStr += 'T';
		else if (x[i] == 'C')
			outStr += 'G';
		else if (x[i] == 'G')
			outStr += 'C';
		else if (x[i] == 'T')
			outStr += 'A';
	}
	return outStr;
}

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

int main(int argc, char *argv[])
{

	//
	//  PARSE INPUT ARGUMENTS
	//

	int c;
	char *inFile  = NULL;
	char *outFile = NULL;
	char *jfFile  = NULL;
	while ((c = getopt(argc, argv, "i:j:o:")) != -1)
	{
		switch (c)
		{
			case 'i':
				inFile  = optarg;
				break;
			case 'o':
				outFile = optarg;
				break;
			case 'j':
				jfFile  = optarg;
				break;
			case '?':
			default:
				fprintf(stderr, "Usage: %s [-i ref.fa -j jellyfish_dump.fa -o kmerLoc.txt]\n",argv[0]);
				exit(1);
		}
	}

	//
	//  READ INPUT REFERENCE FASTA
	//

	ifstream inStream(inFile);
	string inLine, ref, ref_RC, padding;
	for (unsigned int i=0; i<PADDING_SIZE; i++)
		padding += "N";
	int nSeq = 0;
	unsigned long long int refLen = 0;
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
			refLen += PADDING_SIZE;
			refNames.push_back(inLine.substr(1,inLine.length()-1));
			refIndices.push_back(refLen);
		}
	}
	inStream.close();
	reverse(ref_RC.begin(), ref_RC.end());
	cout << inFile << "\n" << ref.length() - nSeq*padding.length() << "bp\n";

	//
	//  READ INPUT JELLYFISH FILE / WRITE OUTPUT FILE
	//

	ifstream inStream2(jfFile);
	string soi;
	kmerInds kmerLoc;
	unsigned long long int nKmersDone = 0;
	unsigned int kmerLen = 0;

	cout << "Writing outfile..." << endl;
	ofstream outF;

	bool firstWriteOut  = true;
	unsigned int errors = 0;

	while (getline(inStream2, soi))
	{
		if (soi[0] != '>')
		{
			//string soi_RC = revComp(soi);
			//cout << soi << endl << soi_RC << endl;
			//exit(1);
			vector<unsigned long long int> newVec;
			kmerLoc[soi] = newVec;
			kmerLen = soi.size();

			if (kmerLoc.size() >= BATCH_SIZE)
			{
				if (firstWriteOut)
				{
					firstWriteOut = false;
					string fName(outFile);
					char buffer1[20], buffer2[20];
					sprintf(buffer1, "_k%d",kmerLen);
					sprintf(buffer2, "_p%d.txt", PADDING_SIZE);
					fName += buffer1;
					fName += buffer2;
					outF.open(fName.c_str());
				}
				unsigned int currentPercent = 0;
				unsigned long long int add1PercentAfterEach = ((ref.length()-kmerLen)/100)-1;
				for (unsigned long long int i=0; i<ref.length()-kmerLen; i++)
				{
					if (i%add1PercentAfterEach == 0)
					{
						cout << currentPercent << "% " << flush;
						currentPercent += 1;
					}
					string kmer    = ref.substr(i,kmerLen);
					string kmer_RC = ref_RC.substr(ref.length()-i-kmerLen,kmerLen);

					if (kmer.find_first_of("N")!=std::string::npos)
						continue;

					if (kmerLoc.find(kmer) != kmerLoc.end())
						kmerLoc[kmer].push_back(i);
					if (kmerLoc.find(kmer_RC) != kmerLoc.end())
						kmerLoc[kmer_RC].push_back(ref.length()+i);
				}
				for ( kmerInds::iterator it = kmerLoc.begin(); it != kmerLoc.end(); ++it )
				{
					unsigned int sl = it->second.size();
					if (sl > 1)
					{
						for (unsigned int j=0; j<sl-1; j++)
							outF << it->second[j] << "\t";
						outF << it->second[sl-1] << "\n";
					}
					else
						errors += 1;
				}
				nKmersDone += kmerLoc.size();
				cout << nKmersDone << endl;
				kmerLoc.clear();
			}
		}
	}

	//
	//	blah code duplication, I'm sorry!
	//
	if (firstWriteOut)
	{
		firstWriteOut = false;
		string fName(outFile);
		char buffer1[20], buffer2[20];
		sprintf(buffer1, "_k%d",kmerLen);
		sprintf(buffer2, "_p%d.txt", PADDING_SIZE);
		fName += buffer1;
		fName += buffer2;
		outF.open(fName.c_str());
	}
	unsigned int currentPercent = 0;
	unsigned long long int add1PercentAfterEach = ((ref.length()-kmerLen)/100)-1;
	for (unsigned long long int i=0; i<ref.length()-kmerLen; i++)
	{
		if (i%add1PercentAfterEach == 0)
		{
			cout << currentPercent << "% " << flush;
			currentPercent += 1;
		}

		string kmer    = ref.substr(i,kmerLen);
		string kmer_RC = ref_RC.substr(ref.length()-i-kmerLen,kmerLen);

		if (kmer.find_first_of("N")!=std::string::npos)
			continue;

		if (kmerLoc.find(kmer) != kmerLoc.end())
			kmerLoc[kmer].push_back(i);
		if (kmerLoc.find(kmer_RC) != kmerLoc.end())
			kmerLoc[kmer_RC].push_back(ref.length()+i);
	}
	for ( kmerInds::iterator it = kmerLoc.begin(); it != kmerLoc.end(); ++it )
	{
		unsigned int sl = it->second.size();
		if (sl > 1)
		{
			for (unsigned int j=0; j<sl-1; j++)
				outF << it->second[j] << "\t";
			outF << it->second[sl-1] << "\n";
		}
		else
			errors += 1;
	}
	nKmersDone += kmerLoc.size();
	cout << nKmersDone << endl;
	kmerLoc.clear();

	outF.close();
	inStream2.close();

	cout << "errors: " << errors << endl;

	return 0;
}

