// HMMSUM written in CPP
// Author: Haidar Khan
// Date: 3/5/2015
// HMMSTR Work under Prof. Chris Bystroff
// This code was written to integrate HMMSTR and HMMSUM into CLUSTAL
// The goal is to use HMMSUM pairwise alignments in CLUSTAL-W

// HEADER files
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//# of HMMSTR states including "naught" state
#define NUMNODES 282
//# of HMMSTR states
#define NUMSTATES 281
//psuedocount to prevent division by zero
#define PSUEDO 0.000000001

//I heard this was bad but it saves me alot of typing
using namespace std;

//external fortran glue code to get HMMSTR to generate GCA files from SEQ files
extern "C" {
  void seq2gca_(char *seqfile, int *nres, long slen);
}
//CPP wrapper for above fortran subroutine
int seq2gca(string);

//returns the alignment matrix between seq1 and seq2 given the seqnames, obsMat and expMat
vector<vector<double> > getAlignMat(string seqName1, string seqName2,
                                    vector<vector<vector<double > > > obsMat,
                                    vector<vector<double> > expMat);
// returns the alignment matrix between seq1 and seq2 given the two sequences, 
// gamma values for both sequences, obsMat, and expMat
vector<vector<double> > makeMij(string seq1, string seq2,
                                vector<vector<double> > gamma1,
                                vector<vector<double> > gamma2,
                                vector<vector<vector<double> > > obsMat,
                                vector<vector<double> > expMat);
// Reads a GCA file and returns a gamma matrix
vector<vector<double> > readGamma(string, int);
// Reads Expected frequencies files and 
// returns Expected frequencies (281 x 20) matrix
vector<vector<double> > readExpected(string);
// Reads Observed frequencies file and 
// returns Observed frequencies (281 x 20 x 20) matrix
vector<vector<vector<double> > > readObserved(string);
// Reads a FASTA file and splits it into SEQ files (1 sequence per SEQ file)
// Returns names of all the sequences in the FASTA file in a vector
vector<string> readFasta(string);
//reads a SEQ file and returns the sequence in a string
string readSeq(string);
//cleans up all the generated .gca and .seq files
void cleanupFiles(vector<string>);
