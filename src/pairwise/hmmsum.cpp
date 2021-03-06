// HMMSUM written in CPP
// Author: Haidar Khan
// Date: 3/5/2015
// HMMSTR Work under Prof. Chris Bystroff
// This code was written to integrate HMMSTR and HMMSUM into CLUSTAL
// The goal is to use HMMSUM pairwise alignments in CLUSTAL-W

#include "hmmsum.h"

//testing code for these functions
int fake_main(){
  //load observed and expected matrices;
  string obsfile="observed";
  string expfile="expected";
  vector<vector<double> > expMat = readExpected(expfile);
  vector<vector<vector<double> > > obsMat = readObserved(obsfile);
  //load sequences from FASTA file and split into separate files
  //seqeuence names are used as separate file names (hopefully this is not problematic
  vector<string> SeqsNames = readFasta("small.fas");
  
  //Run sequences through HMMSTR to generate GCA files
  for(int i = 0; i < SeqsNames.size()-1; ++i){
    //read the SEQ file into a string
    string seq1 = readSeq(SeqsNames[i]+".seq");
    //generate GCA file for this sequence
    int nres = seq2gca(SeqsNames[i]+".seq");    
    //read the gamma value for this sequence
    vector<vector<double> > gamma1 = readGamma(SeqsNames[i]+".gca", nres);
    for(int j = i+1; j < SeqsNames.size(); ++j){
      //read the SEQ file into a string
      string seq2 = readSeq(SeqsNames[j]+".seq");
      //generate GCA file for this sequence
      int nres = seq2gca(SeqsNames[j]+".seq");
      //read the gamma value for this sequence
      vector<vector<double> > gamma2 = readGamma(SeqsNames[j]+".gca", nres);
      
      //make the alignment matrix for these two sequences
      vector<vector<double> > alignM = makeMij(seq1, seq2, gamma1, gamma2, obsMat, expMat);
    }
  }

  return 0;
}

vector<vector<double> > getAlignMat(string seqName1, string seqName2,
				    vector<vector<vector<double > > > obsMat,
				    vector<vector<double> > expMat){
  //read the SEQ file into a string
  string seq1 = readSeq(seqName1+".seq");
  //generate GCA file for this sequence
  int nres1 = seq2gca(seqName1+".seq");
  //read the gamma value for this sequence
  vector<vector<double> > gamma1 = readGamma(seqName1+".gca", nres1);

  //read the SEQ file into a string
  string seq2 = readSeq(seqName2+".seq");
  //generate GCA file for this sequence
  int nres2 = seq2gca(seqName2+".seq");
  //read the gamma value for this sequence
  vector<vector<double> > gamma2 = readGamma(seqName2+".gca", nres1);

  //make the alignment matrix for these two sequences
  vector<vector<double> > alignM = makeMij(seq1, seq2, gamma1, gamma2, obsMat, expMat);

  return alignM;

}

vector<vector<double> > makeMij(string seq1, string seq2,
				vector<vector<double> > gamma1, 
				vector<vector<double> > gamma2, 
				vector<vector<vector<double> > > obsMat, 
				vector<vector<double> > expMat){
  //gamma1 - len1 x 282
  //gamma2 - len2 x 282
  //obsMat - 281 x 20 x 20
  //expMat - 281 x 20


  string alphabet = "ACDEFGHIKLMNPQRSTVWYBZX.-~*";
  //initialize the alignment matrix
  vector<vector<double> > alignM;
  int len1 = seq1.length();
  int len2 = seq2.length();
  for(int i = 0; i < NUMSTATES; ++i){  
    for(int j = 1; j < 20; ++j){
      //      for(int k = 0; k < 20; ++k) cout << obsMat[i][j][k] << " ";
    }
    //cout << endl;
  }
  //resize it to (len1 x len2)
  alignM.resize(len1);
  for(int i = 0; i < len1; ++i){
    alignM[i].resize(len2);
  }
  //calculate probability of a residue
  double Pi[20];
  for(int i = 0; i < 20; ++i){
    Pi[i] = 0;
    for(int j = 0; j < NUMSTATES; ++j){
      Pi[i] = Pi[i]+expMat[j][i];
    }
  }
  //make the substitution matrices for these two sequences

  for(int i = 0; i < len1; ++i){
    for(int j = 0; j < len2; ++j){
      alignM[i][j] = 0;      
      size_t p1 = alphabet.find(seq1[i]);
      size_t p2 = alphabet.find(seq2[j]);
      if(p1 >= 20 || p2 >=20 || p1<0 || p2<0){
	/*cerr << "Something went wrong in makeMij...\n";
	cerr << "p1: " << p1 << endl;
	cerr << "p2: " << p2 << endl;
	int x;
	cin >> x;*/
	continue;
      }
      double sumM = 0;
      //gamma's are indexed 1 to 281; hence the (w+1)'s
      for(int w = 0; w < NUMSTATES; ++w){
	if(p1>=p2){
	  alignM[i][j] = alignM[i][j]+gamma1[i][w+1]*obsMat[w][p2][p1]*gamma2[j][w+1];
	}
	else{
	  alignM[i][j] = alignM[i][j]+gamma1[i][w+1]*obsMat[w][p1][p2]*gamma2[j][w+1];
	}
	sumM = sumM+gamma1[i][w+1]*gamma2[j][w+1];
      }
      alignM[i][j] = alignM[i][j]/sumM;
      if(p1 == p2){
	alignM[i][j] = log10(alignM[i][j]/(Pi[p1]*Pi[p2]))/log10(2.0);
      }else{
	alignM[i][j] = log10(alignM[i][j]/(2*Pi[p1]*Pi[p2]))/log10(2.0);
      }
    }
  }
  return alignM;
}

vector<vector<vector<double> > > readObserved(string obsFile){
  //initialize 3d array
  vector<vector<vector<double> > > obsMat;
  // resize it to (281 x 20 x 20)
  obsMat.resize(NUMSTATES);
  for (int i = 0; i < NUMSTATES; ++i) {
    obsMat[i].resize(20);
    for (int j = 0; j < 20; ++j)
      obsMat[i][j].resize(20);
  }
  //open the observed frequencies file
  ifstream obsF(obsFile.c_str());
  if(!obsF){
    cerr<<"Unable to read file: " << obsFile << endl;
    exit(-1);
  }
  //read in the observed frequencies
  for(int i = 0; i < NUMSTATES; ++i){
    for(int j = 0; j < 20; ++j){
      for(int k = 0; k < 20; ++k){
	obsF >> obsMat[i][j][k];
	obsMat[i][j][k] = obsMat[i][j][k]+(PSUEDO*PSUEDO);
      }
    }
  }
  obsF.close();
  //normalize the observed frequencies
  for(int i = 0; i < NUMSTATES; ++i){
    int sum = 0;
    for(int j = 0; j < 20; ++j){
      for(int k = 0; k < 20; ++k){
	sum = sum + obsMat[i][j][k];
      }
    }
    if(sum > 0){
      for(int j = 0; j < 20; ++j){
	for(int k = 0; k < 20; ++k){
	  obsMat[i][j][k] = obsMat[i][j][k]/sum;
	}
      }
    }
  }
  return obsMat;
}

vector<vector<double> >  readExpected(string expfile){
  //allocate expect matrix
  vector<vector<double> > expMat;
  //resize it to (281 x 20)
  expMat.resize(NUMSTATES);
  for(int i = 0; i < NUMSTATES; ++i){
    expMat[i].resize(20);
  }
  //open expected frequencies file
  ifstream expF(expfile.c_str());
  if(!expF){
    cerr << "Unable to read file: " << expfile << endl;
    exit(-1);
  }
  //read the expected freq values
  for(int i = 0; i < NUMSTATES; i++){
    for(int j = 0; j < 20; j++){
      //get the value from the file
      expF >> expMat[i][j];
      //add the psuedo counts
      expMat[i][j] = expMat[i][j]+PSUEDO;
    }
  }
  expF.close();
  //normalize the expected freq
  int sum;
  for(int i = 0; i< NUMSTATES; ++i){
    for(int j = 0; j < 20; ++j){
      sum = sum + expMat[i][j];
    }
  }
  for(int i = 0; i < NUMSTATES; ++i){
    for(int j = 0; j < 20; ++j){
      expMat[i][j] = expMat[i][j]/sum;
    }
  }
  return expMat;
}

vector<vector<double> > readGamma(string gcaName, int nres){
  //Initialize the gamma 2d array
  vector<vector<double> > gamma;
  //resize it to (nres x NUMNODES)
  gamma.resize(nres);
  for(int i = 0; i < nres; ++i){
    gamma[i].resize(NUMNODES);
  }
  
  //open the gamma file
  ifstream gcaF(gcaName.c_str());
  if(!gcaF){
    cerr << "Unable to read GCA file: "<< gcaF << endl;
    exit(-1);
  }
  int res = -1;
  string line;
  //parse loop
  while(res < nres && gcaF){
    getline(gcaF, line);
    if(line.substr(0, 7) == "GAMMA: ") {
      res++;
      //get string of float values
      string gcaValStr = line.substr(7, line.length()-7);
      //make a stream for the values
      istringstream in(gcaValStr);
      //now store the values in the gamma matrix
      for(int state = 0; state < NUMNODES; state++){
	//gamma[res][state]
	//in the gamma matrix, res are rows and states are columns
	in >> gamma[res][state];
	gamma[res][state] = gamma[res][state] + PSUEDO;
      }      
    }
  }
 
  gcaF.close();
  return gamma;
}

int seq2gca(string seqName){
  //generate the filename and put it in char array
  char *seqfile = (char*)malloc((seqName.length()+1)*sizeof(char*));
  strcpy(seqfile, seqName.c_str());
  long seqfLen = strlen(seqfile)+1;

  //Rest of seq2gca arguments 
  int nres;
  //Generates GCA file for this sequence
  seq2gca_(seqfile, &nres, seqfLen);
  free(seqfile);
  return nres;
}

vector<string> readFasta(string fileName){
  //parses a fasta file and stores the sequences in a vector
  //also splits the FASTA files into multiple seqeunce files for input into HMMSTR
  
  vector<string> Seqs;
  vector<string> SeqsNames;
  ifstream fin(fileName.c_str());
  if(!fin){
    cerr << "Could not read the FASTA file:" << fileName << endl;
    exit(-1);
  }
  string line;
  int counter = -1;
  getline(fin, line);
  //ignore any garbage before start of first sequence
  while(line[0]!= '>'){
    getline(fin, line);
  }
  //main parse loop
  while(fin){
    //start of a new sequences
    if(line[0] == '>'){
      line[0] = '_';
      //put the name of the sequence in the SeqsNames vector
      SeqsNames.push_back(line);
      //get ready to store the sequence
      Seqs.push_back("");
      counter++;
    }
    else{
      //keep appending the sequence
      Seqs[counter] += line;
    }
    //get the new line
    getline(fin, line);
  }
  fin.close();

  //write each sequence to a seperate .seq file
  for(int i = 0; i < Seqs.size(); i++){
    string outFileName = SeqsNames[i]+".seq";
    ofstream fout(outFileName.c_str());
    if(!fout){
      cerr << "Could not write SEQ file:" << outFileName << endl;
      exit(-1);
    }
    fout << ">" << SeqsNames[i];
    for(int j = 0; j < Seqs[i].length(); j++){
      //every 80 characters, print a new line (after the first one)
      if(j%80==0) fout << endl;
      fout << Seqs[i][j];
    }
    fout << endl;
    fout.close();
  }

  return SeqsNames;
}

string readSeq(string seqFile){
  string seq;

  ifstream seqF(seqFile.c_str());
  if(!seqF){
    cerr << "Unable to read file: "<< seqFile << endl;
    exit(-1); 
  }

  string line;
  getline(seqF, line);
  //ignore any garbage before start of first sequence
  while(line[0]!= '>'){
    getline(seqF, line);
  }
  //main parse loop
  while(seqF){
    //start of a new sequences
    if(line[0] == '>'){
      //get ready to store the string
      seq = "";
    }
    else{
      //keep appending the sequence
      seq += line;
    }
    //get the new line
    getline(seqF, line);
  }
  seqF.close();
  return seq;
}

void cleanupFiles(vector<string> seqNames){
  // delete the .seq and .gca file for each sequence  
  for(vector<string>::const_iterator it = seqNames.begin(); it != seqNames.end(); ++it){
    string gcaFile = *it+".gca";
    string seqFile = *it+".seq";
    if( remove(gcaFile.c_str()) != 0 ){
      cerr << "Unable to delete " << gcaFile << ". Please delete manually.\n";
    }
    if( remove(seqFile.c_str()) != 0 ){
      cerr << "Unable to delete " << seqFile << ". Please delete manually.\n";
    }
  }
  return;
}
