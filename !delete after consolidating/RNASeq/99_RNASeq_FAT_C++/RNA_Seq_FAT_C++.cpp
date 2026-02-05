#include "quality_trimmer.h"
#include "adapter_trimmer.h"
#include "single_end.h"
#include "paired_end.h"
#include <iostream>
#include <fstream>          // you need this for opening, reading, writing file
#include <string>           // you need this for string operations
#include <algorithm>        // you need this for find()
#include <chrono>           // you need this to calculate time taken by the script


/* This script does adapter, polyG and quality trimming in paired and single end files.
Input parameters for the script in specific order:
(1) type of read (SE or PE)
(2) minimum length of read after trimming (preferably above 36)
(3) minimum quality score of base to trim from read (preferably 20)
(4) base (usually 33, rarely 64)
(5) path to output directory
(6) forward read fastq file
(7) adapter to trim from forward read
(8) reverse read fastq file (if single end, this parameter is skipped automatically)
(9) adapter to trim from reverse read (if single end, this parameter is skipped automatically)

This script will remove entire adapters as well as partial adapters longer than 10 bases from reads.
Output is a fastq.gz file that can be used for STAR etc.

Note: NextSeq & NovaSeq use 2-colour chemistry for base detection. Green=T, Red=C, Green+Red=A, no signal=G.
2 channel system can't distinguish "G" and "no signal" because both result in no channel emission.
PolyG is an usual artifact of 2-channel systems. 20 or more consecutive G's are considered as artifact and removed.*/

using namespace std;
using namespace std::chrono;

int main (int argc, char** argv) {

  string t = argv[1];
  unsigned int l = atoi(argv[2]);          								// store minimum length of trimmed read in variable l
  unsigned int qscore = atoi(argv[3]);     								// store Quality score threshold in variable qscore
  unsigned int b = atoi(argv[4]);          								// store Phred base in variable b
  string o = argv[5];               									// store output directory where trimmed fastq files will be created in variable o
  string f = argv[6];              										// store input forward read filename with path in variable f
  string f_adap = argv[7];
  
  for (int i = f.length()-1; i >=0; --i){								// extract forward file name from file path
    if (f[i] == '/'){
      f = f.substr(i+1,f.length());
      break;
    }
  }
  string o1 = o + "/trimmed_P_"+ f;
  string o2 = o + "/UP_"+ f;
    
  //cout << "\nSE or PE:" << t 
  //	 << "\nMin Len:" << l 
  //	 << "\nMin QScore:" << qscore 
  //	 << "\nBase:" << b 
  //	 << "\nOutput Directory:" << o 
  //	 << "\nForward Read File:" << f 
  //	 << "\nForward Adapter:" << f_adap 
  //	 << "\nTrimmed Forward File:" << o1 
  //	 << "\nQC Failed Untrimmed Forward File:" << o2;
  
  // As you notice string f = argv[6] converts file path to filename. So, it is best to open files using argv[] and then pass them to se_adapter_trimmer??double check this
  // instead of passing the file names to se_adapter_trimmer and opening the files within se_adapter_trimmer 
  
  std::ifstream f_file;
  f_file.open(argv[6],ios::in);
  std::ofstream p_f_file;
  p_f_file.open(o1,ios::out);
  std::ofstream up_f_file;
  up_f_file.open(o2,ios::out);
  if (f_file.is_open()){cout << "Hello1\n";} else { cout << "Bye1\n";}
  if (p_f_file.is_open()){cout << "Hello2\n";} else { cout << "Bye2\n";}
  if (up_f_file.is_open()){cout << "Hello3\n";} else { cout << "Bye3\n";}
  
  if (t == "PE"){
	string r = argv[8];             					// store input reverse read filename with path in variable r
	string r_adap = argv[9];
	
	for (int i = r.length()-1; i >=0; --i){				// extract reverse file name from file path
      if (r[i] == '/'){
        r = r.substr(i+1,r.length());
        break;
      }
    }
    string o3 = o + "/trimmed_P_"+ r;
    string o4 = o + "/UP_"+ r;
	
	// cout << "\nReverse Read File:" << r 
  	// << "\nReverse Adapter:" << r_adap 
  	// << "\nTrimmed Reverse File:" << o3 
  	// << "\nQC Failed Untrimmed Reverse File:" << o4;
	
	std::ifstream r_file;
    r_file.open(argv[8],ios::in);
    std::ofstream p_r_file;
    p_r_file.open(o3,ios::out);
    std::ofstream up_r_file;
    up_r_file.open(o4,ios::out);
	
	pe_adapter_trimmer(l, qscore, b, f_file, f_adap, p_f_file, up_f_file, r_file, r_adap, p_r_file, up_r_file);
  }
  else if (t == "SE"){
	se_adapter_trimmer(l, qscore, b, f_file, f_adap, p_f_file, up_f_file);
  }
  else{
	printf("Please enter the correct parameters");
  } 
  return 0;
}