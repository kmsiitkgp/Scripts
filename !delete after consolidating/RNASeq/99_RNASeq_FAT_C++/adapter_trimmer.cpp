#include "adapter_trimmer.h"


// This function takes a read sequence, adapter sequence, seed sequence as input
// The output is a trimmed read sequence
std::string adapter_trimmer(std::string seq, std::string adap_seq, std::string seed_seq){

  // Adapters are present in reads due to short insert size. For example, if cDNA insert is 50bp and we are doing 75bp
  // sequencing, there will be 25bp of adapter at 3' end. Ideally, adapter sequences are designed such that they are not
  // identical to sequences in genomes. In worst case, if adapter sequence is identical to a genome sequence, we might 
  // find multiple occurances of the adapter within a read. 
   
  // If we use find(), it will identify the 1st occurance from start of read i.e. 5' end. However, we are interested in
  // finding the 1st occurance from end of read i.e. 3' end as adapters are usually present at 3' end of reads. Hence, 
  // it is better to use rfind() which will identify the 1st occurance from end of read
  
  unsigned int mismatch = 0;
  unsigned int l = adap_seq.length();							// l is length of adapter sequence
  unsigned int s = seed_seq.length();							// s is length of seed sequence
  unsigned int j = 0;											// j tracks position of base from 3' end of read
  
  unsigned int adap_loc = seq.rfind(adap_seq); 					// identify 1st occurance of entire adapter from 3' end                 
    if (adap_loc < seq.length()){    							// if entire adapter is present in a read                  		
      seq = seq.substr(0,adap_loc);                   			// trim the read 
    }
    
	// Ideally, adapter sequences are designed such that they are not identical to sequences in genomes. The probability of a
	// 20bp adapter sequence present in genome is 1 in 4^20 i.e. 1 occurance every 1 trillion bp i.e. very unlikely. So, we can  
	// safely trim adapters that are 20bp long.
	else if ((adap_seq.length() > 20) && (seq.rfind(adap_seq.substr(0,20)) < seq.length())){	            		
      adap_loc = seq.rfind(adap_seq.substr(0,20));    			// identify 1st occurance of first 20 bases of adapter from 3' end	
      if (adap_loc < seq.length()){                     		// if partial adapter is present
        seq = seq.substr(0,adap_loc);                 			// trim the read
      }
	}
	
	// Now, let us consider the possibility that there might be partial adapters of varying length (less than 20 bp) present at 3' 
	// end of read. To remove these partial adapters, we check last base of adapter at 3' end to last base of read at 3' end. If 
	// there is match, we compare penultimate base of adapter to penultimate base of read and so on. If there is mistmatch, we check
    // if mismatched base in adapter matches with last base of read. WE DONT RESET the index of base in adapter but every time there
	// is mismatch we start comparisons again from last base of read. This way we can identify if partial adapter is present at 3' 
	// end of read.
	else {								              			// If entire adapter & partial adapter upto 20 bases is absent,      										
	  for (int i = adap_seq.length()-1; i>=0; --i){				// read the adapter sequence from 3' end
	    if (adap_seq[i] == seq[seq.length()-1-j]){				// if last base of adapter sequence matches with last base of read
		  j = j + 1;											// increase j so that next base of read can be compared to adapter
		}
		else {
		  i = i + j - 1;										// reset i to "penultimate base of previous match"..VERY IMPORTANT..see below
		  j = 0;												// if mismatch, again start comparison from last base of read	  
		  if (adap_seq[i] == seq[seq.length()-1-j]){			// VERY IMPORTANT..see below
          j = j + 1;											// if match, increase j so that next base of read can be compared to adapter
		  }														// if mismatch, compare next base of adapter to last base of read (at start of loop)	  		
		}
	  }
	  seq = seq.substr(0,seq.length()-j);						// if partial adapter is present, then j will be non-zero
	}
	

    /* Example to explain above code (focus from FOR loop iteration #13:
	Adap: A G A T C G G A A G A  G  C  A  C  A  C  G  T  C  T  G  A  A  C  T  C  C  A  G  T  C  A
		  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32	
	
	Read: AGGCGCTTTCTATAGCCTCCTTTACAATGCTGCTCACTTCATCAACAACAAATGCAGTCTCCTCGGAGGCCTGGAAGTCTTCCATCTTTAACGCGTCTCCTTCGCCCGCCGAGTAGAGTGGACTGCGCAGGCGCGG..AGATCGGAAG A G C A C
																																							   4 3 2 1 0	
	
	else {								              			      										
	  for (int i = adap_seq.length()-1; i>=0; --i){			
	    if (adap_seq[i] == seq[seq.length()-1-j]){	(LOOP a)			
		  j = j + 1;								(LOOP a)			
		}
		else {
		  i = i + j - 1;							(LOOP b)	
		  j = 0;									(LOOP b)			  
		  if (adap_seq[i] == seq[seq.length()-1-j]){(LOOP c)			
          j = j + 1;								(LOOP c)			
		  }															  		
		}
	  }
	  seq = seq.substr(0,seq.length()-j);						
	} 
	
    FOR loop (iteration)	i	j	adap_seq[i]		seq[seq.length()-1-j]	Loop executed	i	j	adap_seq[i]		seq[seq.length()-1-j]	Loop executed	i	j		
    1						32	0   A    			C						b				31  0	C               C						Multi			31  1   
	2						30	1	T				A						b				30	0	T 				C
	3						29	0	G				C						b				28	0	A				C 				
	4						27	0 	C   			C                       a				27	1					
	5						26	1	C 				A                       b               26	0	C 				C 						Multi  			26	1
	6						25  1	T 				A						b 				25	0	T 				C 
	7						24	0	C				C 						a 				24	1
	8						23	1	A 				A 						a 				23	2
	9						22	2	A 				C 						b 				23	0	A 				C  {very important to push i up to penultimate base of previous match i.e. A of AC}
	10						22	0	A 				C 						b 				21	0	G 				C
	11						20	0	T 				C 						b 				19	0	C 				C 						Multi			19	1
	12						18	1	T 				A 						b 				18	0	T 				C
	13						17	0	G 				C 						b 				16 	0	C 				C 						Multi 			16	1
	14						15	1	A 				A 						a 				15 	2	
	15						14	2	C 				C 						a 				14  3
    16						13  3   A 				G 						b 				15  0	A 				C 	{very important to push i up to penultimate base of previous match i.e. A of CAC}
	17						14	0	C 				C 						a 				14	1
    18						13	1	A 				A 						a 				13 	2
    19						12	2	C 				C 						a 				12	3
    20						11	3	G 				G 						a 				11	4
    21						10 	4	A 				A 						a 				10 	5	and so on..*/	

	// Finally, let us consider the possibility that there are sequencing or PCR errors in the adapter. Say, we want to allow upto 2
	// such mismatches i.e. all bases of an adapter are identical with a portion of read except for 2 bases and then trim the read.
	// We first identify a seed sequence within adapter that is unique. Then, search for this seed sequence from 3' end of read. 
	// If match is found, we check next base of adapter to next base of read and so on allowing for 2 mismatches. If more than 2 
	// mismatches are found, we stop analysis. If a match <2, we trim the read.
		  
	// Using seed sequence identify if adapter with 2 mismatch is present in read
	unsigned int seed_loc = seq.rfind(seed_seq); 					 
	if (seed_loc < seq.length()-seed_seq.length()-1){	//need to do at least 2 comparisons to check for 2 mismatches	
      for (int k = 0; k < l-s; ++k){
		if (seed_loc+s+k == seq.length()){				// break loop if end of read is reached
		  break;
		}
		if (adap_seq[s+k] != seq[seed_loc+s+k]){
		  mismatch = mismatch + 1;
		}
		else if (mismatch > 2){							//break loop is more than 2 mismatch
		  break;
		}
	  }
	}
	if ((seed_loc < seq.length()) && (mismatch <=2)){
	  seq = seq.substr(0, seed_loc);
	}   
		
  return seq;
}
