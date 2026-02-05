#include "quality_trimmer.h"

// This function takes a quality string, minimum quality score cutoff and base as input
// The output is a trimmed quality string
// Constant parameters (Eg: min_q_score, base) MUST BE DECLARED ONLY after ALL variable parameters (eg: qual) have been declared
// Else, c++ will give error during compilation.
// The quality trimmer will continue trimming until 5 consecutive bases have q_score >= min_q_score
std::string quality_trimmer(std::string qual, unsigned int min_q_score=20, unsigned int base=33){

  unsigned int base_count = 0;										// base_count tracks number of bases we analyzed from 3' end
  unsigned int good_base_count = 0;									// good_base_count tracks number of consecutive bases that are good
   
  for (unsigned int i = qual.length()-1; i>0; --i){					// read the quality string from 3' end
    base_count = base_count + 1; 									// for every base analyzed, increase base_count by 1
    int ord = qual[i];												// find ord value of qual[i]
    if (ord-base < min_q_score){									// check if base is bad i.e. ord-base < min_q_score
      good_base_count = 0;											// if bad base, reset consecutive number of good bases to 0
	}
    else{															// check if base is good i.e. ord-base >= min_q_score
      good_base_count = good_base_count + 1;						// if good base, increase consecutive number of good bases by 1
    }
	if (good_base_count == 5){										// if 5 consecutive good bases are found, stop analysis
      break;
	}
  }
  qual = qual.substr(0,(qual.length()-base_count+good_base_count));	// trim the quality string using count and good_base_count
  return qual;														// return quality string
}