#include "position_of_read_in_reference.h"
#include "position_of_seed_in_reference.h"

// This function takes a read sequence, reference sequence (i.e. concatenated guide or guide_rc sequences), the read number, the seed length and threshold as input
// Threshold is maximum number of insertion+deletion+mutation that is allowed. If threshold=3, only 2 insertion+deletion+mutation is allowed.
// This function (i) identifies which guide is present in the read (ii) returns the true position of guide in reference sequence as output
int position_of_read_in_reference(std::string read, std::string reference, int read_number, int seed_len, int threshold=3){
			
	int legit_match_pos = -1;								// legit_match_pos is the location of the guide (that is present within a read) in the reference 
					
	for(int i = 0; i<read.length()/seed_len; ++i){			// We split the read into bins of size=seed length and check if any of the bins match with reference
		std::string seed = read.substr(i*seed_len,seed_len);// In this FOR loop, we only check the bins for matches. 
				
	    int match_pos = reference.find(seed);				// Once a match is found, we do base by base comparison in while loop below
		while (match_pos != std::string::npos){ 
			legit_match_pos = position_of_seed_in_reference(read, reference, read_number, i, match_pos, seed_len, threshold);					
			if (legit_match_pos == -1){
				match_pos=reference.find(seed,match_pos+1);
			}
			else if (legit_match_pos != -1){
				break;
			}								
		}
		
		if (legit_match_pos != -1){
			break;
		}
	}
	
	// Uncomment these 3 lines FIRST while troubleshooting.THESE ARE THE MOST IMPORTANT TROUBLESHOOTING HELP.
	std::cout << "reference:\t\t\t" << reference.substr(reference.rfind("X",legit_match_pos)+1,reference.find("X",legit_match_pos)-reference.rfind("X",legit_match_pos)) << "\n";
	std::cout << "read #" << read_number << ":\t\t\t" << read << "\n";
	std::cout << "legit match position:" << legit_match_pos << "\n\n";
	
	return legit_match_pos;
}	