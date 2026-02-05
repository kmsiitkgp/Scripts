#include "position_of_seed_in_reference.h"
	
// This function takes a read sequence, reference sequence (i.e. concatenated guide or guide_rc sequences), the read number, the bin number i, 
// the initial matching postion of seed in reference i.e. match_pos, the seed length and threshold as input
// Threshold is maximum number of insertion+deletion+mutation that is allowed. If threshold=3, only 2 insertion+deletion+mutation is allowed.
// This function returns the position of guide in reference sequence as output
int position_of_seed_in_reference(std::string read, std::string reference, int read_number, int i, int match_pos, int seed_len, int threshold=3){
	
	int potential_match_pos = -1;
	
	// If any of the bin matches with reference, we check if next base in read & reference match. 
	// If they dont match, we figure out if it is a mismatch or insertion or deletion by comparing the next 5 bases
	
	int j = 0;			// j keeps track of number of bases that have been compared downstream of match position
	int mis_down = 0;	// mismatch tracks number of mismatches between read and reference downstream of match position
	int ins_down = 0;	// insertion tracks number of insertions between read and reference downstream of match position
	int del_down = 0;	// deletion tracks number of deletions between read and reference downstream of match position
	
	//std::cout << "read:\t\t\t" << read << "\n";
	//std::cout << "seed:\t\t\t " << seed << "\n";
	//std::cout << "reference:\t\t\t" << reference.substr(reference.rfind("X",match_pos)+1,reference.find("X",match_pos)-reference.rfind("X",match_pos)) << "\n";

	while ((mis_down+ins_down+del_down < threshold) && (reference[match_pos+seed_len+j] != 'X')){
		//std::cout << "read[" << i*seed_len+seed_len+j+ins_down-del_down << "]:" << read[i*seed_len+seed_len+j+ins_down-del_down] << "\n";
		//std::cout << "reference[" << match_pos+seed_len+j << "]:" << reference[match_pos+seed_len+j] << "\n";
		
		// if end of read is NOT reached before end of guide,
		if (i*seed_len+seed_len+j+ins_down-del_down < read.length()){	
			//if nulceotide in read doesnt match with nucleotide in reference, check if mismatch or insertion or deletion
			if (read[i*seed_len+seed_len+j+ins_down-del_down] != reference[match_pos+seed_len+j]){
				
				int mis=0;						
				int ins=0;
				int del=0;
				int l=1;	// l keeps track of number of bases that have been compared downstream of mismatch position
				while (l<=5 && (reference[match_pos+seed_len+j+l] != 'X') && (i*seed_len+seed_len+j+ins_down-del_down+l < read.length())){
					if (read[i*seed_len+seed_len+j+ins_down-del_down+l] == reference[match_pos+seed_len+j+l]){
						mis=mis+1;
					}
							
					if (read[i*seed_len+seed_len+j+ins_down-del_down+l] == reference[match_pos+seed_len+j+l-1]){
						ins=ins+1;
					}
					
					if (read[i*seed_len+seed_len+j+ins_down-del_down+l-1] == reference[match_pos+seed_len+j+l]){
						del=del+1;
					}
					l=l+1;
				}
					
				if (mis >= ins && mis >= del && l-mis <= threshold){
					mis_down=mis_down+1;
				}						
				else if (ins > mis && ins > del && l-ins <= threshold){
					ins_down=ins_down+1;
				}
				else if (del > mis && del > ins && l-del <= threshold){
					del_down=del_down+1;
				}
				else{
					mis_down=mis_down+1;
				}
				//std::cout << "mis" << mis << "\t ins" << ins << "\t del" << del << "\n"; 	
				//std::cout << "mismatch:" << mis_down << "\t insertion:" << ins_down << "\t deletion:" << del_down << "\n";
			}
		}
			
		// if end of read is reached before end of guide, increase mismatch
		else if (i*seed_len+seed_len+j+ins_down-del_down >= read.length()){
			//std::cout << "index location:" << i*seed_len+seed_len+j+ins_down-del_down << "\n";
			//std::cout << "read len:"<< read.length() << "\n";
			mis_down = mis_down + 1;
		}
		j=j+1;
	}
	
	int k = 0;				// k keeps track of number of bases that have been compared upstream of match position
	int mis_up = 0;			// mismatch tracks number of mismatches between read and reference upstream of match position
	int ins_up = 0;			// insertion tracks number of insertions between read and reference upstream of match position
	int del_up = 0;			// deletion tracks number of deletions between read and reference upstream of match position
	
	while ((mis_down+mis_up+ins_down+ins_up+del_down+del_up < threshold) && (reference[match_pos-k-1] != 'X')){				
		//std::cout << "read[" << i*seed_len-k-1-ins_up+del_up << "]:" << read[i*seed_len-k-1-ins_up+del_up] << "\n";
		//std::cout << "reference[" << match_pos-k-1 << "]:" << reference[match_pos-k-1] << "\n";
		//std::cout << "i:" << i << "\t k:" << k << "\t insertion_up" << ins_up << "\t deletion_up" << del_up << "\n;
		
		// if start of read is NOT reached before start of guide,
		if (i*seed_len-k-1-ins_up+del_up >= 0){
			//if nulceotide in read doesnt match with nucleotide in reference, check if mismatch or insertion or deletion
			if (read[i*seed_len-k-1-ins_up+del_up] != reference[match_pos-k-1]){
				
				int mis=0;
				int ins=0;
				int del=0;
				int l_up=0;						
				
				while (l_up<5 && (reference[match_pos-k-1-l_up] != 'X') && (i*seed_len-k-1-ins_up+del_up-l_up >=0)){
					if (read[i*seed_len-k-1-ins_up+del_up-l_up] == reference[match_pos-k-1-l_up]){
						mis=mis+1;
					}				
				
					if (read[i*seed_len-k-1-ins_up+del_up-l_up] == reference[match_pos-k-l_up]){
						ins=ins+1;
					}						
				
					if (read[i*seed_len-k-1-ins_up+del_up-l_up] == reference[match_pos-k-2-l_up]){
						del=del+1;
					}
					l_up=l_up+1;
				}
					
				if (mis >= ins && mis >= del && l_up-mis <= threshold){
					mis_up=mis_up+1;
				}						
				else if (ins > mis && ins > del && l_up-ins <= threshold){
					ins_up=ins_up+1;
				}
				else if (del > mis && del > ins && l_up-del <= threshold){
					del_up=del_up+1;
				}
				else{
					mis_up=mis_up+1;
				}
				//std::cout << "match_mis" << mis << "\t match_ins" << ins << "\t match_del" << del << "\n";
				//std::cout << "mismatch:" << mis_up << "\t insertion_up:" << ins_up << "\t deletion_up:" << del_up << "\n";
			}
		}
		
		// if start of read is reached before start of guide, increase mismatch
		else if (i*seed_len-k-1-ins_up+del_up < 0){
			//std::cout << "index location:" << i-k-1-ins_up+del_up << "\n"; 
			mis_up = mis_up + 1;
		}
		k=k+1;
	}
	
	if ((j+k>2) && (mis_down + mis_up + ins_down + ins_up + del_down + del_up < threshold)){
		potential_match_pos = match_pos;
	}
	
	return potential_match_pos;
}
	