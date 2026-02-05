#include "position_of_seed_in_reference.h"
	
// This function takes a read sequence, reference sequence (i.e. concatenated guide or guide_rc sequences), the read number, the bin number i, 
// the initial matching postion of seed in reference i.e. match_pos, the seed length and threshold as input
// Threshold is maximum number of insertion + deletion + mutation that is allowed. If threshold = 3, only 2 insertion + deletion + mutation is allowed.
// This function returns the position of guide in reference sequence as output

// I have noticed that this function is sometimes not called by find_orientation() when using cutadapt results but works fine on full fq.gz file.
// The find_orientation() calls this function for each read twice:
// once on reference_s and once on reference_rc. Then, it calls them again for 2nd read ans so on.
// For example, if I run find_orientation using both reference_s and reference_rc, only 1st read is run on both reference_s and reference_rc.
// The find_orientation() is unable to call this function for the 2nd read on reference_s.
// Similarly, if I run find_orientation() using only reference_s, then only the first 9 reads are aligned. find_orientation() is unable to call
// this function on the 10th read.
// I have observed that this happens when the seed aligns towards the end of a guide and only a few bases are remaining on the guide. While no 
// error messages are being reported, I suspect this somehow breaks the position_of_seed_in_reference() and hence find_orientation() is unable
// to call it again. 
//Also, when read length is 16 just 1 base more than min match this happens.

std::string position_of_seed_in_reference(std::string read, std::string reference, int read_n, int seed_len, int max_errors, int min_match, std::string troubleshoot_mode){
	
	std::string result;				// this stores all the SAM paramters in comma separated format and returns as output
	std::string initial_seed;		// initial seed sequence
	int initial_seed_pos_ref;		// position of initial seed in reference; output from find()	
	int start_site_ref;				// position in reference starting from which seed is searched; used as input in find() to search for seed beyond 1st match position
	int start_site_read = 0;		// position in read starting from which the next seed is generated; used as input in substr()
	
	int n_match;					// number of matches in current alignment
	int n_insertions;				// number of insertions in current alignment
	int n_deletions;				// number of deletions in current alignment
	int n_mismatch;					// number of mismatches in current alignment
	int n_errors;					// number of errors in current alignment = n_insertions + n_deletions + n_mismatch
	int n_tandem_match;				// maximum number of consecutive matches in current alignment before error
	
	int n_read;						// nth base in read after either downstream of (start_site_read + seed_len) or upstream of (start_site_read) being compared in base-by-base comparisons
	int n_ref;						// nth base in reference either downstream of (initial_seed_pos_ref + seed_len) or upstream of (initial_seed_pos_ref) being compared in base-by-base comparisons
	int current_base_pos_read_down;	// position of current base in read being compared downstream of alignment = start_site_read + seed_len + n_read
	int current_base_pos_read_up;	// position of current base in read being compared upstream of alignment = start_site_read-n_read 
	int current_base_pos_ref_down;	// position of current base in reference being compared downstream of alignment = initial_seed_pos_ref + seed_len + n_ref
	int current_base_pos_ref_up;	// position of current base in reference being compared upstream of alignment = initial_seed_pos_ref-n_ref
	int base_pos_ref_X_down; 		// number of bases between 1st base on reference that seed matches and nearest X i.e. end of guide
	
	int best_pos_ref_X_down = -1;		// end position of guide in reference corresponding to best alignment
	int best_pos_ref_X_up = -1;		// start position of guide in reference corresponding to best alignment
	int best_pos_ref = -1;			// start position of best alignment in reference
	int best_pos_read = 0;			// start position of best alignment in read
	int best_match = -1;				// match score of best alignment
	int best_match_second = -2;		// match score of 2nd best alignment
	int best_ins = 0;				// number of insertions in best alignment
	int best_del = 0;				// number of deletions in best alignment
	int best_mis = 0;				// number of mismatches in best alignment
	int best_error = 0;				// number of errors in best alignment
	int best_tandem_match = 0;		// number of tanden matches in best alignment
	
	std::string QNAME;				// variable to store the name of the current read
	int FLAG;						// variable to store the FLAG values of the current read
	int POS;						// variable to store the mapping start position of the current read in reference
	int MAPQ;						// variable to store the MAPQ values of the current read
	std::string CIGAR;				// variable to store the "modified" CIGAR alignment values of the current read
	std::string RNEXT;				// variable to store the RNEXT values of the current read
	int PNEXT;						// variable to store the PNEXT values of the current read
	int TLEN;						// variable to store the TLEN values of the current read
	std::string SEQ;				// variable to store the sequence of the current read	
	std::string QUAL;				// variable to store the quality scores of the current read	
	std::string RSEQ;				// variable to store the sequence of the current read
	
	//std::cout << "Read Number\t\t:" << read_n << "\n";
	
	// We keep finding a new seed until end of read is reached
	while ((start_site_read + min_match) <= read.length()){
		
		// We use first seed_len bases of a read as initial_seed
		initial_seed = read.substr(start_site_read,seed_len);
		
		// Important to reset start_site_ref to 0 for every new seed
		start_site_ref = 0;	
		
		// We keep finding a new alignment until end of reference in reached
		while ((start_site_ref + min_match) <= reference.length()){
			
			// Find matching position of initial_seed in reference
			// NOTE: find() only returns the 1st match of inital seed. However, there could be additional matches beyond the 1st match. 
			// So, we repeat the search for initial seed again beyond the first match position using while loop by adjusting start_site_ref in each iteration
			initial_seed_pos_ref = reference.find(initial_seed, start_site_ref);
			
			// If a match for the current seed is found, start comparing base by base and find its alignment score
			// At the end of comparison, reset start_site_ref so current seed is now searched for matches beyond the 1st matching position
			if (initial_seed_pos_ref != -1){
				
				// Since a match is found, increase match score to number of matching bases = seed length
				// Set all other values as 0 for each iteration i.e. each new alignment
				n_match = seed_len;
				n_read = 0;
				n_ref = 0;
				n_insertions = 0;
				n_deletions = 0;
				n_mismatch = 0;
				n_errors = 0;
				n_tandem_match = 0;				
				current_base_pos_read_down = start_site_read + seed_len + n_read;
				current_base_pos_ref_down = initial_seed_pos_ref + seed_len + n_ref;
				base_pos_ref_X_down = reference.find('X', initial_seed_pos_ref);
				
				// Continue base by base comparison 				
				// (i) end of guide is reached 
				// (ii) end of read is reached									
				
				// ALIGNING DOWNSTREAM OF THE SEED IN GUIDE i.e REFERENCE
				while (reference[current_base_pos_ref_down] != 'X' && current_base_pos_read_down < read.length()){
										
					// If the bases match, increase match score
					if (reference[current_base_pos_ref_down] == read[current_base_pos_read_down]){
						/*if(read_n == 818){
							std::cout << reference[current_base_pos_ref_down] << "\n";
							std::cout << read[current_base_pos_read_down] <<"\n";
						}*/
						n_match = n_match + 1;
						n_read = n_read + 1;
						n_ref = n_ref + 1;						
						current_base_pos_read_down = start_site_read + seed_len + n_read;
						current_base_pos_ref_down = initial_seed_pos_ref + seed_len + n_ref;
						if (n_errors == 0){
							n_tandem_match = n_read + seed_len;							
						}	
					}
					
					// When an error is reached, we need to check for deletion, insertion as well as mismatch and then decide type of error.					
					// If we first check for deletion and DO NOT check for insertion or mismatch, we will have incorrect alignment score.
					
					// Read with REAL insertion							:TTTTTACACACTTT
					//													 *****X********			
					// Reference										:TTTTT CACACTTT
					
					// Assume we first check for deletion and DO NOT check for insertion or mismatch
					// Read with REAL insertion checked for deletion	:TTTTT ACACACTTT
					//													 *****X****XX***
					// Reference										:TTTTTCACAC  TTT
					// Hence, the algorithm will classify it as deletion and the alignment will have an incorrect score.					
					
					// If there is error, it could be due to insertion, deletion, mismatch in read.
					// We check for all 3 types of error until the next error
					// Each instance of else loop will account for ONLY 1 instance of insertion or deletion or mismatch
					else {						
						int deletion_base_pos_ref = current_base_pos_ref_down + 1;
						int deletion_base_pos_read = current_base_pos_read_down;
						int deletion_match_down = 0;
						
						int insertion_base_pos_ref = current_base_pos_ref_down;
						int insertion_base_pos_read = current_base_pos_read_down + 1;
						int insertion_match_down = 0;
						
						int mismatch_base_pos_ref = current_base_pos_ref_down + 1;					
						int mismatch_base_pos_read = current_base_pos_read_down + 1;
						int mismatch_match_down = 0;
						
						// Checking for deletion
						while (reference[deletion_base_pos_ref] == read[deletion_base_pos_read]){
							/*if (read_n == 818){
								std::cout <<"Deletion" << "\n";
								std::cout << reference[deletion_base_pos_ref] << "\n";
								std::cout << read[deletion_base_pos_read] <<"\n";
							}*/
							deletion_match_down = deletion_match_down + 1;
							deletion_base_pos_ref = deletion_base_pos_ref + 1;
							deletion_base_pos_read = deletion_base_pos_read + 1;							
						}		

						// Checking for insertion
						while (reference[insertion_base_pos_ref] == read[insertion_base_pos_read]){
							/*if (read_n == 818){
								std::cout <<"Insertion" << "\n";
								std::cout << reference[insertion_base_pos_ref] << "\n";
								std::cout << read[insertion_base_pos_read] <<"\n";
							}*/
							insertion_match_down = insertion_match_down + 1;
							insertion_base_pos_ref = insertion_base_pos_ref + 1;
							insertion_base_pos_read = insertion_base_pos_read + 1;							
						}

						// Checking for mismatch
						while (reference[mismatch_base_pos_ref] == read[mismatch_base_pos_read]){
							/*if(read_n == 818){
								std::cout <<"Mismatch" << "\n";
								std::cout << reference[mismatch_base_pos_ref] << "\n";
								std::cout << read[mismatch_base_pos_read] <<"\n";
							}*/
							mismatch_match_down = mismatch_match_down + 1;
							mismatch_base_pos_ref = mismatch_base_pos_ref + 1;
							mismatch_base_pos_read = mismatch_base_pos_read + 1;
						}	
							
						// Adjust n_match, n_read and n_ref based on type of error
						// First preference given to mismatch			
						if (mismatch_match_down >= std::max(deletion_match_down, insertion_match_down)){
							n_mismatch = n_mismatch + 1;
							n_match = n_match + mismatch_match_down;
							n_read = n_read + mismatch_match_down + 1;
							n_ref = n_ref + mismatch_match_down + 1;								
						}
						// Second preference given to deletion over insertion
						else if (deletion_match_down >= std::max(insertion_match_down, mismatch_match_down)){
							n_deletions = n_deletions + 1;
							n_match = n_match + deletion_match_down;
							n_read = n_read + deletion_match_down;
							n_ref = n_ref + deletion_match_down + 1;						
						}
						// Third preference given to insertion 						
						else if (insertion_match_down > std::max(deletion_match_down, mismatch_match_down)){
							n_insertions = n_insertions + 1;
							n_match = n_match + insertion_match_down;
							n_read = n_read + insertion_match_down + 1;
							n_ref = n_ref + insertion_match_down;						
						} 					
						
						// Update current_base_pos_read_down, current_base_pos_ref_down and n_errors
						current_base_pos_read_down = start_site_read + seed_len + n_read;
						current_base_pos_ref_down = initial_seed_pos_ref + seed_len + n_ref;
						n_errors = n_insertions + n_deletions + n_mismatch;
					}
					
					// If the total errors already crossed max errors, skip upstream alignment
					if (n_errors > max_errors){
						goto END;
					}
				}		
				
				// ALIGNING UPSTREAM OF THE SEED IN GUIDE i.e REFERENCE
				// So far we only aligned downstream of seed. There could be bases upstream of seed in the guide.
				// We need to check if these bases also match with the read and if not increase n_errors.
				// Reset upstream_base_ref and upstream_base_read to 1 for every alignment
				
				n_read = 1;
				n_ref = 1;		
				current_base_pos_read_up = start_site_read - n_read;
				current_base_pos_ref_up = initial_seed_pos_ref - n_ref;
				//std::cout << "current_base_pos_ref_up:" << current_base_pos_ref_up << "\n";
				//std::cout << "current_base_pos_ref_up:" << reference[current_base_pos_ref_up] << "\n";
				
				// NOTE: If n_errors already = 2 and there are bases upstream in guide (OR) 
				// if start fo read is already reached and there are bases upstream in guide, 
				// we still need to enter loop until there are bases upstream in guide/
				// Hence, we do not set n_errors limit or read pos limit in this while loop
				while(reference[current_base_pos_ref_up] != 'X'){			

					// If there are bases upstream in read, do comparison
					if (current_base_pos_read_up >= 0){						
						
						// If the bases match, increase match score
						if (reference[current_base_pos_ref_up] == read[current_base_pos_read_up]){
							/*if(read_n == 818){
								std::cout << reference[current_base_pos_ref_up] << "\n";
								std::cout << read[current_base_pos_read_up] <<"\n";
							}*/
							n_match = n_match + 1;
							n_read = n_read + 1;
							n_ref = n_ref + 1;
							current_base_pos_read_up = start_site_read - n_read;
							current_base_pos_ref_up = initial_seed_pos_ref - n_ref;						
						}
						//else if (reference[current_base_pos_ref_up] != 'X'){
						else{
							int deletion_base_pos_ref = current_base_pos_ref_up - 1;
							int deletion_base_pos_read = current_base_pos_read_up;
							int deletion_match_up = 0;
							
							int insertion_base_pos_ref = current_base_pos_ref_up;
							int insertion_base_pos_read = current_base_pos_read_up - 1;
							int insertion_match_up = 0;
							
							int mismatch_base_pos_ref = current_base_pos_ref_up - 1;					
							int mismatch_base_pos_read = current_base_pos_read_up - 1;
							int mismatch_match_up = 0;	
						
							// Checking for deletion
							while(deletion_base_pos_read >= 0 && deletion_base_pos_ref >= 0 && reference[deletion_base_pos_ref] != 'X'){
								if (reference[deletion_base_pos_ref] == read[deletion_base_pos_read]){
									/*if(read_n == 818){
										std::cout <<"Deletion_up" << "\n";
										std::cout << reference[deletion_base_pos_ref] << "\n";
										std::cout << read[deletion_base_pos_read] <<"\n";
									}*/
									deletion_match_up = deletion_match_up + 1;
									deletion_base_pos_ref = deletion_base_pos_ref - 1;
									deletion_base_pos_read = deletion_base_pos_read - 1;
								}
								else{								
									break;
								}
							}

							// Checking for insertion
							while(insertion_base_pos_read >= 0 && insertion_base_pos_ref >= 0 && reference[insertion_base_pos_ref] != 'X'){
								if (reference[insertion_base_pos_ref] == read[insertion_base_pos_read]){
									/*if(read_n == 818){
										std::cout <<"Insertion_up" << "\n";
										std::cout << reference[insertion_base_pos_ref] << "\n";
										std::cout << read[insertion_base_pos_read] <<"\n";
									}*/
									insertion_match_up = insertion_match_up + 1;
									insertion_base_pos_ref = insertion_base_pos_ref - 1;
									insertion_base_pos_read = insertion_base_pos_read - 1;														
								}
								else{
									break;
								}
							}

							// Checking for mismatch
							while(mismatch_base_pos_read >= 0 && mismatch_base_pos_ref >= 0 && reference[mismatch_base_pos_ref] != 'X'){
								if (reference[mismatch_base_pos_ref] == read[mismatch_base_pos_read]){
									/*if(read_n == 818){
										std::cout <<"Mismatch_up" << "\n";
										std::cout << reference[mismatch_base_pos_ref] << "\n";
										std::cout << read[mismatch_base_pos_read] <<"\n";
									}*/
									mismatch_match_up = mismatch_match_up + 1;
									mismatch_base_pos_ref = mismatch_base_pos_ref - 1;
									mismatch_base_pos_read = mismatch_base_pos_read - 1;							
								}
								else{
									break;
								}
							}
							
							// Adjust n_match, n_read and n_ref based on type of error
							// First preference given to mismatch			
							if (mismatch_match_up >= std::max(deletion_match_up, insertion_match_up)){
								n_mismatch = n_mismatch + 1;
								n_match = n_match + mismatch_match_up;
								n_read = n_read + mismatch_match_up + 1;
								n_ref = n_ref + mismatch_match_up + 1;								
							}
							// Second preference given to deletion over insertion
							else if (deletion_match_up >= std::max(insertion_match_up, mismatch_match_up)){
								n_deletions = n_deletions + 1;
								n_match = n_match + deletion_match_up;
								n_read = n_read + deletion_match_up;
								n_ref = n_ref + deletion_match_up + 1;						
							} 
							// Third preference given to insertion							
							else if (insertion_match_up > std::max(deletion_match_up, mismatch_match_up)){
								n_insertions = n_insertions + 1;
								n_match = n_match + insertion_match_up;
								n_read = n_read + insertion_match_up + 1;
								n_ref = n_ref + insertion_match_up;						
							}
							
							// Update current_base_pos_read_up, current_base_pos_ref_up and n_errors
							current_base_pos_read_up = start_site_read - n_read;
							current_base_pos_ref_up = initial_seed_pos_ref - n_ref;
							n_errors = n_insertions + n_deletions + n_mismatch;							
						}
					}
					// If there are NO bases upstream in read, increase error
					else{
						n_errors = n_errors + 1;
						current_base_pos_ref_up = current_base_pos_ref_up - 1;	
					}

					// If the total errors already crossed max errors, skip upstream alignment
					if (n_errors > max_errors){
						goto END;
					}
				}

				/*std::cout << "\nInitial Seed:" << initial_seed << "\n";				
				std::cout << "n match:" << n_match << "\n";					
				std::cout << "Best match:" << best_match << "\n";
				std::cout << "Best match 2nd:" << best_match_second << "\n";
				std::cout << "n tandem match:" << n_tandem_match << "\n";
				std::cout << "Best tandem match:" << best_tandem_match << "\n";
				std::cout << "Best tandem match 2nd:" << n_tandem_match << "\n";
				std::cout << "n insertion:" << n_insertions << "\n";
				std::cout << "Best insertion:" << best_ins << "\n";
				std::cout << "n deletion:" << n_deletions << "\n";
				std::cout << "Best deletion:" << best_del << "\n";
				std::cout << "n mismatch:" << n_mismatch << "\n";
				std::cout << "Best mismatch:" << best_mis << "\n";
				std::cout << "n errors:" << n_errors << "\n";
				std::cout << "Best error:" << best_error << "\n";
				std::cout << "Current seed pos ref:" << initial_seed_pos_ref << "\n";
				std::cout << "Best pos ref:" << best_pos_ref << "\n";*/
				
				END:
				// CHECKING QUALITY OF ALIGNMENT	
				// If alignment using current seed has match score greater than min match score and also has the best score for current read, save it 
				if (n_match > best_match && n_match >= min_match && n_errors <= max_errors){
					best_match = n_match;
					best_ins = n_insertions;
					best_del = n_deletions;
					best_mis = n_mismatch;
					best_error = n_errors;
					best_tandem_match = n_tandem_match;
					best_pos_ref = initial_seed_pos_ref;
					best_pos_read = start_site_read;
					best_pos_ref_X_down = base_pos_ref_X_down;
				}
				else if ((n_match == best_match) && (n_insertions + n_deletions == best_ins + best_del) && (n_mismatch == best_mis) && (n_errors == best_error) && (n_tandem_match == best_tandem_match) && (base_pos_ref_X_down != best_pos_ref_X_down)){
					best_match_second = n_match;					
				} 
				else if ((n_match >= best_match) && (n_insertions + n_deletions <= best_ins + best_del) && (n_mismatch <= best_mis) && (n_errors <= best_error)&& (n_tandem_match >= best_tandem_match)){
					best_match = n_match;
					best_ins = n_insertions;
					best_del = n_deletions;
					best_mis = n_mismatch;
					best_error = n_errors;
					best_tandem_match = n_tandem_match;
					best_pos_ref = initial_seed_pos_ref;
					best_pos_read = start_site_read;
					best_pos_ref_X_down = base_pos_ref_X_down;
				} 		
				
				// Reset start_site_ref so seed is now searched for matches beyond the 1st matching position
				start_site_ref = initial_seed_pos_ref + seed_len;
				//std::cout << "start_Site_ref:" << start_site_ref << "\n";
			}
			
			// If a match for the current seed is NOT found, break out of loop and use new seed
			else{
				break;
			}
		}
		
		// Use new seed and to see if there are better alignments
		if (best_match >= min_match && best_error <= max_errors && best_tandem_match > 0){
			start_site_read = start_site_read + best_tandem_match;		
		}
		else{
			start_site_read = start_site_read + 1;
		}	
	}
	
	// After exhausting all possible seed sequences,
	// If NO proper alignment is found, mark read as unmapped
	if(best_match == -1){
		best_pos_ref = -1;
		//SEQ = read.substr(start_site_read,min_match);					// If proper alignment is NOT found, return read sequence upto min_match positions
		SEQ = read;
		RSEQ = "*";
		CIGAR = "*";
	}
	
	// If MORE THAN ONE proper alignment is found, mark read as multimapped
	else if(best_match == best_match_second){
		best_pos_ref = -10;
		//SEQ = read.substr(start_site_read,min_match);					// If proper alignment is NOT found, return read sequence upto min_match positions
		SEQ = read;
		RSEQ = "*";
		CIGAR = "*";		
	}	 
	
	// If A SINGLE BEST alignment is found, mark read as mapped
	else if (best_pos_ref > 0){												
		//SEQ = read.substr(best_pos_read,(best_match + n_errors));		// If proper alignment is found, return read sequence upto (best_match + n_errors) positions
		SEQ = read;	
		best_pos_ref_X_up = reference.rfind('X', best_pos_ref);		
		RSEQ = reference.substr(best_pos_ref_X_up + 1,(best_pos_ref_X_down - best_pos_ref_X_up));	// If proper alignment is found, return reference sequence
		CIGAR = std::to_string(best_match) + 'M' + std::to_string(best_ins) + 'I' + std::to_string(best_del) + 'D' + std::to_string(best_mis) + 'm';	
	}
		
	QNAME = "Read_" + std::to_string(read_n);		// For simplicity, we assign read number as QNAME
	FLAG = 0;												// Not yet developed to assign FLAG values. Currently, all set to 0	
	POS = best_pos_ref;	
	MAPQ = 0;												// Not yet developed to assign MAPQ values. Currently, all set to 0
	RNEXT = "*";											// Since we ONLY handle single end reads, there is no mate pair. Hence, set as *
	PNEXT = 0;											// Since we ONLY handle single end reads, there is no mate pair. Hence, set as 0
	TLEN = 0;												// Not yet developed to assign TLEN values. Currently, all set to 0
	QUAL = "*"; 											// Not yet developed to assign QUAL values. Currently, all set to *
	
	result = QNAME + "," + std::to_string(FLAG) + "," + std::to_string(POS) + "," + std::to_string(MAPQ) + "," + CIGAR + "," + RNEXT + "," + std::to_string(PNEXT) + "," + std::to_string(TLEN) + "," + SEQ + "," + QUAL + "," + RSEQ;
	
	return result;
}
	