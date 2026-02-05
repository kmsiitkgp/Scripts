#include "position_of_seed_in_reference.h"
#include "find_orientation.h"

// This function takes a fasta file, head_trim, tail_trim, 
// reference (i.e. concatenated guide sequences), reverse complemented reference (i.e. concatenated guide_rc sequences), 
// seed length, max errors allowed, minimum match score and troubleshooting mode as input.
// This function calls position_of_seed_in_reference() internally.
// This function returns the reference genome in an orientation where reads map without reverse complementation.

std::string find_orientation(std::string fasta, int head_trim, int tail_trim, std::string reference_s, std::string reference_rc, int seed_len, int max_errors, int min_match, std::string troubleshoot_mode){
	
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "IDENTIFYING GUIDE ORIENTATION\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;

	
	std::string read_id;				// create a string to store the name of the current read
	std::string read_seq;               // create a string to store the sequence of the current read
	std::string read_third_line;        // create a string to store the 3rd line (usually a +) in the current read
	std::string read_qual;              // create a string to store the quality scores of the current read
	int read_number = 0;                // create an integer varible to store the number of the current read 
	std::string read;                   // create a string to store the trimmed sequence of the current read
	std::string alignment_f;
	std::string alignment_r;
	int unmapped_f = 0;             	// create an integer variable to store the number of unmapped reads
	int unmapped_r = 0;             	// create an integer variable to store the number of unmapped reads
	
	std::ifstream file_in;              // create a variable "file_in" of type ifstream
	file_in.open(fasta,std::ios::in);   // connect the ifstream variable "file_in" to fasta file using the open()
	
	if (file_in.is_open()){		
		std::getline(file_in,read_id);                 						// read 1st line from file_in, store it as read_id	
		while (read_id.find("@") == 0){               						// WHILE loop will end once "end of file" is reached
			read_number = read_number + 1;             						// increase read_number by 1
			std::getline(file_in,read_seq);              					// read 2nd line from file_in, store it as read_seq
			std::getline(file_in,read_third_line);         					// read 3rd line from file_in (it will be +)
			std::getline(file_in,read_qual);								// read 4th line from file_in, store it as read_qual
			read = read_seq.substr(head_trim,read_seq.length()-tail_trim); 	// trim the read based on user input		
			
			// Map using forward orientated reference			
			alignment_f = position_of_seed_in_reference(read, reference_s, read_number, seed_len, max_errors, min_match, troubleshoot_mode);		// CALLING position_of_seed_in_reference()
			std::vector <std::string> result_f;
			boost::split(result_f, alignment_f, boost::is_any_of(","));	
			if (std::stoi(result_f[2]) == -1){									// If a proper match of read in reference is NOT found,		
				unmapped_f = unmapped_f+1;						// increase unmapped read count by 1
			}
			
			// Map using reverse complemented reference
			alignment_r = position_of_seed_in_reference(read, reference_rc, read_number, seed_len, max_errors, min_match, troubleshoot_mode);	// CALLING position_of_seed_in_reference()
			std::vector <std::string> result_r;
			boost::split(result_r, alignment_r, boost::is_any_of(","));			
			if (std::stoi(result_r[2]) == -1){									// If a proper match of read in reference is NOT found,
				unmapped_r = unmapped_r+1;						// increase unmapped read count by 1
			}			
			
			std::getline(file_in,read_id);                		// read 5th line from file, store it as read_id
			if (read_number==100){								// Exit loop once 1000 reads have been analyzed
				break;											
			}
		}
	}
	file_in.close();

	std::cout << "Number of Unmapped Reads using forward orientated reference: " << unmapped_f << "\n" << std::flush;
	std::cout << "Number of Unmapped Reads using reverse complemented reference: " << unmapped_r << "\n" << std::flush;

	// Determine orientation based on unmapped read counts
	std::string reference;
	if (unmapped_f < unmapped_r){		
		reference = reference_s;
		std::cout << "Using forward orientated reference\n" << std::flush;		
	}
	else {
		reference = reference_rc;
		std::cout << "Using reverse complemented reference\n" << std::flush;
	}

	return reference;
}
