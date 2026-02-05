#include "count_guides_in_file.h"
#include "reverse_complement.h"
#include "position_of_read_in_reference.h"
#include "position_of_seed_in_reference.h"
#include "seed_length_calculator.h"
#include <iostream>			
#include <fstream>          // you need this for opening, reading, writing file
#include <string>           // you need this for string operations
#include <algorithm>        // you need this for find()
#include <chrono>           // you need this to calculate time taken by the script

int main (int argc, char** argv) {

	auto begin = std::chrono::high_resolution_clock::now();
	// argv[1] is input fastq.gz file with path
	// argv[2] is output csv file with path
	// argv[3] is csv metafile with guide or shRNA sequence and corresponding gene
	unsigned int head_trim = atoi(argv[4]);   		// number of bases to ignore at 5' end of read before starting search for guide or shRNA
	unsigned int tail_trim = atoi(argv[5]);   		// number of bases to ignore at 3' end of read before starting search for guide or shRNA
	unsigned int threshold = atoi(argv[6]);   		// maximum number of mismatches allowed between guide or shRNA sequence and a read; DEFAULT = 3 
	std::string orientation = argv[7];				// can ONLY take values F or R.             	
	// F means guide sequence present in reads. R indicates reverse complement of guide sequence present in reads
  	
  	std::ifstream metafile;
	metafile.open(argv[3], std::ios::in);
	int line_count = count_guides_in_file(metafile);	// CALLING count_guides_in_file()
	
	std::string gene[line_count];						// create a list to store gene names
	std::string guide[line_count];						// create a list to store guide sequences
	std::string guide_rc[line_count];					// create a list to store reverse complementary guide sequences
	int start[line_count];								// create a list to store start position of matching guide sequence
	int end[line_count];								// create a list to store end position of matching guide sequence
	int count[line_count];								// create a list to store the counts of each guide sequence
	std::string reference_s = "X";						// initialize reference as X
	std::string reference_rc = "X";						// initialize reverse complementary reference as X
	
	metafile.open(argv[3], std::ios::in);
	for (int i=0; i < line_count; ++i){ 
		std::getline(metafile, gene[i], ',');
		std::getline(metafile, guide[i]);
		int end_of_line = guide[i].find('\r');
		guide[i] = guide[i].substr(0,end_of_line);
		guide_rc[i] = reverse_complement(guide[i]);		// CALLING reverse_complement()
		reference_s = reference_s + guide[i] + 'X';		// Add every guide sequence to reference and end each it with an X
		reference_rc = reference_rc + guide_rc[i] + 'X';// Add every reverse complementary guide sequence to reverse complementary reference and end it with an X
		start[i] = reference_s.find(guide[i]);
		end[i] = start[i] + guide[i].length() - 1;
		count[i] = 0;
	}
	metafile.close();
	
	int number_of_guides = sizeof(guide)/sizeof(guide[0]);
	int seed_length = seed_length_calculator(guide, number_of_guides);	// CALLING seed_length_calculator ()
		
	std::string pseudogenome;
	if (orientation == "F" || orientation == "f" || orientation == "FORWARD" || orientation == "forward" || orientation == "Forward"){
		pseudogenome = reference_s;
	}
	else if (orientation == "R" || orientation == "r" || orientation == "REVERSE" || orientation == "reverse" || orientation == "Reverse"){
		pseudogenome = reference_rc;
	}
	else{
		std::cout << "Please enter F or R for orientation\n";
	}
	
	std::ifstream file_in;
    file_in.open(argv[1],std::ios::in);
	
	int unmapped_read = 0;
	std::string read_id;
	std::string read_seq;
	std::string read_third_line;
	std::string read_qual;
	int read_number = 0;
	std::string read;
	
	//auto begin = std::chrono::high_resolution_clock::now();
	
	if (file_in.is_open()){
		std::getline(file_in,read_id);                 						// read 1st line from file_in, store it as read_id
		
		while (read_id.find("@") == 0){               						// WHILE loop will end once "end of file" is reached
			read_number = read_number + 1;             						// increase read_number by 1
			std::getline(file_in,read_seq);              					// read 2nd line from file_in, store it as read_seq
			std::getline(file_in,read_third_line);         					// read 3rd line from file_in (it will be +)
			std::getline(file_in,read_qual);								// read 4th line from file_in, store it as read_qual
			read = read_seq.substr(head_trim,read_seq.length()-tail_trim); 	// trim the read based on user input
			//std::cout << read_id << "\n";									// USE THIS LINE FOR TROUBLESHOOTING
			//std::cout << read_seq << "\n";								// USE THIS LINE FOR TROUBLESHOOTING
			//std::cout << read_third_line << "\n";							// USE THIS LINE FOR TROUBLESHOOTING
			//std::cout << read_qual << "\n";								// USE THIS LINE FOR TROUBLESHOOTING
			
			int match = position_of_read_in_reference(read, pseudogenome, read_number, seed_length, threshold);	// CALLING position_of_read_in_reference()
			if (match != -1){												// If a proper match of read in reference is found,
				for (int i=0; i<line_count; ++i){							// find which guide it maps to, based on its position in reference
					if ((start[i] <= match) && (end[i] > match)){
						count[i] = count[i]+1;								// once guide is identified, increase its count by 1
						break;												// exit from FOR loop
					}
				}
			}
			else{
				unmapped_read = unmapped_read+1;							// If a proper match of read in reference is NOT found, increase unmapped read count by 1
			}
			std::getline(file_in,read_id);                					// read 5th line from file, store it as read_id
			//if (read_number==1000){											// USE THIS LINE FOR TROUBLESHOOTING
			//	break;														// USE THIS LINE FOR TROUBLESHOOTING
			//}																// USE THIS LINE FOR TROUBLESHOOTING
		}
	}
	file_in.close();
	
	std::ofstream countfile;												// create a countfile to store the count data
	countfile.open(argv[2], std::ios::out);
	countfile << "GENE" << "," << "GUIDE" << "," << "REVERSE COMPLEMENT GUIDE" << "," << "START" << "," << "END" << "," <<  "COUNT" << "\n"; 
	for (int i=0; i < line_count; ++i){
		countfile << gene[i] << "," << guide[i] << "," << guide_rc[i] << "," << start[i] << "," <<  end[i] << "," <<  count[i] << "\n"; 
	}
	countfile.close();
	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - begin);
	
	std::cout << read_number << " Reads Analyzed in " << duration.count() << " seconds\n";
	std::cout << "Number of Unmapped Reads: " << unmapped_read << "\n";
	std::cout << "Seed Length: " << seed_length << "\n";
	
	return 0;
}




