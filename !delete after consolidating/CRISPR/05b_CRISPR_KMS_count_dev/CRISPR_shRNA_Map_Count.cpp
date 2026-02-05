#include "count_guides.h"
#include "count_reads.h"
#include "smallest_guide.h"
#include "reverse_complement.h"
#include "find_orientation.h"
#include "position_of_seed_in_reference.h"
#include <iostream>			
#include <fstream>          // you need this for opening, reading, writing file
#include <string>           // you need this for string operations
#include <algorithm>        // you need this for find()
#include <chrono>           // you need this to calculate time taken by the script
#include <cmath> 
#include <vector>
#include <boost/algorithm/string.hpp>

// Use https://www.jdoodle.com/online-compiler-c++ for testing
// If you do not use std::flush, std::cout will be printed ONLY after program ends.
// Since we want to see the progress of the script, use std::flush 

// The 1st argument is a fasta file with path. Fastq files wont work with this script
// The 2nd argument is a csv file with path in which results will be saved by this script
/* The 3rd argument is a csv file with path that has (i) no headers (ii) sgRNA ids in 1st column (iii) guide sequence in 2nd column
Ptprh	GTATGGTTGTTCTGTGGTC
Ptprh	GTGGGAACTGTGGGGTTTGG
Cnr1	GTTGGTTGTGTCTCCTGC*/
// The 4th argument is the number of bases to clip at 5' end before mapping
// The 5th argument is the number of bases to clip at 3' end before mapping
// The 6th argument is the number of maximum mismatches allowed
// The 7th argument determines whether to run script on all reads or first 1000 reads

int main (int argc, char** argv) {
	
	// argv[1] is input fasta file with path
	// argv[2] is filename with path but without file extension like .csv, .sam
	// argv[3] is csv metafile with guide or shRNA sequence and corresponding sgRNA or shRNA id
	std::string fasta = argv[1];
	unsigned int head_trim = atoi(argv[4]);   		// number of bases to ignore at 5' end of read before starting search for guide or shRNA
	unsigned int tail_trim = atoi(argv[5]);   		// number of bases to ignore at 3' end of read before starting search for guide or shRNA
	unsigned int max_errors = atoi(argv[6]);   		// maximum number of mismatches allowed between guide/shRNA sequence and a read; DEFAULT = 3 
	std::string troubleshoot_mode = argv[7];     	// "T" (test code on first 1000 reads in fasta file), "F" (run code on entire fasta file)
  	//std::cout << fasta << "\n";
	//std::cout << argv[1] << "\n";
	
	// Identify number of guides present in library	
  	std::ifstream metafile;                             // create a variable "metafile" of type ifstream
	metafile.open(argv[3], std::ios::in);               // connect the ifstream variable "metafile" to a csv file with guide details using the open()
	int line_count = count_guides(metafile);			// CALLING count_guides_in_file(). We need to know the number of guides to create the arrays below	
	
	// Identify number of reads present in fasta file	
  	std::ifstream fastafile;                            // create a variable "fastafile" of type ifstream
	fastafile.open(argv[1], std::ios::in);              // connect the ifstream variable "metafile" to a csv file with guide details using the open()
	int read_count = count_reads(fastafile);			// CALLING count_reads(). We need to know the number of reads to create the arrays below	
	
	//read_count = 1000;
	
	// Declare variables
	std::string sgRNA[line_count];						// create an array to store sgRNA ids
	std::string guide[line_count];						// create an array to store guide sequences
	std::string guide_rc[line_count];					// create an array to store reverse complementary guide sequences
	int start[line_count];								// create an array to store start position of matching guide sequence in pseudoreference
	int end[line_count];								// create an array to store end position of matching guide sequence in pseudoreference
	int count[line_count];								// create an array to store the counts of each guide sequence
	std::string reference_s = "X";						// create a string to store pseudoreference; initialize it as X
	std::string reference_rc = "X";						// create a string to store reverse complementary pseudoreference; initialize it as X
	
	// Generate pseudoreference and prepare csv file to store counts
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "  PREPARING PSEUDOREFERENCE  \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
    
			
	metafile.open(argv[3], std::ios::in);               	// connect the ifstream variable "metafile" to a file using the open()
	for (int i=0; i < line_count; ++i){ 
		std::getline(metafile, sgRNA[i], ',');          	// read a line from ifstream variable "metafile" until delimiter "," and store in array "sgRNA" at ith position
		std::getline(metafile, guide[i]);               	// read a line from ifstream variable "metafile" until default delimiter "\n" and store in array "guide" at ith position
		int end_of_line = guide[i].find('\r');          	// find length of guide using index position of carriage return "\r"
		guide[i] = guide[i].substr(0,end_of_line);      	// HOW DOES THIS WORK??  why not just use guide[i]??
		guide_rc[i] = reverse_complement(guide[i]);			// CALLING reverse_complement()
		reference_s = reference_s + guide[i] + 'X';			// Add every guide sequence to reference and end each it with an X
		reference_rc = reference_rc + guide_rc[i] + 'X';	// Add every reverse complementary guide sequence to reverse complementary reference and end it with an X
		if (i==0){
			start[i] = reference_s.find(guide[i]);		
		}
		else{
			start[i] = reference_s.find(guide[i], end[i-1]);// Store starting position of each guide in pseudoreference
		}	
		           
		end[i] = start[i] + guide[i].length() - 1;      	// Store ending postion of each guide in pseudoreference
		count[i] = 0;                                   	// Set counts for all guides as 0 initially
	}
	metafile.close();                                   	// disconnect the ifstream variable "metafile" using the close()
	std::cout << "Pseudoreference generated\n" << std::flush;
	
	// Identify length of smallest guide
	//int number_of_guides = sizeof(guide)/sizeof(guide[0]);     // HOW DOES THIS WORK?? Each element in array "guide" has same size irrespective of guide length??
	int min_guide_length = smallest_guide(guide, line_count);	// CALLING count_guides_in_file(). We need to know the number of guides to create the arrays below
	
	// Calculate seed length	
	//int seed_len=std::ceil(min_guide_length - 2)/3;
	int seed_len = min_guide_length/3;	
	seed_len = std::max(seed_len, 8);
	seed_len = std::min(seed_len, 20);	
	std::cout << "Seed Length: " << seed_len << "\n" << std::flush;
	
	// Calculate minimum score allowed
	int min_match=min_guide_length-max_errors;
	std::cout << "Minimum Score Required: " << min_match << "\n" << std::flush;
	
	// Map 1K reads to determine orientation of reference	
	std::string pseudogenome;
	pseudogenome = find_orientation(argv[1], head_trim, tail_trim, reference_s, reference_rc, seed_len, max_errors, min_match, troubleshoot_mode);			// CALLING find_orientation()
 	//pseudogenome = reference_s;
	
	// Start timer before starting mapping
	auto begin = std::chrono::high_resolution_clock::now();
	
	//Map the reads using proper orientation of reference identified above
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "        MAPPING READS        \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;	
	
	int read_number = 0;                    // create an integer varible to store the number of the current read	
	std::string read_id;				// create a string to store the name of the current read
	std::string read_seq;               // create a string to store the sequence of the current read
	std::string read_third_line;        // create a string to store the 3rd line (usually a +) in the current read
	std::string read_qual;              // create a string to store the quality scores of the current read	
	std::string trimmed_read;
	int multimapped = 0; 					// create an integer variable to store the number of multimapped reads
	int unmapped = 0;               		// create an integer variable to store the number of unmapped reads	
	
	std::string alignment;
	std::string QNAME[read_count];						// create an array to store the name of the current read
	int FLAG[read_count];                				// create an array to store the FLAG values of the current read
	std::string RNAME[read_count];						// create an array to store the gene/guide mapped to the current read
	int POS[read_count];                				// create an array to store the mapping start position of the current read in reference
	int MAPQ[read_count];                				// create an array to store the MAPQ values of the current read
	std::string CIGAR[read_count];						// create an array to store the CIGAR alignment values of the current read
	std::string RNEXT[read_count];                		// create an array to store the RNEXT values of the current read
	int PNEXT[read_count];                				// create an array to store the PNEXT values of the current read
	int TLEN[read_count];                				// create an array to store the TLEN values of the current read
	std::string SEQ[read_count];        				// create an array to store the sequence of the current read	
	std::string QUAL[read_count];       				// create an array to store the quality scores of the current read	
	std::string RSEQ[read_count];        				// create an array to store the sequence of the current read
	
	std::ifstream file_in;                  // create a variable "file_in" of type ifstream
    file_in.open(argv[1],std::ios::in);  	// connect the ifstream variable "file_in" to fasta file using the open()	
	if (file_in.is_open()){
		
		std::getline(file_in,read_id);                 			// read 1st line from file_in, store it as QNAME		
		while (read_id.find("@") == 0){               			// WHILE loop will end once "end of file" is reached			
			std::getline(file_in,read_seq);              			// read 2nd line from file_in, store it as SEQ
			std::getline(file_in,read_third_line);         					// read 3rd line from file_in (it will be +)
			std::getline(file_in,read_qual);						// read 4th line from file_in, store it as QUAL
			trimmed_read = read_seq.substr(head_trim,(read_seq.length()-tail_trim)); 	// trim the read based on user input
			
			/*if (troubleshoot_mode=="T"){
				std::cout << read_id << "\n" << std::flush;		// USE THIS LINE FOR TROUBLESHOOTING
				std::cout << SEQ[read_number] << "\n" << std::flush;		// USE THIS LINE FOR TROUBLESHOOTING
				std::cout << read_third_line << "\n" << std::flush;			// USE THIS LINE FOR TROUBLESHOOTING
				std::cout << read_qual << "\n" << std::flush;		// USE THIS LINE FOR TROUBLESHOOTING
			}*/
			//if (read_number > 24000000){										// USE THIS LINE FOR TROUBLESHOOTING
			alignment = position_of_seed_in_reference(trimmed_read, pseudogenome, read_number, seed_len, max_errors, min_match, troubleshoot_mode);			// CALLING position_of_read_in_reference()	
			/*} 
			else{
				alignment = "Read_" + std::to_string(read_number) + "," + "0" + "," + "0" + "," + "0" + "," + "*" + "," + "*" + "," + "0" + "," + "0" + "," + "*" + "," + "*" + "," + "*";
			}*/
			
			std::vector <std::string> result;
			boost::split(result, alignment, boost::is_any_of(","));
			QNAME[read_number] = result[0];
			FLAG[read_number] = std::stoi(result[1]);
			POS[read_number] = std::stoi(result[2]);
			MAPQ[read_number] = std::stoi(result[3]);
			CIGAR[read_number] = result[4];
			RNEXT[read_number] = result[5];
			PNEXT[read_number] = std::stoi(result[6]);
			TLEN[read_number] = std::stoi(result[7]);
			SEQ[read_number] = result[8];
			QUAL[read_number] = result[9];
			RSEQ[read_number] = result[10];						
			
			if (POS[read_number] > 0){										// If a proper match of read in reference is found,
				for (int i=0; i<line_count; ++i){							// find which guide it maps to, based on its position in reference
					if ((start[i] <= POS[read_number]) && (end[i] > POS[read_number])){
						count[i] = count[i]+1;								// once guide is identified, increase its count by 1
						RNAME[read_number] = sgRNA[i];							
						break;												// exit from FOR loop
					}
				}
			}
			else if(POS[read_number] == -1){
				unmapped = unmapped+1;										// If a proper match of read in reference is NOT found, increase unmapped read count by 1
				RNAME[read_number] = "Unmapped";				
			}
			else if(POS[read_number] ==  -10){
				multimapped = multimapped+1;
				RNAME[read_number] = "Multimapped";				
			}
			
			read_number = read_number + 1;             						// increase read_number by 1
			std::getline(file_in,read_id);                		// read 5th line from file, store it as QNAME

			if (read_number%1000000==0){
				int quotient = read_number/1000000; 
				std::cout << "Mapped " << quotient << "M reads\n" << std::flush;			
			}	
		}	
	}
	file_in.close();
	
	// Stop timer after finishing mapping
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - begin);
	
	// Save SAM like file
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "     SAVING SAM LIKE FILE    \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	
	std::string samfilename = argv[2] + std::string(".sam");
	std::ofstream samfile;									// create a sam like file to store alignment info
	samfile.open(samfilename, std::ios::out);
	samfile << "QNAME" << "," << "FLAG" << "," << "RNAME" << "," << "POS" << "," << "MAPQ" << "," << "CIGAR" << "," << "RNEXT" << "," << "PNEXT" << "," << "TLEN" << "," << "SEQ" << "," << "QUAL" << "," << "RSEQ" << "\n"; 
	for (int i=0; i < read_count; ++i){
		samfile << QNAME[i] << "," << FLAG[i] << "," << RNAME[i] << "," << POS[i] << "," << MAPQ[i] << "," << CIGAR[i] << "," << RNEXT[i] << "," << PNEXT[i] << "," << TLEN[i] << "," << SEQ[i] << "," <<  QUAL[i] << "," << RSEQ[i] << "\n"; 
	}
	samfile.close();	
	
	// Save counts file
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "    SAVING COUNTS TO FILE    \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	
	std::string countfilename = argv[2] + std::string(".csv");
	std::ofstream countfile;									// create a countfile to store the count data
	countfile.open(countfilename, std::ios::out);
	countfile << "sgRNA" << "," << "GUIDE" << "," << "REVERSE COMPLEMENT GUIDE" << "," << "START" << "," << "END" << "," <<  "COUNT" << "\n"; 
	for (int i=0; i < line_count; ++i){
		countfile << sgRNA[i] << "," << guide[i] << "," << guide_rc[i] << "," << start[i] << "," <<  end[i] << "," <<  count[i] << "\n"; 
	}
	countfile.close();	
	
	std::cout << read_number << " Reads Analyzed in " << duration.count() << " seconds\n" << std::flush;
	std::cout << "Number of Multimapped Reads: " << multimapped << "\n" << std::flush;
	std::cout << "Number of Unmapped Reads: " << unmapped << "\n" << std::flush;			
	
	return 0;
}




