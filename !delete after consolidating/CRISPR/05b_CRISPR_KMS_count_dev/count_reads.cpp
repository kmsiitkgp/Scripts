#include "count_reads.h"

// This function returns the number of lines (=number of reads) present in fastafile

int count_reads (std::ifstream& fastafile){
	
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << " COUNTING READS IN FASTA FILE\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	
	int number_of_lines = 0;
    std::string line;
    if (fastafile.is_open()){
	    while (std::getline(fastafile, line)){
			number_of_lines = number_of_lines + 1;
		}
	}
	fastafile.close();	
	
	number_of_lines=number_of_lines/4;	
	std::cout << "Total Number of Reads: " << number_of_lines << "\n" << std::flush;
	
	return number_of_lines;
}