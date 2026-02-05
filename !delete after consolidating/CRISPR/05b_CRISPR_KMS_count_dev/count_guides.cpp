#include "count_guides.h"

// This function returns the number of lines (=number of guides) present in metafile

int count_guides (std::ifstream& metafile){
	
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "  COUNTING GUIDES IN LIBRARY \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	
	int number_of_lines = 0;
    std::string line;
    if (metafile.is_open()){
	    while (std::getline(metafile, line)){
			number_of_lines = number_of_lines + 1;
		}
	}
	metafile.close();
	std::cout << "Total Number of Guides: " << number_of_lines << "\n" << std::flush;
	
	return number_of_lines;
}