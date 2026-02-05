#include "count_guides_in_file.h"

// This function takes a metafile (a csv file) with path as input
// This function returns the number of lines (=number of guides) present in metafile
// In the metafile, each guide & its sequence MUST be present in a single row
// The metafile must (i) have no headers (ii) have gene names in 1st column (iii) have guide sequence in 2nd column
/* This is an example metafile file with no headers, with gene names in 1st column and guide sequence in 2nd column
Ptprh	GTATGGTTGTTCTGTGGTC
Ptprh	GTGGGAACTGTGGGGTTTGG
Ptprh	GTACCAGTGTCACAGTGGAC
Ptprh	GAAAACTCGTATGAAGAAGC
Ptprh	GCGGTGGTTACAGGAATCTG
Ptprh	GCACTTCTGGGGGATATG
Ptprh	GACCTGGTCCTGCCCATC
Ptprh	GCACTGGGCCCAGTAGATG
Ptprh	GAACATTTATATAAAAACCC
Ptprh	GTCTCAAATAACGGCAT
Cnr1	GTGGGAAGTATCCTAATT
Cnr1	GTCTGTGGTGATGGTACGGA
Cnr1	GTTGGTTGTGTCTCCTGC*/

int count_guides_in_file (std::ifstream& metafile){
	
	int number_of_lines = 0;
    std::string line;
    if (metafile.is_open()){
	    while (std::getline(metafile, line)){
			number_of_lines = number_of_lines + 1;
		}
	}
	metafile.close();
	std::cout << "Total Number of Guides: " << number_of_lines << "\n";
	
	return number_of_lines;
}