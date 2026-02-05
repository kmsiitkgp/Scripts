#include "smallest_guide.h"

// This function takes a list of guides as input and identifies length of smallest guide

int smallest_guide(std::string guide[], int number_of_guides){ 
/* You cannot pass an array to a function in C++. Only the address of 1st element &guide[0] 
is being passed to the function. So, you need to indicate the size of the array as separate parameter
as you cannot calculate size of array within the function using int number_of_guides = sizeof(guide)/sizeof(guide[0]); */
	
	std::cout << "\n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	std::cout << "  IDENTIFYING SMALLEST GUIDE \n" << std::flush;
	std::cout << "*****************************\n" << std::flush;
	
	int max = 0;
	int min = 1000;	
	
    for(int i=0;i<number_of_guides;++i){
		if (guide[i].length() > max){
			max = guide[i].length();
		}
	
		if (guide[i].length() < min){
			min=guide[i].length();
		}
	}	
	
	std::cout << "Smallest guide: " << min << "\n" << std::flush;
	std::cout << "Longest guide: " << max << "\n" << std::flush;

	
	return min;
}