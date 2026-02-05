#include "reverse_complement.h"

// This function takes a guide sequence as input
// This function returns the reverse complement of guide sequence as output
std::string reverse_complement (std::string guide){
	
	std::string reverse_complement_guide = "";		// initialize reverse complement of guide as empty string
    for (int j=0; j < guide.length(); ++j){			// loop through each letter of seq
		if (guide[j] == 'A'){
			reverse_complement_guide = 'T' + reverse_complement_guide;
		}
		else if (guide[j] == 'T'){
			reverse_complement_guide = 'A' + reverse_complement_guide;
		}
		else if (guide[j] == 'G'){
			reverse_complement_guide = 'C' + reverse_complement_guide;
		}
		else if (guide[j] == 'C'){
			reverse_complement_guide = 'G' + reverse_complement_guide;
		}
		else if (guide[j] == 'N'){
			reverse_complement_guide = 'N' + reverse_complement_guide;
		}
		else {
			reverse_complement_guide = reverse_complement_guide;
		}
    }
	//std::cout << reverse_complement_guide << "\n";
	return reverse_complement_guide;
}