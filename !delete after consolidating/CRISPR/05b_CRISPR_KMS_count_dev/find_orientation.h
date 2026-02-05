#ifndef FIND_ORIENTATION_H
#define FIND_ORIENTATION_H

#include <string>  
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
std::string find_orientation(std::string fasta, int head_trim, int tail_trim, std::string reference_s, std::string reference_rc, int seed_len, int max_errors, int min_match, std::string troubleshoot_mode);
#endif