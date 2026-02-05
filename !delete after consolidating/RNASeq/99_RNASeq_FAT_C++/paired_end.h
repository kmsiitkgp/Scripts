#ifndef PAIRED_END_H
#define PAIRED_END_H

#include <iostream>
#include <fstream>          // you need this for opening, reading, writing file
#include <string>           // you need this for string operations
#include <algorithm>        // you need this for find()
#include <chrono>           // you need this to calculate time taken by the script

void pe_adapter_trimmer(unsigned int min_len, unsigned int min_q_score, unsigned int base, std::ifstream& file_1, std::string f_adap_seq, std::ofstream& p_f_file, std::ofstream& up_f_file, std::ifstream& file_2, std::string r_adap_seq, std::ofstream& p_r_file, std::ofstream& up_r_file);

#endif


