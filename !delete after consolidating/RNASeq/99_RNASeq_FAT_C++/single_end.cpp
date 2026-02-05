#include "single_end.h"
#include "quality_trimmer.h"
#include "adapter_trimmer.h"

// This function takes minimum length of trimmed read, minimum quality score, base, input fastq filestream, adapter sequence, output trimmed read filestream and output untrimmed read filestream as input
// It creates 2 files: one with trimmed reads and one with untrimmed reads (i.e. reads that if trimmed were smaller than minimum length specified
void se_adapter_trimmer(unsigned int min_len, unsigned int min_q_score, unsigned int base, std::ifstream& file_1, std::string f_adap_seq, std::ofstream& p_f_file, std::ofstream& up_f_file){

  auto start = std::chrono::high_resolution_clock::now();

  unsigned int f_read_count = 0;  	// number of reads before trimming
  unsigned int f_adap_trim = 0;   	// number of adapter trimmed reads
  unsigned int f_polyG_trim = 0;  	// number of polyG trimmed reads
  unsigned int f_qual_trim = 0;   	// number of quality trimmed reads
  unsigned int read_lost = 0;   	// number of paired reads lost

  std::string f_id;
  std::string f_seq;
  std::string f_third_line;
  std::string f_qual;
  std::string f_seq_orig;
  std::string f_qual_orig;
  std::string f_seed_seq;
  
  unsigned int f_polyG_loc;
  unsigned int f_adap_loc;
  
  // Find the shortest seed sequence of adapter that is unique within the adapter
  for (int i = 2; i < f_adap_seq.length(); ++i){	// if i=1, adap_seq.substr(0,1) will be 1 base, which will always be present in adap_seq	
	if (f_adap_seq.substr(i,f_adap_seq.length()).find(f_adap_seq.substr(0,i)) > f_adap_seq.length()){
	  f_seed_seq = f_adap_seq.substr(0,i);
	  break;
	}
  }

  if (file_1.is_open()){
    std::getline(file_1,f_id);                 							// read 1st line from file_1, store it as f_id
    while (f_id.find("@") == 0){               							// loop will end once "end of file" is reached
      //std::cout << f_id << "\n";											// USE THIS LINE FOR TROUBLESHOOTING
	  f_read_count = f_read_count + 1;             						// increase f_read_count by 1
      std::getline(file_1,f_seq);              							// read 2nd line from file_1, store it as f_seq
      std::getline(file_1,f_third_line);         						// read 3rd line from file_1 (it will be +)
      std::getline(file_1,f_qual);             							// read 4th line from file_1, store it as f_qual
      f_seq_orig = f_seq;                                       		// store original f_seq in f_seq_orig
      f_qual_orig = f_qual;                                     		// store original f_qual in f_qual_orig

      f_polyG_loc = f_seq.find("GGGGGGGGGGGGGGGGGGGG");         		// find 1st location where polyG is present
      if (f_polyG_loc < f_seq.length()){                        		// if polyG is "present"
        f_polyG_trim = f_polyG_trim + 1;                            		// increase f_polyG_trim count by 1
        f_seq = f_seq.substr(0,f_polyG_loc);                    		// store the trimmed seq to f_seq
        f_qual = f_qual.substr(0,f_polyG_loc);                  		// store the trimmed quality to f_qual
        //std::cout << "\nOriginal Seq:\t\t" << f_seq_orig;          		// USE THIS LINE FOR TROUBLESHOOTING
        //std::cout << "\npolyG trimmed Seq:\t" << f_seq;            		// USE THIS LINE FOR TROUBLESHOOTING
        //break;                                                		// USE THIS LINE FOR TROUBLESHOOTING
      }

      if (f_seq.length() > 0){											//To avoid segmentation fault when f_seq.length==0 
		f_seq = adapter_trimmer(f_seq, f_adap_seq, f_seed_seq);
        if (f_qual.length() > f_seq.length()){                       	// if adapter trimmed, length of f_seq < f_qual
          f_adap_trim = f_adap_trim + 1;                                    // increase f_adap_trim count by 1
          f_qual = f_qual.substr(0,f_seq.length());                     // store the trimmed qual to f_qual
	      //std::cout << "\nOriginal Seq:\t\t" << f_seq_orig;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nOriginal Qual:\t\t" << f_qual_orig;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nAdapter trimmed Seq:\t" << f_seq;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nAdapter trimmed Qual:\t" << f_qual;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\n";                                               // USE THIS LINE FOR TROUBLESHOOTING
          //break;                                                      // USE THIS LINE FOR TROUBLESHOOTING
	    }
	  }
	  
	  if (f_seq.length() > 0){											//To avoid segmentation fault when f_seq.length==0
	    f_qual = quality_trimmer(f_qual, min_q_score, base);            
        if (f_seq.length() > f_qual.length()){                          // if quality trimmed, length of f_qual < f_seq
          f_qual_trim = f_qual_trim + 1;                                    // increase f_qual_trim count by 1
          f_seq = f_seq.substr(0,f_qual.length());                      // store the trimmed seq to f_seq
          //std::cout << "\nOriginal Seq:\t\t" << f_seq_orig;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nOriginal Qual:\t\t" << f_qual_orig;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Seq:\t" << f_seq;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Qual:\t" << f_qual;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\n";                                               // USE THIS LINE FOR TROUBLESHOOTING
          //break;                                                      // USE THIS LINE FOR TROUBLESHOOTING
        }
      }
	  
      if (f_seq.length() < min_len){                                    // if length of f_seq < min_len
        read_lost = read_lost + 1;                                      // increase read_lost count by 1
        up_f_file << f_id;                                         		// write read_id
        up_f_file << "\n";                                         		// write new line
        up_f_file << f_seq_orig;                                   		// write f_seq_orig
        up_f_file << "\n";                                         		// write new_line
        up_f_file << "+";                                          		// write +
        up_f_file << "\n";                                         		// write new_line
        up_f_file << f_qual_orig;                                  		// write f_qual_orig
        up_f_file << "\n";                                         		// write new_line
        //break;                                                       	// USE THIS LINE FOR TROUBLESHOOTING
      }
      else{
        p_f_file << f_id;
        p_f_file << "\n";
        p_f_file <<f_seq;
        p_f_file << "\n";
        p_f_file << "+";
        p_f_file << "\n";
        p_f_file << f_qual;
        p_f_file << "\n";
      }

      std::getline(file_1,f_id);                						// read 5th line from file, store it as read_id
    }
  }
  else std::cout << "Unable to open file\n";

  file_1.close();														// close file_1
  up_f_file.close();													// close up_f_file
  p_f_file.close();														// close p_f_file

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

  //In %10s,  10 indicates 10 blank spaces and s indicates the variable isd a string. The string will be stored in these 10 spaces and right-justified
  //In %10.2f , 10 indicates 10 blank spaces and .2f indicates variable must be truncated to 2 decimal places.
  // The decimal number will be stored in these 10 spaces and right-justified.
  printf("Trimming Statistics      | %10s\t%10s\n", "# Forward", "% Forward");	
  printf("Reads Before Trimming    | %10i\t%10.2f\n", f_read_count, 				100.00);
  printf("Reads After Trimming     | %10i\t%10.2f\n", (f_read_count - read_lost), 	double(f_read_count-read_lost)*100/double(f_read_count));
  printf("Reads Lost               | %10i\t%10.2f\n", read_lost, 					double(read_lost)*100/double(f_read_count));
  printf("Adapter Trimmed Reads    | %10i\t%10.2f\n", f_adap_trim, 					double(f_adap_trim)*100/double(f_read_count));
  printf("polyG Trimmed Reads      | %10i\t%10.2f\n", f_polyG_trim, 				double(f_polyG_trim)*100/double(f_read_count));
  printf("Quality Trimmed Reads    | %10i\t%10.2f\n", f_qual_trim, 					double(f_qual_trim)*100/double(f_read_count));
  std::cout << f_read_count << " Reads Analyzed in " << duration.count() << " seconds\n";
}
