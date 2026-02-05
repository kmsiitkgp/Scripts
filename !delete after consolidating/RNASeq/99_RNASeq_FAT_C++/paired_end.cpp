#include "paired_end.h"
#include "quality_trimmer.h"
#include "adapter_trimmer.h"

// This function takes minimum length of trimmed read, minimum quality score, base, input forward fastq filestream, forward adapter sequence, output forward trimmed read filestream and output forward untrimmed read filestream,
// input reverse fastq filestream, reverse adapter sequence, output reverse trimmed read filestream and output reverse untrimmed read filestream as input
// It creates 4 files: forward trimmed reads,  foward untrimmed & unpaired reads (i.e. reads that if trimmed were smaller than minimum length specified), reverse trimmed reads,  reverse untrimmed & unpaired reads (i.e. reads that if trimmed were smaller than minimum length specified)

void pe_adapter_trimmer(unsigned int min_len, unsigned int min_q_score, unsigned int base, std::ifstream& file_1, std::string f_adap_seq, std::ofstream& p_f_file, std::ofstream& up_f_file, std::ifstream& file_2, std::string r_adap_seq, std::ofstream& p_r_file, std::ofstream& up_r_file){
  
  auto start = std::chrono::high_resolution_clock::now();

  unsigned int f_read_count = 0;  // number of reads before trimming
  unsigned int f_adap_trim = 0;   // number of adapter trimmed reads
  unsigned int f_polyG_trim = 0;  // number of polyG trimmed reads
  unsigned int f_qual_trim = 0;   // number of quality trimmed reads
  unsigned int read_lost = 0;   // number of paired reads lost

  std::string f_id;
  std::string f_seq;
  std::string f_third_line;
  std::string f_qual;
  std::string f_seq_orig;
  std::string f_qual_orig;
  std::string f_seed_seq;
  
  unsigned int f_polyG_loc;
  unsigned int f_adap_loc;
  
  unsigned int r_read_count = 0;  	// number of reads before trimming
  unsigned int r_adap_trim = 0;   	// number of adapter trimmed reads
  unsigned int r_polyG_trim = 0;  	// number of polyG trimmed reads
  unsigned int r_qual_trim = 0;   	// number of quality trimmed reads
    
  std::string r_id;
  std::string r_seq;
  std::string r_third_line;
  std::string r_qual;
  std::string r_seq_orig;
  std::string r_qual_orig;
  std::string r_seed_seq;
  
  unsigned int r_polyG_loc;
  unsigned int r_adap_loc;
  
  // Find the shortest seed sequence of adapter that is unique within the adapter
  for (int i = 2; i < f_adap_seq.length(); ++i){	// if i=1, adap_seq.substr(0,1) will be 1 base, which will always be present in adap_seq	
	if (f_adap_seq.substr(i,f_adap_seq.length()).find(f_adap_seq.substr(0,i)) > f_adap_seq.length()){
	  f_seed_seq = f_adap_seq.substr(0,i);
	  break;
	}
  }
  
  for (int i = 2; i < r_adap_seq.length(); ++i){		
	if (r_adap_seq.substr(i,r_adap_seq.length()).find(r_adap_seq.substr(0,i)) > r_adap_seq.length()){
	  r_seed_seq = r_adap_seq.substr(0,i);
	  break;
	}
  }  
  
  if (file_1.is_open() && file_2.is_open()){
    std::getline(file_1,f_id);                 							// read 1st line from file_1, store it as f_id
    std::getline(file_2,r_id);
  
    while (f_id.find("@") == 0 && r_id.find("@") == 0){     			// loop will end once "end of file" is reached
      //std::cout << f_id << "\n";											// USE THIS LINE FOR TROUBLESHOOTING
	  //std::cout << r_id << "\n";											// USE THIS LINE FOR TROUBLESHOOTING
	  f_read_count = f_read_count + 1;             							
      r_read_count = r_read_count + 1;
	  
	  std::getline(file_1,f_seq);              							
      std::getline(file_1,f_third_line);         							
      std::getline(file_1,f_qual);             							
      f_seq_orig = f_seq;                                       		
      f_qual_orig = f_qual;
                 							
      std::getline(file_2,r_seq);              							
      std::getline(file_2,r_third_line);         							
      std::getline(file_2,r_qual);             							
      r_seq_orig = r_seq;                                       		
      r_qual_orig = r_qual;          	  
	  
	  f_polyG_loc = f_seq.find("GGGGGGGGGGGGGGGGGGGG");         		
      if (f_polyG_loc < f_seq.length()){                        		
        f_polyG_trim = f_polyG_trim + 1;                            		
        f_seq = f_seq.substr(0,f_polyG_loc);                    		
        f_qual = f_qual.substr(0,f_polyG_loc);                  		
        //std::cout << "\nOriginal Seq:\t\t" << f_seq_orig;          		// USE THIS LINE FOR TROUBLESHOOTING
        //std::cout << "\npolyG trimmed Seq:\t" << f_seq;            		// USE THIS LINE FOR TROUBLESHOOTING
        //break;                                                		// USE THIS LINE FOR TROUBLESHOOTING
      }
	  
	  r_polyG_loc = r_seq.find("GGGGGGGGGGGGGGGGGGGG");         		
      if (r_polyG_loc < r_seq.length()){                        		
        r_polyG_trim = r_polyG_trim + 1;                            		
        r_seq = r_seq.substr(0,r_polyG_loc);                    		
        r_qual = r_qual.substr(0,r_polyG_loc);                  		
        //std::cout << "\nOriginal Seq:\t\t" << r_seq_orig;          		// USE THIS LINE FOR TROUBLESHOOTING
        //std::cout << "\npolyG trimmed Seq:\t" << r_seq;            		// USE THIS LINE FOR TROUBLESHOOTING
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
	  
	  if (r_seq.length() > 0){											//To avoid segmentation fault when f_seq.length==0 
		r_seq = adapter_trimmer(r_seq, r_adap_seq, r_seed_seq);
        if (r_qual.length() > r_seq.length()){                       	// if adapter trimmed, length of f_seq < f_qual
          r_adap_trim = r_adap_trim + 1;                                    // increase f_adap_trim count by 1
          r_qual = r_qual.substr(0,r_seq.length());                     // store the trimmed qual to f_qual
	      //std::cout << "\nOriginal Seq:\t\t" << r_seq_orig;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nOriginal Qual:\t\t" << r_qual_orig;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nAdapter trimmed Seq:\t" << r_seq;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nAdapter trimmed Qual:\t" << r_qual;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\n";                                               // USE THIS LINE FOR TROUBLESHOOTING
          //break;                                                      // USE THIS LINE FOR TROUBLESHOOTING
	    }
	  }
	  
	  if (f_seq.length() > 0){
	    f_qual = quality_trimmer(f_qual, min_q_score, base);            
        if (f_seq.length() > f_qual.length()){                          
          f_qual_trim = f_qual_trim + 1;                                    
          f_seq = f_seq.substr(0,f_qual.length());                      
          //std::cout << "\nOriginal Seq:\t\t" << f_seq_orig;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nOriginal Qual:\t\t" << f_qual_orig;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Seq:\t" << f_seq;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Qual:\t" << f_qual;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\n";                                               // USE THIS LINE FOR TROUBLESHOOTING
          //break;                                                      // USE THIS LINE FOR TROUBLESHOOTING
        }
      }

      if (r_seq.length() > 0){
	    r_qual = quality_trimmer(r_qual, min_q_score, base);            
        if (r_seq.length() > r_qual.length()){                          
          r_qual_trim = r_qual_trim + 1;                                    
          r_seq = r_seq.substr(0,r_qual.length());                      
          //std::cout << "\nOriginal Seq:\t\t" << r_seq_orig;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nOriginal Qual:\t\t" << r_qual_orig;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Seq:\t" << r_seq;                // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\nQuality trimmed Qual:\t" << r_qual;              // USE THIS LINE FOR TROUBLESHOOTING
          //std::cout << "\n";                                               // USE THIS LINE FOR TROUBLESHOOTING
          //break;                                                      // USE THIS LINE FOR TROUBLESHOOTING
        }
      }

      if (f_seq.length() < min_len || r_seq.length() < min_len){        // if length of f_seq < min_len
        read_lost = read_lost + 1;                                      // increase read_lost count by 1
        up_f_file << f_id;                                         		
        up_f_file << "\n";                                         		
        up_f_file << f_seq_orig;                                   		
        up_f_file << "\n";                                         		
        up_f_file << "+";                                          		
        up_f_file << "\n";                                         		
        up_f_file << f_qual_orig;                                  		
        up_f_file << "\n";                                         		
        
		up_r_file << r_id;                                         		
        up_r_file << "\n";                                         		
        up_r_file << r_seq_orig;                                   		
        up_r_file << "\n";                                         		
        up_r_file << "+";                                          		
        up_r_file << "\n";                                         		
        up_r_file << r_qual_orig;                                  		
        up_r_file << "\n";
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
        
		p_r_file << r_id;
        p_r_file << "\n";
        p_r_file << r_seq;
        p_r_file << "\n";
        p_r_file << "+";
        p_r_file << "\n";
        p_r_file << r_qual;
        p_r_file << "\n";
	  }
      
	  std::getline(file_1,f_id);
      std::getline(file_2,r_id); 	  
    }
  }
  else std::cout << "Unable to open file\n";
	 
  file_1.close();														
  up_f_file.close();													
  p_f_file.close();
  file_2.close();														
  up_r_file.close();													
  p_r_file.close();  
	  
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  
  int t_read_count = f_read_count + r_read_count;
    
  double f_lost_percent = double(read_lost)*100/double(f_read_count);
  double r_lost_percent = double(read_lost)*100/double(r_read_count);
  double t_lost_percent = double(read_lost)*2*100/double(t_read_count);
  
  int f_trim = f_read_count - read_lost;
  int r_trim = r_read_count - read_lost;
  int t_trim = t_read_count - read_lost*2;
  double f_trim_percent = double(f_trim)*100/double(f_read_count);
  double r_trim_percent = double(r_trim)*100/double(r_read_count);
  double t_trim_percent = double(t_trim)*100/double(t_read_count);
  
  
  int t_adap_trim = f_adap_trim + r_adap_trim;
  double f_adap_trim_percent = double(f_adap_trim)*100/double(f_read_count);
  double r_adap_trim_percent = double(r_adap_trim)*100/double(r_read_count);
  double t_adap_trim_percent = double(t_adap_trim)*100/double(t_read_count);
  
  int t_polyG_trim = f_polyG_trim + r_polyG_trim;
  double f_polyG_trim_percent = double(f_polyG_trim)*100/double(f_read_count);
  double r_polyG_trim_percent = double(r_polyG_trim)*100/double(r_read_count);
  double t_polyG_trim_percent = double(t_polyG_trim)*100/double(t_read_count);
  
  int t_qual_trim = f_qual_trim + r_qual_trim;
  double f_qual_trim_percent = double(f_qual_trim)*100/double(f_read_count);
  double r_qual_trim_percent = double(r_qual_trim)*100/double(r_read_count);
  double t_qual_trim_percent = double(t_qual_trim)*100/double(t_read_count);
  
  printf("Trimming Statistics      | %10s\t%10s\t\t\t%10s\t%10s\t\t\t%10s\t%10s\n", "# Forward", "% Forward", "# Reverse", "% Reverse", "# Total", "% Total");
  printf("Reads Before Trimming    | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", f_read_count, 100.00, 				r_read_count, 	100.00, 				t_read_count,	100.00);
  printf("Reads After Trimming     | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", f_trim, 		f_trim_percent, 		r_trim, 		r_trim_percent, 		t_trim, 		t_trim_percent);
  printf("Reads Lost               | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", read_lost, 	f_lost_percent, 		read_lost, 		r_lost_percent, 		read_lost*2, 	t_lost_percent);
  printf("Adapter Trimmed Reads    | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", f_adap_trim, 	f_adap_trim_percent, 	r_adap_trim, 	r_adap_trim_percent, 	t_adap_trim, 	t_adap_trim_percent);
  printf("polyG Trimmed Reads      | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", f_polyG_trim, f_polyG_trim_percent, 	r_polyG_trim, 	r_polyG_trim_percent, 	t_polyG_trim, 	t_polyG_trim_percent);
  printf("Quality Trimmed Reads    | %10i\t%10.2f\t\t\t%10i\t%10.2f\t\t\t%10i\t%10.2f\n", f_qual_trim, 	f_qual_trim_percent, 	r_qual_trim, 	r_qual_trim_percent, 	t_qual_trim, 	t_qual_trim_percent);
  std::cout << t_read_count << " Reads Analyzed in " << duration.count() << " seconds\n";
}