#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>

//Boost 
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

//Function to compute mut_freq from ATCG data for a line
double compute_mut_freq(string ATCG_data);

//Function to compute mut_freq from ATCG data for a line
int compute_coverage(string ATCG_data);



int main(int argc, char** argv){

	//*******************************************
	//          COMMAND LINE ARGS
	//*******************************************

	if (argc != 5){
		cerr << "Usage: " << argv[0] << " processed_pileup_ancestor_1 processed_pileup_ancestor_2" 
		<< " common_ancestor_file_1 common_ancestor_file_2" << endl;
		return -1;
	}
	//filenames args
	string pseudo_pileup_1 = argv[1];
	string pseudo_pileup_2 = argv[2];

	string common_ancesestor_1 = argv[3];
	string common_ancesestor_2 = argv[4];
	
	//*******************************************
	//          COMPRESSED INPUT
	//*******************************************
	
	//(IN) FILE 1 --------------------------------------------------
	//assume it's compressed (.gz)
    	ifstream file_1(pseudo_pileup_1, ios_base::in | ios_base::binary);
    	
    	if (file_1.is_open() == false) {
		cerr << "File 1 not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> file_1_inbuf;
    	file_1_inbuf.push(boost::iostreams::gzip_decompressor());
    	file_1_inbuf.push(file_1);
    	//Convert streambuf to istream
    	istream file_1_instream(&file_1_inbuf);
	
	cerr << endl << "File 1 opened." << endl;

	//(IN) FILE 2 --------------------------------------------------
	//assume it's compressed (.gz)
    	ifstream file_2(pseudo_pileup_2, ios_base::in | ios_base::binary);
    	
    	if (file_2.is_open() == false) {
		cerr << "File 2 file not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> file_2_inbuf;
    	file_2_inbuf.push(boost::iostreams::gzip_decompressor());
    	file_2_inbuf.push(file_2);
    	//Convert streambuf to istream
    	istream file_2_instream(&file_2_inbuf);
	
	cerr << "File 2 opened." << endl << endl;
	
	//*******************************************
	//          COMPRESSED OUTPUT
	//*******************************************
	// Out for file 1
	ofstream out_file_1(common_ancesestor_1, ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_1;
	outbuf_1.push(boost::iostreams::gzip_compressor());
	outbuf_1.push(out_file_1);
	//Convert streambuf to ostream
	ostream out_1(&outbuf_1);

	// Out for file 2
	ofstream out_file_2(common_ancesestor_2, ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_2;
	outbuf_2.push(boost::iostreams::gzip_compressor());
	outbuf_2.push(out_file_2);
	//Convert streambuf to ostream
	ostream out_2(&outbuf_2);
	
	//*******************************************
	//          PARSING VARIABLES
	//*******************************************
	//File 1 parsing variables-------------
	string file_1_line; //line of the file to parse
	string file_1_chromosome; //pos
	int file_1_chromosome_number = 0;
	long int file_1_base_number = 0; //1-based
	string file_1_ATCG_data; //not used here
	
	//File 2 parsing variables-------------
	string file_2_line; //line of the file to parse
	string file_2_chromosome; 
	int file_2_chromosome_number = 0;
	long int file_2_base_number = 0; 
	string file_2_ATCG_data; //not used here
	
	//Useful variables
	long int file_1_line_counter= 0;
	long int file_2_line_counter= 0;
	long int found_lines = 0;
	long int not_found_lines = 0;
	int current_chr = 0;
	long int found_in_chr = 0;
	bool file_1_ended = false;
	
	double mut_freq_ancestor_1;
	double mut_freq_ancestor_2;
	int coverage;
	
	
	//******************************************
	//--------------------------------------
	// 		    START
	// 
	// for all lines in file 2 search the 
	// corrisponding line in file 1
	//--------------------------------------
	//******************************************
	cerr<< "Start reading files." << endl << endl;
	
	//read first line of file 1
	getline(file_1_instream, file_1_line);
	file_1_line_counter++;
	
	stringstream file_1_linestream(file_1_line);
	//get file_1 pos finding the separator ('\t')
	getline(file_1_linestream, file_1_chromosome, '\t');
	// "chrN", N integer
	//     ^        
	file_1_chromosome_number = atoi(&file_1_chromosome[3]);
	//get file_1 loc
	file_1_linestream >> file_1_base_number;
	//get file_1 ATCG data finding the separator ('\n')
	getline(file_1_linestream, file_1_ATCG_data, '\n');
	
	//Loop for all line in file 2
	while (getline(file_2_instream, file_2_line)) { //until file_2 EOF
		
		//Parse file 2 input line
		stringstream file_2_linestream(file_2_line);
		file_2_line_counter++;
		//get file_2 line pos finding the separator ('\t')
		getline(file_2_linestream, file_2_chromosome, '\t');
		//get file_2 loc
		file_2_linestream >> file_2_base_number;
		//we need chromosome numeber and loc
		//extract chromosome number from chromosome string
		// "chrN", N integer
		//     ^        
		file_2_chromosome_number = atoi(&file_2_chromosome[3]);
		//get file_2 ATCG data finding the separator ('\n')
		getline(file_2_linestream, file_2_ATCG_data, '\n');
		
		//cout<< "File 2 line red: " <<file_2_chromosome << " " << file_2_base_number << " "<< file_2_ATCG_data<< endl;
		
		//Add line to file 2 histogram
		mut_freq_ancestor_2 = compute_mut_freq(file_2_ATCG_data);
		
		if(mut_freq_ancestor_2 < 0.)  {
			//rarely there are "M","N", or "R" in reference
			//we need to skip that line
			cerr << "Freq error! Skip this line from file 2:"<< endl
			<< file_2_chromosome << "\t" << file_2_base_number << file_2_ATCG_data<< endl;
		} 
		
		//Go on with file 1 until we find or exceed the line 
		
		while((file_1_chromosome_number < file_2_chromosome_number) 
			|| (
				(file_1_chromosome_number == file_2_chromosome_number)
				&& (file_1_base_number < file_2_base_number) 
			   )
			){
			
			//cout <<"Discarded from file 1 "<<file_1_chromosome << " " << file_1_base_number << " "<< file_1_ATCG_data<< endl;
		     
			
			// read next line on file_1
			if(getline(file_1_instream, file_1_line)){//if not file 1 EOF
				file_1_line_counter++;
	
				stringstream file_1_linestream(file_1_line);
				//get file_1 pos finding the separator ('\t')
				getline(file_1_linestream, file_1_chromosome, '\t');
	
				// "chrN", N integer
				//     ^        
				file_1_chromosome_number = atoi(&file_1_chromosome[3]);
	
				//get file_1 loc
				file_1_linestream >> file_1_base_number;
	
				//get file_1 ATCG data finding the separator ('\n')
				getline(file_1_linestream, file_1_ATCG_data, '\n'); 
			} //if not file 1 EOF
			else{ // file 1 finished before file 2
				if (!file_1_ended){
					cout << "File 1 finished before file 2." <<endl;
					cout << "Last line was: " << file_1_chromosome << " " << file_1_base_number << endl;
					file_1_ended = true;
				}
				break; //break file 1 loop
				// we don't break also file 2 loop just to count it's lines
			}
		
		} //file 1 loop
		
		//cout<< "Last line after file 1 loop" << endl;
		//cout<<file_1_chromosome << " " << file_1_base_number << " "<< file_1_ATCG_data<< endl;
		
		//count common lines per chromosome
		if (file_2_chromosome_number != current_chr){
			//check if chromosomes are ordered
			assert(current_chr < file_2_chromosome_number);
			
			if (current_chr != 0) {
				cout << "---> Found in chr" << current_chr <<" "<< found_in_chr << " common lines"<< endl;
				found_in_chr = 0;
			}
			current_chr = file_2_chromosome_number;
		}
			
		if ((file_1_chromosome_number == file_2_chromosome_number)
			&& (file_1_base_number == file_2_base_number)) {
			//cout << "Found line!" << endl;
			// We found a common line!			
			mut_freq_ancestor_2 = compute_mut_freq(file_2_ATCG_data);
			mut_freq_ancestor_1 = compute_mut_freq(file_1_ATCG_data);
			
			//if both ancestor confirm the reference (zero muted lectures)
			if(mut_freq_ancestor_1 == 0. && mut_freq_ancestor_2 == 0.)  {
				
				found_lines++;
				found_in_chr++;
				
				//output relative lines to outfiles 
				//formatted correcly
				out_1 << file_1_chromosome << "\t" 
				<< file_1_base_number << file_1_ATCG_data << endl;
				
				out_2 << file_2_chromosome << "\t" 
				<< file_2_base_number << file_2_ATCG_data << endl;
			}
			
		} else {
			//line not found
			
			//cout << "Not common line: " <<file_1_chromosome << " " << file_1_base_number << " "<< file_1_ATCG_data<< endl; 
			//cout<< "Looking for: " <<file_2_chromosome << " " << file_2_base_number << endl<<endl;
			if (!file_1_ended){
				// the last line red from file 1
				// is after the one we are looking for
				assert( ( (file_1_chromosome_number == file_2_chromosome_number)
					   && (file_1_base_number > file_2_base_number) )
					 || (file_1_chromosome_number > file_2_chromosome_number) );
			}
			not_found_lines++;
			//cout << "Not found line!" << endl; 
		}
		
	} //file 2 loop, until file 2 EOF
	
	//print last chromosome counter
	cout << "---> Found in chr" << current_chr <<" "<< found_in_chr << " common lines"<< endl;
	cout << endl;
		
	//print summary 
	cout << "Terminated: " << endl << "total common lines found " << found_lines << " of " 
	<< file_2_line_counter << " red from file 2." << endl;
	
	cout << "Lines not found " << not_found_lines << endl;
	
	if (file_1_ended) {
		cout << "File 1 ended before file 2 was totally red." << endl;
		cout << "File 1 total lines: " << file_1_line_counter << endl;
	} else {
		cout << "File 2 was totally red, but file 1 was not." << endl;
		cout << "File 1 lines red: " << file_1_line_counter << endl;
	}
	
	
	//*******************************************
	//            CLEANUP
	//*******************************************
	//close input files 
	file_1.close();
	file_2.close();
	
	//close output
	boost::iostreams::close(outbuf_1); 
	boost::iostreams::close(outbuf_2);
    	out_file_1.close();
    	out_file_2.close();
    	cout << "Otput written on files."<<endl;
	
	//end
	return 0;
} //main() end


//Function to compute mut_freq from ATCG data for a line
double compute_mut_freq(string ATCG_data){
	char reference;
	int A_counter;
	int T_counter;
	int C_counter;
	int G_counter;
	int a_counter;
	int t_counter;
	int c_counter;
	int g_counter;
	
	stringstream ATCG_linestream(ATCG_data);
	
	int coverage = 0;
	double mut_freq;
	
	ATCG_linestream >> reference;
	//rarely there are "M","N", or "R" in reference
	//we need to discard that line
	if(reference != 'A' && reference != 'T' 
		&& reference != 'C' && reference != 'G') {
		//reference is not A, T, C, G
		// return -1. as impossible frequency
		return -1.;
	}
	ATCG_linestream >> A_counter; coverage += A_counter;
	ATCG_linestream >> T_counter; coverage += T_counter;
	ATCG_linestream >> C_counter; coverage += C_counter;
	ATCG_linestream >> G_counter; coverage += G_counter;
	
	ATCG_linestream >> a_counter; coverage += a_counter;
	ATCG_linestream >> t_counter; coverage += t_counter;
	ATCG_linestream >> c_counter; coverage += c_counter;
	ATCG_linestream >> g_counter; coverage += g_counter;
	
	if (reference == 'A'){
		int muted_reads = coverage - ( A_counter + a_counter);
		mut_freq = double(muted_reads)/double(coverage);
		}
	if (reference == 'T'){
		int muted_reads = coverage - ( T_counter + t_counter);
		mut_freq = double(muted_reads)/double(coverage);
		}
	if (reference == 'C'){
		int muted_reads = coverage - ( C_counter + c_counter);
		mut_freq = double(muted_reads)/double(coverage);
		}
	if (reference == 'G'){
		int muted_reads = coverage - ( G_counter + g_counter);
		mut_freq = double(muted_reads)/double(coverage);
		}
		
	return mut_freq;
}

//Function to compute mut_freq from ATCG data for a line
int compute_coverage(string ATCG_data){
	char reference;
	int A_counter;
	int T_counter;
	int C_counter;
	int G_counter;
	int a_counter;
	int t_counter;
	int c_counter;
	int g_counter;
	
	stringstream ATCG_linestream(ATCG_data);
	
	int coverage = 0;
	double mut_freq;
	
	ATCG_linestream >> reference;
	
	ATCG_linestream >> A_counter; coverage += A_counter;
	ATCG_linestream >> T_counter; coverage += T_counter;
	ATCG_linestream >> C_counter; coverage += C_counter;
	ATCG_linestream >> G_counter; coverage += G_counter;
	
	ATCG_linestream >> a_counter; coverage += a_counter;
	ATCG_linestream >> t_counter; coverage += t_counter;
	ATCG_linestream >> c_counter; coverage += c_counter;
	ATCG_linestream >> g_counter; coverage += g_counter;
		
	return coverage;
}

