/*

Program that compares lines from ancestor and endpoint common lines file and estimates mut rate.
Note that these input files are preprocessed and derived from pileups. (ATCG data)
We have two files with the same number of lines and the two files correspond line by line.

This program takes 3 argument: extant ancestor_common_processed_pileup endpoint_common_processed_pileup

	extant: for defining minimum mut_freq to consider endpoint as muted
	
	extant must be in {8, 12, 16, 20, 24, 28, 32, 48, 75, 100}
	
This program, for all couples of lines in the two files:

1) Reads one line from ancestor file and one from endpoint; 
2) Parses them to get muted frequency and coverage;
3) Applies thresholds and in case skips lines;
4) Checks if lines that are fine with thresholds are muted or not. 

Tresholds:

	- For both files discards lines with reference base different from A, T, C, G
	
	- For ancestor discard lines with
		 coverage < 100
		 mut_freq > 0.0
	
	- For endpoint discard lines with:
		mut_freq > 0.2 
		unbalanced forward and reverse muted reads.



Counts base as muted if endpoint frequency 

	f >= Cutoff(Coverage, f_min) e Coverage = TotCoverage/3

	cutoff (Coverage, f_min) = f_min + alpha/Sqrt[Coverage]

based on extant f_min and corresponding alphas are

	f_min = {1/8, 1/12, 1/16, 1/20, 1/24, 1/28, 1/32, 1/48, 1/75, 1/100}

	alpha = {0.52, 0.45, 0.36, 0.32, 0.3, 0.28, 0.25, 0.2, 0.18, 0.15} 
*/


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm> // std::max(a,b)

//Boost 
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

double compute_mut_rate(double P_hat, double death_prob, int extant);

bool endpoint_muted(float mut_freq, int total_coverage, int extant);

int main(int argc, char** argv){

	//*******************************************
	//          COMMAND LINE ARGS
	//*******************************************
	
	if (argc != 4){
		cerr << "Usage: " << argv[0] << " extant ancestor_common_lines_file endpoint_common_lines_file" << endl;
		return -1;
	}
	
	
	int extant = atoi(argv[1]);
	
	string ancestor_file = argv[2];
	string endpoint_file = argv[3];
	
	
	if (extant != 8 && extant != 12 &&extant != 16 && extant != 20 && extant != 24 && extant != 28 && extant != 32 && extant != 48 && extant != 75 && extant != 100){
		cerr << "Extant must be of following: {8, 12, 16, 20, 24, 28, 32, 48, 75, 100}" 
		<< endl;
		return -1;
	}
	
	
	//*******************************************
	//          FIXED THRESHOLDS 
	//*******************************************
	//ancestor
	int coverage_min_ancestor = 100;
	float freq_max_ancestor_confirm_reference = 0.;
	//endpoint
	float endpoint_freq_max = 0.2;
	
	//*******************************************
	//          COMPRESSED INPUT
	//*******************************************
	

	//open and unzip ancestor file
	ifstream file_1(ancestor_file, ios_base::in | ios_base::binary);
	
	if (file_1.is_open() == false) {
		cerr << "Ancestor not opened: exit!" <<endl;
		return -1;
	}
	
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> file_1_inbuf;
    	file_1_inbuf.push(boost::iostreams::gzip_decompressor());
    	file_1_inbuf.push(file_1);
    	//Convert streambuf to istream
    	istream file_1_instream(&file_1_inbuf);
	
	cout << endl << "Ancestor file opened." << endl;
	
	//open and unzip endpoint file
	ifstream file_2(endpoint_file, ios_base::in | ios_base::binary);
	
	if (file_2.is_open() == false) {
		cerr << "Endopint file not opened: exit!" <<endl;
		return -1;
	}
	
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> file_2_inbuf;
    	file_2_inbuf.push(boost::iostreams::gzip_decompressor());
    	file_2_inbuf.push(file_2);
    	//Convert streambuf to istream
    	istream file_2_instream(&file_2_inbuf);
	
	cout << "Endpoint file opened." << endl << endl;
	
	//*******************************************
	//          PARSING VARIABLES
	//*******************************************
	//Ancestor parsing variables-------------
	string ancestor_line; //line of the file to parse
	string ancestor_chromosome; //pos
	int ancestor_chromosome_number = 0;
	long int ancestor_base_number = 0; //1-based
	string ancestor_ATCG_data; 
	char ancestor_reference;
	int ancestor_A_counter = 0;
	int ancestor_T_counter = 0;
	int ancestor_C_counter = 0;
	int ancestor_G_counter = 0;
	int ancestor_a_counter = 0;
	int ancestor_t_counter = 0;
	int ancestor_c_counter = 0;
	int ancestor_g_counter = 0;
	int ancestor_coverage = 0;
	double ancestor_mut_freq = 0;
	
	//Endpoint parsing variables-------------
	string endpoint_line; //line of the file to parse
	string endpoint_chromosome; 
	int endpoint_chromosome_number = 0;
	long int endpoint_base_number = 0; 
	string endpoint_ATCG_data; 
	char endpoint_reference;
	int endpoint_A_counter = 0;
	int endpoint_T_counter = 0;
	int endpoint_C_counter = 0;
	int endpoint_G_counter = 0;
	int endpoint_a_counter = 0;
	int endpoint_t_counter = 0;
	int endpoint_c_counter = 0;
	int endpoint_g_counter = 0;
	int endpoint_coverage = 0;
	double endpoint_mut_freq = 0;
	char endpoint_muted_in;
	
	//*******************************************
	//          OTHERS 
	//*******************************************
	long int count_line = 0; //all lines in file
	long int discarded_lines = 0; //not matching thresholds
	long int muted_bases = 0;
	long int not_muted_bases = 0;
	double P_hat;
	
	//******************************************
	//--------------------------------------
	// 		    START
	// 
	// read line by line both ancestor and endpoint
	// check thresholds and count muted lines
	//--------------------------------------
	//******************************************
	
	cout << "Start reading files." << endl << endl;
	
	//Loop for all line in file 2 (endpoint)
	while (getline(file_2_instream, endpoint_line)) { //until file_2 EOF
		
		//reset counters
		endpoint_coverage = 0;
		ancestor_coverage = 0;
		
		//*******************************************
		// READ AND PARSE A LINE FROM ENDPOINT FILE
		//*******************************************
		
		//Parse file 2 input line
		stringstream file_2_linestream(endpoint_line);
		//get file_2 line pos finding the separator ('\t')
		getline(file_2_linestream, endpoint_chromosome, '\t');
		//get file_2 loc
		file_2_linestream >> endpoint_base_number;
		//we need chromosome numeber and loc
		//extract chromosome number from chromosome string
		// "chrN", N integer
		//     ^        
		endpoint_chromosome_number = atoi(&endpoint_chromosome[3]);
		//get file_2 ATCG data finding the separator ('\n')
		getline(file_2_linestream, endpoint_ATCG_data, '\n');
		
		//parse endpoint ATCG data
		stringstream endpoint_ATCG_linestream(endpoint_ATCG_data);
		//endpoint reference
		endpoint_ATCG_linestream >> endpoint_reference;
		//counters
		endpoint_ATCG_linestream >> endpoint_A_counter; endpoint_coverage += endpoint_A_counter;
		endpoint_ATCG_linestream >> endpoint_T_counter; endpoint_coverage += endpoint_T_counter;
		endpoint_ATCG_linestream >> endpoint_C_counter; endpoint_coverage += endpoint_C_counter;
		endpoint_ATCG_linestream >> endpoint_G_counter; endpoint_coverage += endpoint_G_counter;
	
		endpoint_ATCG_linestream >> endpoint_a_counter; endpoint_coverage += endpoint_a_counter;
		endpoint_ATCG_linestream >> endpoint_t_counter; endpoint_coverage += endpoint_t_counter;
		endpoint_ATCG_linestream >> endpoint_c_counter; endpoint_coverage += endpoint_c_counter;
		endpoint_ATCG_linestream >> endpoint_g_counter; endpoint_coverage += endpoint_g_counter;
		
		//**********************************************
		// COMPUTE MUTED READS AND MUT FREQ FOR ENDPOINT
		//**********************************************
		int A = endpoint_A_counter + endpoint_a_counter;
		int T = endpoint_T_counter + endpoint_t_counter;
		int C = endpoint_C_counter + endpoint_c_counter;
		int G = endpoint_G_counter + endpoint_g_counter;
		
		if (endpoint_reference == 'A'){
			int muted_reads = max(T,C);
			muted_reads = max(muted_reads,G);
			endpoint_mut_freq = double(muted_reads)/double(endpoint_coverage);
			if (muted_reads != 0) {
				if (muted_reads == T) endpoint_muted_in = 'T';
				if (muted_reads == C) endpoint_muted_in = 'C';
				if (muted_reads == G) endpoint_muted_in = 'G';
			} else endpoint_muted_in = ' ';
		}
		if (endpoint_reference == 'T'){
			int muted_reads = max(A,C);
			muted_reads = max(muted_reads,G);
			endpoint_mut_freq = double(muted_reads)/double(endpoint_coverage);
			if (muted_reads != 0) {
				if (muted_reads == A) endpoint_muted_in = 'A';
				if (muted_reads == C) endpoint_muted_in = 'C';
				if (muted_reads == G) endpoint_muted_in = 'G';
			} else endpoint_muted_in = ' ';
		}
		if (endpoint_reference == 'C'){
			int muted_reads = max(A,T);
			muted_reads = max(muted_reads,G);
			endpoint_mut_freq = double(muted_reads)/double(endpoint_coverage);
			if (muted_reads != 0) {
				if (muted_reads == A) endpoint_muted_in = 'A';
				if (muted_reads == T) endpoint_muted_in = 'T';
				if (muted_reads == G) endpoint_muted_in = 'G';
			} else endpoint_muted_in = ' ';
		}
		if (endpoint_reference == 'G'){
			int muted_reads = max(A,T);
			muted_reads = max(muted_reads,C);
			endpoint_mut_freq = double(muted_reads)/double(endpoint_coverage);
			if (muted_reads != 0) {
			if (muted_reads == A) endpoint_muted_in = 'A';
			if (muted_reads == T) endpoint_muted_in = 'T';
			if (muted_reads == C) endpoint_muted_in = 'C';
			} else endpoint_muted_in = ' ';
		}
		
		//*******************************************
		// READ AND PARSE A LINE FROM ENDPOINT FILE
		//*******************************************
		//take line from file 1
		getline(file_1_instream, ancestor_line);
		stringstream file_1_linestream(ancestor_line);
		getline(file_1_linestream, ancestor_chromosome, '\t');
		// "chrN", N integer
		//     ^        
		ancestor_chromosome_number = atoi(&ancestor_chromosome[3]);
		file_1_linestream >> ancestor_base_number;
		getline(file_1_linestream, ancestor_ATCG_data, '\n'); 
		
		//parse ancestor ATCG data
		stringstream ancestor_ATCG_linestream(ancestor_ATCG_data);
		//ancestor reference
		ancestor_ATCG_linestream >> ancestor_reference;
		//counters
		ancestor_ATCG_linestream >> ancestor_A_counter; ancestor_coverage += ancestor_A_counter;
		ancestor_ATCG_linestream >> ancestor_T_counter; ancestor_coverage += ancestor_T_counter;
		ancestor_ATCG_linestream >> ancestor_C_counter; ancestor_coverage += ancestor_C_counter;
		ancestor_ATCG_linestream >> ancestor_G_counter; ancestor_coverage += ancestor_G_counter;
	
		ancestor_ATCG_linestream >> ancestor_a_counter; ancestor_coverage += ancestor_a_counter;
		ancestor_ATCG_linestream >> ancestor_t_counter; ancestor_coverage += ancestor_t_counter;
		ancestor_ATCG_linestream >> ancestor_c_counter; ancestor_coverage += ancestor_c_counter;
		ancestor_ATCG_linestream >> ancestor_g_counter; ancestor_coverage += ancestor_g_counter;
		
		//**********************************************
		// COMPUTE MUTED READS AND MUT FREQ FOR ANCESTOR
		//**********************************************
		A = ancestor_A_counter + ancestor_a_counter;
		T = ancestor_T_counter + ancestor_t_counter;
		C = ancestor_C_counter + ancestor_c_counter;
		G = ancestor_G_counter + ancestor_g_counter;
		
		if (ancestor_reference == 'A'){
			int muted_reads = max(T,C);
			muted_reads = max(muted_reads,G);
			ancestor_mut_freq = double(muted_reads)/double(ancestor_coverage);
		}
		if (ancestor_reference == 'T'){
			int muted_reads = max(A,C);
			muted_reads = max(muted_reads,G);
			ancestor_mut_freq = double(muted_reads)/double(ancestor_coverage);
		}
		if (ancestor_reference == 'C'){
			int muted_reads = max(A,T);
			muted_reads = max(muted_reads,G);
			ancestor_mut_freq = double(muted_reads)/double(ancestor_coverage);
		}
		if (ancestor_reference == 'G'){
			int muted_reads = max(A,T);
			muted_reads = max(muted_reads,C);
			ancestor_mut_freq = double(muted_reads)/double(ancestor_coverage);
		}

		
		count_line++; //all lines in files
		
		//**********************************************
		// APPLY THRESHOLDS (and do checks)
		//**********************************************
		
		//rarely there are "M","N", or "R" in reference
		//we need to discard that line
		if(ancestor_reference != 'A' && ancestor_reference != 'T' 
		&& ancestor_reference != 'C' && ancestor_reference != 'G') {
			//reference is not A, T, C, G
			discarded_lines++;
			continue; //skip this line
		}
		
		if(endpoint_reference != 'A' && endpoint_reference != 'T' 
		&& endpoint_reference != 'C' && endpoint_reference != 'G') {
			//reference is not A, T, C, G
			discarded_lines++;
			continue; //skip this line
		}
		
		// check that we loaded corresponding lines
		// (common lines file are ok!)
		assert(endpoint_chromosome_number == ancestor_chromosome_number);
		assert(endpoint_base_number == ancestor_base_number);
		assert(endpoint_reference == ancestor_reference);
		
		// apply fixed threshold on ancestor's coverage min
		if (ancestor_coverage < coverage_min_ancestor) {
			discarded_lines++;
			continue; //skip this line
		}
		
		// apply theshold on ancestor mut_freq
		//same as check if == zero 
		if(ancestor_mut_freq > freq_max_ancestor_confirm_reference){
			discarded_lines++;
			continue; //skip this line
		}
		
		// apply theshold on endpoint mut_freq max
		if(endpoint_mut_freq > endpoint_freq_max) {
			discarded_lines++;
			continue; //skip this line
		}
		
		//are forward and reverse reads balanced?
		// example: endpoint muted in A
		// accept line if 
		// (1/2 - 1/sqrt(#(A + a))) < #(A) / #(A + a) < (1/2 + 1/sqrt(#(A + a)))
		
		if(endpoint_muted_in == 'A') {
			float forward_reads = float(endpoint_A_counter);
			float total = forward_reads + float(endpoint_a_counter);
			float x = forward_reads / total ;
			
			if ( (x < 0.5 - 1./sqrt(total)) || (x > 0.5 + 1./sqrt(total)) ){
				discarded_lines++;
				continue; //skip this line
			}
		}
			
		if(endpoint_muted_in == 'T') {
			float forward_reads = float(endpoint_T_counter);
			float total = forward_reads + float(endpoint_t_counter);
			float x = forward_reads / total ;
			
			if ( (x < 0.5 - 1./sqrt(total)) || (x > 0.5 + 1./sqrt(total)) ){
				discarded_lines++;
				continue; //skip this line
			}
		}	
		if(endpoint_muted_in == 'C') {
			float forward_reads = float(endpoint_C_counter);
			float total = forward_reads + float(endpoint_c_counter);
			float x = forward_reads / total ;
			
			if ( (x < 0.5 - 1./sqrt(total)) || (x > 0.5 + 1./sqrt(total)) ){
				discarded_lines++;
				continue; //skip this line
			}
		}
		
		if(endpoint_muted_in == 'G') {
			float forward_reads = float(endpoint_G_counter);
			float total = forward_reads + float(endpoint_g_counter);
			float x = forward_reads / total ;
			
			if ( (x < 0.5 - 1./sqrt(total)) || (x > 0.5 + 1./sqrt(total)) ){
				discarded_lines++;
				continue; //skip this line
			}
		}
		
		//**********************************************
		// DECIDE IF MUTED OR NOT
		//**********************************************
		// finally if here theese lines are fine with all thresholds
		//check if to be considered muted or not 
		
		if( endpoint_muted(endpoint_mut_freq, endpoint_coverage, extant)){
			//muted -> Numerator & denominator
			muted_bases++; 
		} else {
			//not muted -> Denominator
			not_muted_bases++;
		}
	} //EOF
	
	//Print result summary
	
	cout << "File readings ended: total lines " << count_line << endl;
	cout << "Discarded: " << discarded_lines << endl;
	cout << "Parameter used extant = " << extant << endl;
	
	assert(count_line - discarded_lines == muted_bases + not_muted_bases);
	
	// compute muted/(muted + not_muted)	
	cout << "Muted bases:" << muted_bases << endl;
	cout << "Not muted bases: " << not_muted_bases << endl;
	
	P_hat = double (muted_bases) / double (muted_bases + not_muted_bases);
	
	cout << "muted/(muted + not_muted) = " << P_hat <<endl; 
	
	cout << "Mutatior rate min [d/(b+d) = 0.45] = " << compute_mut_rate(P_hat, 0.45, extant) << endl;
	
	cout << "Mutation rate max [d/(b+d) = 0.1] = " << compute_mut_rate(P_hat, 0.1, extant) << endl;
	
	//Cleanup
	file_1.close();
	file_2.close();
	
	return 0;
} // main
		
		
double compute_mut_rate(double P_hat, double death_prob, int extant){

	// gen = log_{2*(1-death_prob)} extant
	float generations = log(extant)/log(2.*(1.-death_prob));
	
	//integral estimate with continuous time
	float attempts = (1. - death_prob*death_prob) * (
				( pow((2*(1-death_prob)), generations - 1.) - 1.) /
			  	log(2.*(1.- death_prob))
			  	) + extant;
			  			  	
	double mut_rate = -1. *( log(1. - P_hat)/attempts);
	
	return mut_rate;
}

bool endpoint_muted(float mut_freq, int total_coverage, int extant){
/*
Count base as muted if endpoint frequency 

	f >= Cutoff(Coverage, f_min) e Coverage = TotCoverage/3

	cutoff (Coverage, f_min) = f_min + alpha/Sqrt[Coverage]

based on extant f_min and corresponding alphas are

	f_min = {1/8, 1/12, 1/16, 1/20, 1/24, 1/28, 1/32, 1/48, 1/75, 1/100}

	alpha = {0.52, 0.45, 0.36, 0.32, 0.3, 0.28, 0.25, 0.2, 0.18, 0.15} 
*/
	float f_min= 0.;
	float coverage = float(total_coverage)/3.;
	double alpha = 0.;
	
	if (extant == 8){
		f_min = 1./extant;
		alpha = 0.52;
	}
	if (extant == 12){
		f_min = 1./extant;
		alpha = 0.45;
	}
	
	if (extant == 16){
		f_min = 1./extant;
		alpha = 0.36;
	}
	if (extant == 20){
		f_min = 1./extant;
		alpha = 0.32;
	}
	if (extant == 24){
		f_min = 1./extant;
		alpha = 0.3;
	}
	if (extant == 28){
		f_min = 1./extant;
		alpha = 0.28;
	}
	if (extant == 32){
		f_min = 1./extant;
		alpha = 0.25;
	}
	if (extant == 48){
		f_min = 1./extant;
		alpha = 0.2;
	}
	
	if (extant == 75){
		f_min = 1./extant;
		alpha = 0.18;
	}
	if (extant == 100){
		f_min = 1./extant;
		alpha = 0.15;
	}
	
	if ( mut_freq >= f_min + alpha/sqrt(coverage) ) {
		//muted
		return true;
	} else {
		//not muted
		return false;
	}
}
