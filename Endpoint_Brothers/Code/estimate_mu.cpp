
#include <string>
#include <sstream>
#include <cassert>
#include <algorithm> // std::max(a,b)

//Boost 
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

//Root
#include "TApplication.h"
#include "TComplex.h"
#include <TGraph.h>
#include "TAxis.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "THStack.h"
#include "TLegend.h"
#include "TH2D.h"


using namespace std;

char muted_in(string ATCG_data);

double compute_mut_freq(string ATCG_data, char muted_in);

int compute_coverage(string ATCG_data);

double compute_mut_rate(double P_hat, double death_prob, int extant);

bool endpoint_muted(double mut_freq, int total_coverage, int extant);

/*
Count base as muted if endpoint frequency 

	f >= Cutoff(Coverage, f_min) e Coverage = TotCoverage/3

	cutoff (Coverage, f_min) = f_min + alpha/Sqrt[Coverage]

based on extant f_min and corresponding alphas are

	f_min = {1/8, 1/12, 1/16, 1/20, 1/24, 1/28, 1/32, 1/48, 1/75, 1/100}

	alpha = {0.52, 0.45, 0.36, 0.32, 0.3, 0.28, 0.25, 0.2, 0.18, 0.15} 
	
	extant: for defining minimum mut_freq to consider endpoint as muted
	
	extant must be in {8, 12, 16, 20, 24, 28, 32, 48, 75, 100}
			    |                                   |
			    |-----------10 elements-------------|
*/

int main(){

	//*******************************************
	//          THRESHOLDS
	//*******************************************
	int ancestor_coverage_min = 30;
	long int discarded_by_threshold = 0;
	double endpoint_max_mut_freq = 0.2;

	//*******************************************
	//          FILENAMES
	//*******************************************
	string ancestor_1 = "CRC1307-02-0.all_common.gz"; 
	string ancestor_2 = "CRC1307-09-0.all_common.gz";
	
	string endpoint_1 = "CRC1307-02-1-E.all_common.gz";
	string endpoint_2 = "CRC1307-09-1-B.all_common.gz";
	
	
	//*******************************************
	//          COMPRESSED INPUT
	//*******************************************
	//(IN) ANCESTOR FILE 1 --------------------------------------------------
	ifstream file_ancestor_1(ancestor_1, ios_base::in | ios_base::binary);
    	
    	if (file_ancestor_1.is_open() == false) {
		cerr << ancestor_1 << " not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> ancestor_1_inbuf;
    	ancestor_1_inbuf.push(boost::iostreams::gzip_decompressor());
    	ancestor_1_inbuf.push(file_ancestor_1);
    	//Convert streambuf to istream
    	istream ancestor_1_instream(&ancestor_1_inbuf);
	
	cerr << endl << ancestor_1 << " opened." << endl;
	
	//(IN) ANCESTOR FILE 2 --------------------------------------------------
	ifstream file_ancestor_2(ancestor_2, ios_base::in | ios_base::binary);
    	
    	    	if (file_ancestor_2.is_open() == false) {
		cerr << ancestor_2 << " not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> ancestor_2_inbuf;
    	ancestor_2_inbuf.push(boost::iostreams::gzip_decompressor());
    	ancestor_2_inbuf.push(file_ancestor_2);
    	//Convert streambuf to istream
    	istream ancestor_2_instream(&ancestor_2_inbuf);
	
	cerr << endl << ancestor_2 << " opened." << endl;
	
	//(IN) ENDPOINT FILE 1 --------------------------------------------------
	ifstream file_endpoint_1(endpoint_1, ios_base::in | ios_base::binary);
    	
    	if (file_endpoint_1.is_open() == false) {
		cerr << endpoint_1 << " not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> endpoint_1_inbuf;
    	endpoint_1_inbuf.push(boost::iostreams::gzip_decompressor());
    	endpoint_1_inbuf.push(file_endpoint_1);
    	//Convert streambuf to istream
    	istream endpoint_1_instream(&endpoint_1_inbuf);
	
	cerr << endl << endpoint_1 << " opened." << endl;
	
	//(IN) ENDPOINT FILE 2 --------------------------------------------------
	ifstream file_endpoint_2(endpoint_2, ios_base::in | ios_base::binary);
    	
    	if (file_endpoint_2.is_open() == false) {
		cerr << endpoint_2 << " not opened: exit!" <<endl;
		return -1;
	}
	//uncompress
    	boost::iostreams::filtering_streambuf<boost::iostreams::input> endpoint_2_inbuf;
    	endpoint_2_inbuf.push(boost::iostreams::gzip_decompressor());
    	endpoint_2_inbuf.push(file_endpoint_2);
    	//Convert streambuf to istream
    	istream endpoint_2_instream(&endpoint_2_inbuf);
	
	cerr << endl << endpoint_2 << " opened." << endl;
	
	
	//*******************************************
	//          PARSING VARIABLES
	//*******************************************
	//Ancestor 1 parsing variables-------------
	string ancestor_1_line; //line of the file to parse
	string ancestor_1_chromosome; //pos
	int ancestor_1_chromosome_number = 0;
	long int ancestor_1_base_number = 0; //1-based
	string ancestor_1_ATCG_data;
	
	//Ancestor 2 parsing variables-------------
	string ancestor_2_line; //line of the file to parse
	string ancestor_2_chromosome; //pos
	int ancestor_2_chromosome_number = 0;
	long int ancestor_2_base_number = 0; //1-based
	string ancestor_2_ATCG_data;
	
	//Endpoint 1 parsing variables-------------
	string endpoint_1_line; //line of the file to parse
	string endpoint_1_chromosome; //pos
	int endpoint_1_chromosome_number = 0;
	long int endpoint_1_base_number = 0; //1-based
	string endpoint_1_ATCG_data;
	
	
	//Endpoint 2 parsing variables-------------
	string endpoint_2_line; //line of the file to parse
	string endpoint_2_chromosome; //pos
	int endpoint_2_chromosome_number = 0;
	long int endpoint_2_base_number = 0; //1-based
	string endpoint_2_ATCG_data;
	
	//*******************************************
	//          USEFUL VARIABLES
	//*******************************************
	long int total_lines = 0;
	long int both_not_muted = 0;
	long int both_muted = 0;
	long int both_same_mutation = 0;
	long int one_muted_and_one_not = 0;
	
	char endpoint_1_muted_in;
	char endpoint_2_muted_in;
	double endpoint_1_mut_freq;
	double endpoint_2_mut_freq;
	
	int extants [10] = {8, 12, 16, 20, 24, 28, 32, 48, 75, 100};
	
	long int endpoint_1_muted_bases [10] = {}; //10 long int values, each initialized with a value of zero
	long int endpoint_1_not_muted_bases [10] = {};
	double endpoint_1_P_hat [10] = {};
	
	long int endpoint_2_muted_bases [10] = {}; 
	long int endpoint_2_not_muted_bases [10] = {};
	double endpoint_2_P_hat [10] = {};
	
	
	//ROOT
	TApplication myApp("myApp",0,0);
	
	//Histograms (frequecy_1, frequency_2) 
	TH2D *h_all = new TH2D("All mutations", "All mutations", 201, 0., 1.005, 201, 0., 1.005);
	
	TH2D *h_same = new TH2D("Same mutation for both endpoints", "Same mutation for both endpoints", 201, 0., 1.005, 201, 0., 1.005);
	
	TH2D *h_different = new TH2D("Different mutation for the two endpoints", "Different mutation for the two endpoints", 201, 0., 1.005, 201, 0., 1.005);
	
	//*******************************************
	//          START
	//*******************************************
	cout<< endl << "Start reading files..." << endl;
	//Loop for all line 
	while (getline(ancestor_1_instream, ancestor_1_line)) {
		//Read a line from each file
		getline(ancestor_2_instream, ancestor_2_line);
		getline(endpoint_1_instream, endpoint_1_line);
		getline(endpoint_2_instream, endpoint_2_line);
	
		// Parse ancestor_1 line 
		stringstream ancestor_1_linestream(ancestor_1_line);
		//get file_2 line pos finding the separator ('\t')
		getline(ancestor_1_linestream, ancestor_1_chromosome, '\t');
		//get file_2 loc
		ancestor_1_linestream >> ancestor_1_base_number;
		//we need chromosome numeber and loc
		//extract chromosome number from chromosome string
		// "chrN", N integer
		//     ^        
		ancestor_1_chromosome_number = atoi(&ancestor_1_chromosome[3]);
		//get file_2 ATCG data finding the separator ('\n')
		getline(ancestor_1_linestream, ancestor_1_ATCG_data, '\n');
		
		// Parse ancestor_2 line
		stringstream ancestor_2_linestream(ancestor_2_line);
		getline(ancestor_2_linestream, ancestor_2_chromosome, '\t');
		ancestor_2_linestream >> ancestor_2_base_number;
		ancestor_2_chromosome_number = atoi(&ancestor_2_chromosome[3]);
		getline(ancestor_2_linestream, ancestor_2_ATCG_data, '\n');
		
		// Parse endpoint_1 line 
		stringstream endpoint_1_linestream(endpoint_1_line);
		getline(endpoint_1_linestream, endpoint_1_chromosome, '\t');
		endpoint_1_linestream >> endpoint_1_base_number;
		endpoint_1_chromosome_number = atoi(&endpoint_1_chromosome[3]);
		getline(endpoint_1_linestream, endpoint_1_ATCG_data, '\n');
		
		// Parse endpoint_2 line
		stringstream endpoint_2_linestream(endpoint_2_line);
		getline(endpoint_2_linestream, endpoint_2_chromosome, '\t');
		endpoint_2_linestream >> endpoint_2_base_number;
		endpoint_2_chromosome_number = atoi(&endpoint_2_chromosome[3]);
		getline(endpoint_2_linestream, endpoint_2_ATCG_data, '\n');
		
		//check that we loaded corrisponding lines (=> input files are ok)
		assert((ancestor_1_chromosome_number == ancestor_2_chromosome_number) 
			&& (ancestor_2_chromosome_number == endpoint_1_chromosome_number) 
			&& (ancestor_2_chromosome_number == endpoint_2_chromosome_number));
			
		assert((ancestor_1_base_number == ancestor_2_base_number)
			&& (ancestor_2_base_number == endpoint_1_base_number)
			&& (ancestor_2_base_number == endpoint_2_base_number));
		
		//Count loaded lines
		total_lines++;
		if (total_lines % 100000000 == 0){
			cout << total_lines/100000000 << "00 millions of lines read." << endl;
		}
		
		//Check ancestor coverage and discard according to threshold
		if (compute_coverage(ancestor_1_ATCG_data) < ancestor_coverage_min || 
			compute_coverage(ancestor_2_ATCG_data) < ancestor_coverage_min) {
			//one or both of the ancestors have a coverage under threshold
			//we discard the corresponding base 
			discarded_by_threshold++;
			//skip this base and don't put it in the histograms
			//load next lines from files files
			continue; 
				
		}
		
		//Check if endpoint 1 muted
		endpoint_1_muted_in = muted_in(endpoint_1_ATCG_data);
		if (endpoint_1_muted_in != ' '){
			endpoint_1_mut_freq = compute_mut_freq(endpoint_1_ATCG_data, endpoint_1_muted_in);
		} else endpoint_1_mut_freq = 0.;
		
		//Check if endpoint 2 muted
		endpoint_2_muted_in = muted_in(endpoint_2_ATCG_data);
		if (endpoint_2_muted_in != ' '){
			endpoint_2_mut_freq = compute_mut_freq(endpoint_2_ATCG_data, endpoint_2_muted_in);
		} else endpoint_2_mut_freq = 0.;
		
		//Check endpointand discard according to threshold
		if (endpoint_1_mut_freq > endpoint_max_mut_freq 
		    || endpoint_2_mut_freq > endpoint_max_mut_freq) {
		    	discarded_by_threshold++;
		    	continue;
		    }
		
		//compute enpoint 1 and 2 coverages 
		int endpoint_1_coverage = compute_coverage(endpoint_1_ATCG_data);
		int endpoint_2_coverage = compute_coverage(endpoint_2_ATCG_data);
		
		if (endpoint_1_muted_in == ' ' && endpoint_2_muted_in == ' '){
			//Both endpoints totally confim reference 
			both_not_muted++;
			//count these lines in both ancestors' denominator
			//for all extants 
			for (int i = 0; i < 10; i++){
				endpoint_1_not_muted_bases [i] ++;
				endpoint_2_not_muted_bases [i] ++;
			}
		} else {
			//at least one muted
			h_all -> Fill(endpoint_1_mut_freq, endpoint_2_mut_freq);
			
			if(endpoint_1_muted_in == ' ' || endpoint_2_muted_in == ' '){
				//one muted and one not
				one_muted_and_one_not++;
				h_different -> Fill(endpoint_1_mut_freq, endpoint_2_mut_freq);
				//Need to count that line based on the cutoff(extant)
				//for one of the two endpoints
				
				if (endpoint_2_muted_in == ' '){
					//endpoint 1 muted and 2 not
					
					//label endpoint_1 as muted or not
					//using cutoff varying extant
					for (int i = 0; i < 10; i++){
						//base is not muted for endpoint 2
						endpoint_2_not_muted_bases [i] ++;
						//check cutoff for endpoint 1
						if (endpoint_muted(endpoint_1_mut_freq, endpoint_1_coverage, extants[i]) ){
							endpoint_1_muted_bases [i] ++;
						} else{
							endpoint_1_not_muted_bases [i] ++;
						}
					}
				} else {
					//endpoint 2 muted and 1 not
					//label endpoint_2 as muted or not
					//using cutoff varying extant
					for (int i = 0; i < 10; i++){
						//base not muted for endpoint 1
						endpoint_1_not_muted_bases [i] ++;
						//check cutoff for endpoint 2
						if (endpoint_muted(endpoint_2_mut_freq, endpoint_2_coverage, extants[i]) ){
							endpoint_2_muted_bases [i] ++;
						} else{
							endpoint_2_not_muted_bases [i] ++;
						}
					}
				}
			} else {
				both_muted++;
				if(endpoint_1_muted_in == endpoint_2_muted_in) {
					//Same mutation for both endpoints
					both_same_mutation++;
					h_same -> Fill(endpoint_1_mut_freq, endpoint_2_mut_freq);
					//Skip these lines for the mut rate computation
				} else {
					//both muted but bifferent mutation
					h_different -> Fill(endpoint_1_mut_freq, endpoint_2_mut_freq);
					//Need to count these lines based on the cutoff(extant)
					//for both the endpoints
					
					for (int i = 0; i < 10; i++){
						//check cutoff for endpoint 1
						if (endpoint_muted(endpoint_1_mut_freq, endpoint_1_coverage, extants[i]) ){
							endpoint_1_muted_bases [i] ++;
						} else{
							endpoint_1_not_muted_bases [i] ++;
						}
						//check cutoff for endpoint 2
						if (endpoint_muted(endpoint_2_mut_freq, endpoint_2_coverage, extants[i]) ){
							endpoint_2_muted_bases [i] ++;
						} else{
							endpoint_2_not_muted_bases [i] ++;
						}
					}					
					
				}
			}
		}
		
	} //EOF
	
	assert(total_lines - discarded_by_threshold == both_not_muted + both_muted + one_muted_and_one_not);
	assert(both_same_mutation <= both_muted);
	
	cout <<"Total bases in files: " << total_lines << endl<<endl;
	cout <<"Total lines discarded by threshold (ancestor coverage min): " << discarded_by_threshold << endl;
	cout <<"Both endpoints not muted: " << both_not_muted <<endl;
	cout <<"Both endpoints muted: " << both_muted << endl;
	cout <<"Both endpoints with the same mutation: " << both_same_mutation << endl;
	cout <<"One endpoint muted and one not: "<< one_muted_and_one_not <<endl;
	
	//*******************************************
	//          FILES CLEANUP
	//*******************************************
	file_ancestor_1.close();
	file_ancestor_2.close();
	file_endpoint_1.close();
	file_endpoint_2.close();
	
	//*******************************************
	//         COMPUTE MUT RATE
	//*******************************************
	cout << "Total lines considered for the estimation of mut rate (denominator): "<< total_lines - discarded_by_threshold - both_same_mutation <<endl;
	double min_mu, max_mu;
	
	
	ofstream endpoint_1_results;
	ofstream endpoint_2_results;
	endpoint_1_results.open("CRC1307-02-1-E.all_common_results.csv");
	endpoint_2_results.open("CRC1307-09-1-B.all_common_results.csv");
	
	cout << endl << endl << "Computing mutation rate for endpoint CRC1307-02-1-E"<< endl<< endl;
	cout << "extant\tmuted_bases\tp_hat\tmin_mu\tmax_mu" << endl;
	
	endpoint_1_results << "extant,p_hat,min_mu,max_mu" << endl;
	
	for (int i = 0; i < 10; i++){
		assert(total_lines - discarded_by_threshold - both_same_mutation == endpoint_1_muted_bases[i] + endpoint_1_not_muted_bases[i]);
		endpoint_1_results << extants [i] << ",";
		cout << extants [i] << "\t";
		cout << endpoint_1_muted_bases[i] << "\t";
		
		endpoint_1_P_hat [i] = double(endpoint_1_muted_bases[i]) / double(endpoint_1_muted_bases[i] + endpoint_1_not_muted_bases[i]);
		
		endpoint_1_results << endpoint_1_P_hat [i] << ",";
		cout << endpoint_1_P_hat [i] << "\t";
		
		min_mu = compute_mut_rate(endpoint_1_P_hat [i], 0.45, extants [i]);
		max_mu = compute_mut_rate(endpoint_1_P_hat [i], 0.1, extants [i]);
		
		endpoint_1_results << min_mu << "," << max_mu << endl;
		cout << min_mu << "\t" << max_mu << endl;
	}
	
	cout << endl << endl << "Computing mutation rate for endpoint CRC1307-09-1-B"<< endl <<endl;
	cout << "extant\tmuted_bases\tp_hat\tmin_mu\tmax_mu" << endl;
	
	endpoint_2_results << "extant,p_hat,min_mu,max_mu" << endl;
	
	for (int i = 0; i < 10; i++){
		assert(total_lines - discarded_by_threshold - both_same_mutation == endpoint_2_muted_bases[i] + endpoint_2_not_muted_bases[i]);
		endpoint_2_results << extants [i] << ",";
		cout << extants [i] << "\t";
		cout << endpoint_2_muted_bases[i] << "\t";
		
		endpoint_2_P_hat [i] = double(endpoint_2_muted_bases[i]) / double(endpoint_2_muted_bases[i] + endpoint_2_not_muted_bases[i]);
		
		endpoint_2_results << endpoint_2_P_hat [i] << ",";
		cout << endpoint_2_P_hat [i] << "\t";
		
		min_mu = compute_mut_rate(endpoint_2_P_hat [i], 0.45, extants [i]);
		max_mu = compute_mut_rate(endpoint_2_P_hat [i], 0.1, extants [i]);
		
		endpoint_2_results << min_mu << "," << max_mu << endl;
		cout << min_mu << "\t" << max_mu << endl;
	}
	
	
	endpoint_1_results.close();
	endpoint_2_results.close();
	
	cout << "Results written on .csv files." << endl << "Starting ROOT graphs." << endl;
	
	//Start Root graphics
	TCanvas c1("Same mutation", "Same mutation",800,800);
   	h_same -> GetXaxis() -> SetTitle("mut freq endpoint 1");
   	h_same -> GetYaxis() -> SetTitle("mut freq endpoint 2");
   	h_same -> GetZaxis() -> SetTitle("Counts");
   	c1.SetLogx();
   	c1.SetLogy();
   	c1.SetLogz();
   	h_same -> Draw("colz");
   	
   	TCanvas c2("Different mutation", "Different mutation",800,800);
   	h_different -> GetXaxis() -> SetTitle("mut freq endpoint 1");
   	h_different -> GetYaxis() -> SetTitle("mut freq endpoint 2");
   	h_different -> GetZaxis() -> SetTitle("Counts");
   	c2.SetLogx();
   	c2.SetLogy();
   	c2.SetLogz();
   	h_different -> Draw("colz");
   	
   	TCanvas c3("All mutations", "All mutations",800,800);
   	h_all -> GetXaxis() -> SetTitle("mut freq endpoint 1");
   	h_all -> GetYaxis() -> SetTitle("mut freq endpoint 2");
   	h_all -> GetZaxis() -> SetTitle("Counts");
   	c3.SetLogx();
   	c3.SetLogy();
   	c3.SetLogz();
   	h_all -> Draw("colz");
   	
   	myApp.Run(1);
	
	return 0;
}


char muted_in(string ATCG_data){

	char reference;
	int A_counter = 0;
	int a_counter = 0;
	int T_counter = 0;
	int t_counter = 0;
	int C_counter = 0;
	int c_counter = 0;
	int G_counter = 0;
	int g_counter = 0;
	int coverage = 0;
	char muted_in;
	
	stringstream ATCG_linestream(ATCG_data);
	ATCG_linestream >> reference;
	//counters
	ATCG_linestream >> A_counter; coverage += A_counter;
	ATCG_linestream >> T_counter; coverage += T_counter;
	ATCG_linestream >> C_counter; coverage += C_counter;
	ATCG_linestream >> G_counter; coverage += G_counter;
	ATCG_linestream >> a_counter; coverage += a_counter;
	ATCG_linestream >> t_counter; coverage += t_counter;
	ATCG_linestream >> c_counter; coverage += c_counter;
	ATCG_linestream >> g_counter; coverage += g_counter;
	
	//total counts (forward+reverse)
	int A = A_counter + a_counter;
	int T = T_counter + t_counter;
	int C = C_counter + c_counter;
	int G = G_counter + g_counter;
	
	if (reference == 'A'){
		int muted_reads = max(T,C);
		muted_reads = max(muted_reads,G);
		if (muted_reads != 0) {
			if (muted_reads == T) muted_in = 'T';
			if (muted_reads == C) muted_in = 'C';
			if (muted_reads == G) muted_in = 'G';
		} else muted_in = ' ';
	}
	if (reference == 'T'){
		int muted_reads = max(A,C);
		muted_reads = max(muted_reads,G);
		if (muted_reads != 0) {
			if (muted_reads == A) muted_in = 'A';
			if (muted_reads == C) muted_in = 'C';
			if (muted_reads == G) muted_in = 'G';
		} else muted_in = ' ';
	}
	if (reference == 'C'){
		int muted_reads = max(A,T);
		muted_reads = max(muted_reads,G);
		if (muted_reads != 0) {
			if (muted_reads == A) muted_in = 'A';
			if (muted_reads == T) muted_in = 'T';
			if (muted_reads == G) muted_in = 'G';
		} else muted_in = ' ';
	}
	if (reference == 'G'){
		int muted_reads = max(A,T);
		muted_reads = max(muted_reads,C);
		if (muted_reads != 0) {
			if (muted_reads == A) muted_in = 'A';
			if (muted_reads == T) muted_in = 'T';
			if (muted_reads == C) muted_in = 'C';
		} else muted_in = ' ';
	}
	
	return muted_in;
}
double compute_mut_freq(string ATCG_data, char muted_in){
	assert(muted_in == 'A' || muted_in == 'T' || muted_in == 'C' || muted_in == 'G');
	char reference;
	int A_counter = 0;
	int a_counter = 0;
	int T_counter = 0;
	int t_counter = 0;
	int C_counter = 0;
	int c_counter = 0;
	int G_counter = 0;
	int g_counter = 0;
	int coverage = 0;
	
	double mut_freq;
	int muted_reads=0;
	
	stringstream ATCG_linestream(ATCG_data);
	ATCG_linestream >> reference;
	
	assert(muted_in != reference);
	
	//counters
	ATCG_linestream >> A_counter; coverage += A_counter;
	ATCG_linestream >> T_counter; coverage += T_counter;
	ATCG_linestream >> C_counter; coverage += C_counter;
	ATCG_linestream >> G_counter; coverage += G_counter;
	ATCG_linestream >> a_counter; coverage += a_counter;
	ATCG_linestream >> t_counter; coverage += t_counter;
	ATCG_linestream >> c_counter; coverage += c_counter;
	ATCG_linestream >> g_counter; coverage += g_counter;
	
	if(muted_in == 'A'){
		muted_reads = A_counter + a_counter;
		mut_freq = double(muted_reads)/double(coverage);
		return mut_freq;
	}
	if(muted_in == 'T'){
		muted_reads = T_counter + t_counter;
		mut_freq = double(muted_reads)/double(coverage);
		return mut_freq;
	}
	if(muted_in == 'C'){
		muted_reads = C_counter + c_counter;
		mut_freq = double(muted_reads)/double(coverage);
		return mut_freq;
	}
	if(muted_in == 'G'){
		muted_reads = G_counter + g_counter;
		mut_freq = double(muted_reads)/double(coverage);
		return mut_freq;
	}
	
}

int compute_coverage(string ATCG_data){
	char reference;
	int A_counter = 0;
	int a_counter = 0;
	int T_counter = 0;
	int t_counter = 0;
	int C_counter = 0;
	int c_counter = 0;
	int G_counter = 0;
	int g_counter = 0;
	int coverage = 0;

	stringstream ATCG_linestream(ATCG_data);
	ATCG_linestream >> reference;
	
	//counters
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

bool endpoint_muted(double mut_freq, int total_coverage, int extant){
	/*Count base as muted if endpoint frequency 

	f >= Cutoff(Coverage, f_min) e Coverage = TotCoverage/3

	cutoff (Coverage, f_min) = f_min + alpha/Sqrt[Coverage]

based on extant f_min and corresponding alphas are

	f_min = {1/8, 1/12, 1/16, 1/20, 1/24, 1/28, 1/32, 1/48, 1/75, 1/100}

	alpha = {0.52, 0.45, 0.36, 0.32, 0.3, 0.28, 0.25, 0.2, 0.18, 0.15} 
*/
	assert( extant == 8 || extant == 12 || extant == 16 || extant == 20
		|| extant == 24 || extant == 28 || extant == 32 || extant == 48
		|| extant == 75 || extant == 100);
		 
	double f_min = 0.;
	float coverage = float(total_coverage)/3.;
	double alpha = 0.;
	
	f_min = 1./double(extant);
	
	if (extant == 8){
		alpha = 0.52;
	}
	if (extant == 12){
		alpha = 0.45;
	}
	
	if (extant == 16){
		alpha = 0.36;
	}
	if (extant == 20){
		alpha = 0.32;
	}
	if (extant == 24){
		alpha = 0.3;
	}
	if (extant == 28){
		alpha = 0.28;
	}
	if (extant == 32){
		alpha = 0.25;
	}
	if (extant == 48){
		alpha = 0.2;
	}
	
	if (extant == 75){
		alpha = 0.18;
	}
	if (extant == 100){
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
