/*

This program reads from stdin uncopressed "pseudo_mpileup" (given by collaborators), 
parses every line from pileup format and outputs to stdout formatted as explicit ATCG counters.

COMMAND LINE USAGE:

	- with uncompressed pileup (and uncompressed output file):
		cat pseudo_mpileup_file | ./process_pseudo_mpileup.x > output_file.txt
		
	- with compressed pileup (and output file)
		zcat pseudo_mpileup_file.gz | ./process_pseudo_mpileup.x | gzip > output_file.txt.gz
		
FORMATTED INPUT FILE EXAMPLE:

chr	loc(0-based)	loc(1-based)	ref	#reads	read_chars
chr1	888658	888659	T	14	c$CCCCcccCcCccc	

OUTPUT EXAMPLE:

chr	loc	ref	A	T	C	G	a	t	c	g	
chr1	888659	T	0	9	0	2	0	0	0	0


NOTE THAT:

- input file has 0-based and 1-based loc -> here is used 1-based loc;

- insertions and deletions are omitted, as well as mapping quality.

- In the original format 'N' and 'n' exists: here we just check that they don't occur.

- In the original file the symbol '*' sometimes is counted as a read and sometimes it isn't. 
  So we can't check if the number of reads in the input file matches the sum of our counters 
  (see below the comment for '*' and assert at the end of the program).
  
  "*" is a placeholder for a deleted base in a multiple basepair deletion that was mentioned 
  in a previous line by the -[0-9]+[ACGTNacgtn]+ notation.		
  "Sometimes" they are counted as a read in the original file: when in the line they are mixed 
  with "normal" reads. So evey '*' -> read = read - 1
  but "some other times" they are not (ex. lines with 0 reads and only '*')


*/

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

int main(){
		
	string line;
	
	string chromosome; //loc
	long int base_number; //1-based
	char reference; //ref
	int reads;
	string char_results;
	
	int in_del_length; //length of insertion or deletion in bp
	
	long int line_counter = 0;
	
	
	while (getline(cin, line)){ //until EOF
		
		//parse a line 
	
		stringstream linestream(line);
		
		//get chr finding the separator ('\t')
		getline(linestream, chromosome, '\t');
		
		//get loc
		linestream >> base_number; //discard 0-based 
		linestream >> base_number;
		
		//get ref
		linestream >> reference;
		
		//get reads number
		linestream >> reads;
		
		//discard lines with only '*' and zero reads
		if (reads == 0) {line_counter++; continue; }  
		
		//get sequence of reads' chars
		getline(linestream, char_results);
		int results_length = char_results.length();
		
		//initialize counters 
		int dot_counter = 0; //count matching ref
		int comma_counter = 0;
		
		int A_counter = 0; //count base for forward strand
		int T_counter = 0;
		int C_counter = 0;
		int G_counter = 0;
		int N_counter = 0;
		
		int a_counter = 0; //count base for reverse strand
		int t_counter = 0;
		int c_counter = 0;
		int g_counter = 0;
		int n_counter = 0;	
			
		//parse char_results		
		for(int i = 0; i < results_length; i++){
			
			//A symbol `$' marks the end of a read segment. -> Skip
			if (char_results[i] == '$') continue;
			
			//"*" is a placeholder for a deleted base 
			//in a multiple basepair deletion that was mentioned 
			//in a previous line by the -[0-9]+[ACGTNacgtn]+ notation
			// -> Skip
			
			// "sometimes" they are counted as a read in the original file
			// (when in the line they are mixed with "normal" reads)
			// so evey '*' -> read = read - 1
			//
			// but "some other times" they are not 
			// (ex. lines with 0 reads and only '*')
			if (char_results[i] == '*') {reads--; continue; }

			
			//`^' marks the start of a read segment 
			//The ASCII of the character following `^' 
			//minus 33 gives the mapping quality. 
			// -> Skip both
			if (char_results[i] == '^') {
				i++; //skip this char and the next one
				continue;
			}
			
			// Look for insertion or deletions:
			// \+[0-9]+[ACGTNacgtn]+' indicates there is an insertion 
			//  Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion
			if(char_results[i] == '+' || char_results[i] == '-'){
				
				// check in/del length (need to convert char to int)
				in_del_length = atoi(&char_results[i+1]);
				
				// count how many digits to represent the in/del lenght
				int digits = 0; 
				int number = in_del_length;
				while (number != 0) { number /= 10; digits++; }
				
				// skip symbol, length's digits and chars of the in/del
				i += (digits + in_del_length);
				continue;
			}
			
			//check reads and count them
			if (char_results[i] == '.') {dot_counter++; continue; }
			
			if (char_results[i] == ',') {comma_counter++; continue; }
			
			if (char_results[i] == 'A') {A_counter++; continue; }
			
			if (char_results[i] == 'T') {T_counter++; continue; }	
			
			if (char_results[i] == 'C') {C_counter++; continue; }	
			
			if (char_results[i] == 'G') {G_counter++; continue; }
			
			if (char_results[i] == 'N') {N_counter++; continue; }
			
			if (char_results[i] == 'a') {a_counter++; continue; }
			
			if (char_results[i] == 't') {t_counter++; continue; }	
			
			if (char_results[i] == 'c') {c_counter++; continue; }	
			
			if (char_results[i] == 'g') {g_counter++; continue; }
			
			if (char_results[i] == 'n') {n_counter++; continue; }
			
		} //end to parse char_results
		
		//set dot and commas counted according to ref
		if (reference == 'A') {A_counter = dot_counter; a_counter = comma_counter; }
		
		if (reference == 'T') {T_counter = dot_counter; t_counter = comma_counter; }
			
		if (reference == 'C') {C_counter = dot_counter; c_counter = comma_counter; }
			
		if (reference == 'G') {G_counter = dot_counter; g_counter = comma_counter; }
		
		//output processed line				
		cout << chromosome << '\t' << base_number << '\t'<< reference 
		<< '\t' << A_counter << '\t' <<  T_counter << '\t' << C_counter << '\t' << G_counter 
		<< '\t' << a_counter << '\t' << t_counter << '\t' << c_counter << '\t' << g_counter
		<< endl;
		
		//check for errors:
		//Apparently fails only with lines that have zero reads and only '*'
		assert(A_counter + T_counter + C_counter + G_counter +
			a_counter + t_counter + c_counter + g_counter == reads);
		
		//check that no 'n' or 'N' were present
		assert (n_counter + N_counter == 0);	
    		
    		line_counter++;
    		
	} // each line until EOF
	
	cerr << "Pileup finished, total lines: " << line_counter << endl;
	
	return 0;
}
