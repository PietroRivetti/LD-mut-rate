LIBS:=-lz -lboost_iostreams

all: process_pseudo_mpileup common_lines_to_files estimate_mut_rate

process_pseudo_mpileup: 
	g++ process_pseudo_mpileup.cpp -o process_pseudo_mpileup.x
	
common_lines_to_files:
	g++ common_lines_to_files.cpp -o common_lines_to_files.x ${LIBS}
	
estimate_mut_rate:
	g++ estimate_mut_rate.cpp -o estimate_mut_rate.x ${LIBS}
