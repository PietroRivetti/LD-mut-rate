LIBS:=`root-config --libs` -lz -lboost_iostreams
INCS:=`root-config --cflags` 

all: make_common_ancestor common_lines_to_files estimate_mu

common_lines_to_files:
	clang++ -O3 common_lines_to_files.cpp -o common_lines_to_files.x ${LIBS}
	
make_common_ancestor:
	clang++ -O3 make_common_ancestor.cpp -o make_common_ancestor.x ${LIBS}
	
estimate_mu:
	clang++ -O3 estimate_mu.cpp -o estimate_mu.x ${INCS} ${LIBS}
