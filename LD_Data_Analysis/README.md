# Estimate mutation rate from sequencing data using Luria-Delbrück model

This folder contains the code used to analyse the dataset from the ancestor and endpoint sequecing but the real dataset is not present.

The script `estimate_mu_from_pileup.sh` shows how the code has been used. An example of the terminal output of this script used on example data (just 10k rows for each pileup) is given in `output_example.txt`. 

Here Boost and ZLib libraries have been used to read and write from compressed `.gz` text files (see https://www.boost.org/ and https://zlib.net/). 
