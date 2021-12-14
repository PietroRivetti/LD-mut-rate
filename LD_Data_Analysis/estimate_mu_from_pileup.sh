#!/bin/bash

ancestor_pileup="Data/ancestor.pileup.gz"
endpoint_pileup="Data/endpoint.pileup.gz"

ancestor_ATCG="Output_Files/ancestor.processed_pileup.gz"
endpoint_ATCG="Output_Files/endpoint.processed_pileup.gz"

ancestor_ATCG_common="Output_Files/ancestor.common_processed_pileup.gz"
endpoint_ATCG_common="Output_Files/endpoint.common_processed_pileup.gz"


echo ""
echo "Compile codes: "
cd Code
make 
cd ..

echo ""

echo "Start converting pileup files to ATCG human readable data."
echo ""
echo "Start converting ancestor..."
echo ""

zcat $ancestor_pileup | ./Code/process_pseudo_mpileup.x | gzip > $ancestor_ATCG

echo ""
echo "Ancestor's pileup conversion finished. Start converting endpoint..."
echo ""

zcat $endpoint_pileup | ./Code/process_pseudo_mpileup.x | gzip > $endpoint_ATCG

echo ""
echo "Endpoint's pileup conversion finished."
echo ""
echo "Start reading ancestor and endpoint ATCG data looking for common lines..."

./Code/common_lines_to_files.x $ancestor_ATCG $endpoint_ATCG $ancestor_ATCG_common $endpoint_ATCG_common

echo ""
echo "Start program to estimate mutation rate."
echo ""

./Code/estimate_mut_rate.x 20 $ancestor_ATCG_common $endpoint_ATCG_common

echo ""
echo "Terminated."
echo ""

