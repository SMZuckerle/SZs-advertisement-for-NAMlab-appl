#!/usr/bin/sh

###Define resources to be allocated - Higher resources will usually result in longer waiting time###

#CCS --res=rset=ncpus=8:mem=40g
#CCS --res=place=free:excl
#CCS -t 1h
#CCS --stderr=job_blastN_%reqid.err
#CCS --stdout=job_blastN_%reqid.out
#CCS --trace=job_blastN_%reqid.trace
#CCS -joe

#!/bin/bash

###Open the script and edit paths and or file names###
###Use chmod +x [job_Script.sh] to define Script###
###Copy fastas in blast folder:$PC2PFS/hpc-prf-ubo1/BLAST/ubo10001/ncbi-blast-2.9.0+/bin/ ###
###Use ccsalloc command to queue BLAST job ###

###INPUT library construction###
#prepare library with the following command if not done before: Change specifications in [...]#


makeblastdb -in $PC2PFS/hpc-prf-ubo1/ubo10001/BLAST/ncbi-blast-2.9.0+/bin/[file.fas] -dbtype nucl -parse_seqids

#further information:
#https://www.ncbi.nlm.nih.gov/books/NBK279684/

####INPUT Query#####
#prepare library with the following command if not done before: Change specifications in [...]#
Pfad=$PC2PFS/hpc-prf-ubo1/ubo10001/BLAST/ncbi-blast-2.9.0+/bin/	#path/to/query location/
Database=$PC2PFS/hpc-prf-ubo1/ubo10001/BLAST/ncbi-blast-2.9.0+/bin/[file.fas] #path/to/library ".fas/.fa"#
Query=[file]	#alle .fas
WS=5 #wordsize
eV=0.000000000000001	#eValue  e^-10=0.000045
cores=8	#threads
max_out=100	#num_alignments


#############################DO not alter from here on##############################################

#####blast#####
echo Start...
blastn -query $Pfad$Query.fas -db $Database -out results_$Query.txt -outfmt "6 std" -word_size $WS -evalue $eV -num_threads $cores -num_alignments $max_out	#table output

blastn -query $Pfad$Query.fas -db $Database -out results_$Query.html -html -word_size $WS -evalue $eV -num_threads $cores -num_alignments $max_out	#html output

blastn -query $Pfad$Query.fas -db $Database -outfmt "6 sallacc" -out resultstab.txt -word_size $WS -evalue $eV -num_threads $cores -num_alignments $max_out	#output for extraction

uniq resultstab.txt >resultstab2.txt

echo Extract...
blastdbcmd -db $Database -dbtype nucl -entry_batch resultstab2.txt -outfmt "%f" -out hitcontigs.fas -line_length 1000000	# extract hits from database


#####cleaning up#####
echo cleaning up...

mkdir blastn_"$Query"_WS-"$WS"_eV-"$eV"
cp $Pfad$Query.fas blastn_"$Query"_WS-"$WS"_eV-"$eV"
mv hitcontigs.fas blastn_"$Query"_WS-"$WS"_eV-"$eV"
mv results_$Query.txt blastn_"$Query"_WS-"$WS"_eV-"$eV"
mv results_$Query.html blastn_"$Query"_WS-"$WS"_eV-"$eV"
mv blastn_"$Query"_WS-"$WS"_eV-"$eV" "$Pfad"


rm resultstab.txt
rm resultstab2.txt

echo "done !"

exit $?
