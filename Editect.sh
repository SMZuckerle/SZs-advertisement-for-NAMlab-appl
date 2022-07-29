#!/usr/bin/sh

#CCS --res=rset=ncpus=16:mem=160g
#CCS --res=place=free:excl
#CCS -t 3h
#CCS --stderr=jobEDITECT_%reqid.err
#CCS --stdout=jobEDITECT_%reqid.out
#CCS --trace=jobEDITECT_%reqid.trace
#CCS -joe



#!/bin/bash
#####Input#####

#reference#
ref=$PC2PFS/hpc-prf-ubo1/ubo10001/rawdata/Haplopteris_Vittaria_DNA_RNA/Haplopteris_cpDNA_woIR_20191128.fas	#path/to/reference.fas

name=H2017_cpDNA_woIR_20191128
#raw-reads#
RNA_fwd=$PC2PFS/hpc-prf-ubo1/ubo10001/rawdata/Haplopteris_Vittaria_DNA_RNA/Haplopteris-ensiformis-2_1_RNA.fq	#path/to/reads

RNA_rev=$PC2PFS/hpc-prf-ubo1/ubo10001/rawdata/Haplopteris_Vittaria_DNA_RNA/Haplopteris-ensiformis-2_2_RNA.fq	#path/to/reads

DNA_fwd=$PC2PFS/hpc-prf-ubo1/ubo10001/rawdata/Haplopteris_Vittaria_DNA_RNA/HAPLO_2017_dna_FCHHFMVALXX_L8_8521612005041_1.fq	#path/to/reads

DNA_rev=$PC2PFS/hpc-prf-ubo1/ubo10001/rawdata/Haplopteris_Vittaria_DNA_RNA/HAPLO_2017_dna_FCHHFMVALXX_L8_8521612005041_2.fq	#path/to/reads

###output files###

#1. JACUSA-out - potential edits
#2. RNA.concordant.uniq sam file - for visualisation with tablet
#3. DNA.concordant.uniq sam file - for visualisation with tablet
#4. mapping ref.fas
#5. indexed ref.fas.fai -for visualisation with tablet

###Citations###
#Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N.,. . . Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

#Picardi, E., & Pesole, G. (2013). REDItools: High-throughput RNA editing detection made easy. Bioinformatics (Oxford, England), 29(14), 1813–1814. https://doi.org/10.1093/bioinformatics/btt287

#Piechotta, M., Wyler, E., Ohler, U., Landthaler, M., & Dieterich, C. (2017). JACUSA: Site-specific identification of RNA editing events from replicate sequencing data. BMC bioinformatics, 18(1), 7. https://doi.org/10.1186/s12859-016-1432-8

#Wu, T. D., Reeder, J., Lawrence, M., Becker, G., & Brauer, M. J. (2016). GMAP and GSNAP for Genomic Sequence Alignment: Enhancements to Speed, Accuracy, and Functionality. Methods in molecular biology (Clifton, N.J.), 1418, 283–334. https://doi.org/10.1007/978-1-4939-3578-9_15

###Discription###
###gsnap###

#  --gunzip                      	 Uncompress gzipped input files

#  --bunzip2                    	 Uncompress bzip2-compressed input files

#  -N, --novelsplicing=INT         	 Look for novel splicing (0=no (default), 1=yes)

#  -E, --distant-splice-penalty=INT      Penalty for a distant splice (default 1).  A distant splice is one where
#                                        the intron length exceeds the value of -w, or --localsplicedist, or is an
#                                        inversion, scramble, or translocation between two different chromosomes
#                                        Counts against mismatches allowed

#  -Q, --quiet-if-excessive     	 If more than maximum number of paths are found,
#                                   	 then nothing is printed.

#  -O, --ordered                 	 Print output in same order as input (relevant
#                              	         only if there is more than one worker thread)

#  --nofails                     	 Exclude printing of failed alignments

#  -A, --format=STRING           	 Another format type, other than default.
#                                	 Currently implemented: sam, m8 (BLAST tabular format)

#  --split-output=STRING  	         Basename for multiple-file output, separately for nomapping,
#                                   	 halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult,
#                                    	 paired_uniq, paired_mult, concordant_uniq, and concordant_mult results

#  -t, --nthreads=INT           	 Number of worker threads

###JACUSA##
#usage: jacusa.jar [OPTIONS] BAM1_1[,BAM1_2,BAM1_3,...] BAM2_1[,BAM2_2,BAM2_3,...]
# -a,--pileup-filter <PILEUP-FILTER>             chain of PILEUP-FILTER to apply to pileups:
#                                                B | Filter distance to Read Start/End. Default: 6:0.5 (F:distance:min_ratio)
#                                                R | Rare event filter. Default: 1:2:0.1 (R:pool:reads:level)
#                                                S | Filter distance to Splice Site. Default: 6:0.5 (S:distance:min_ratio)
#                                                D | Filter distance to Read Start/End, Intron, and INDEL position. Default: 5:0.5 (D:distance:min_ratio)
#                                                H | Filter non-homozygous pileup/BAM in sample 1 or 2 (MUST be set to H:1 or H:2). Default: none
#                                                I | Filter distance to INDEL position. Default: 6:0.5 (I:distance:min_ratio)
#                                                Y | Filter wrong variant calls in the vicinity of homopolymers. Default: 7 (Y:length)
#                                                L | Min difference filter. Default: 0:0:0.1 (R:pool:reads:level)
#                                                M | Max allowed alleles per parallel pileup. Default: 2

# -p,--threads <THREADS>                         use # THREADS
#                                                default: 1

# -T,--threshold <THRESHOLD>                     Filter positions based on test-statistic THRESHOLD
#                                                default: DO NOT FILTER

# -f,--output-format <OUTPUT-FORMAT>             Choose output format:
#                                                <*> B: Default
#                                                < > V: VCF Output format. Option -P will be ignored (VCF is unstranded)

# -r,--result-file <RESULT-FILE>                 results are written to RESULT-FILE or STDOUT if empty

#############################DO not alter from here on##############################################

#####etc#####
echo $ref | sed -E -e 's-/-\n-g' | tail -1 | sed -E -e 's/\.fas//g' > ref

ref2=`cat ref`

date=`date +%Y-%m-%d`

#Load Modules#

module load compiler/GCC/7.3.0-2.30 tools/bzip2/1.0.6-GCCcore-7.3.0 devel/ncurses/6.1-GCCcore-7.3.0

#####Mapping#####
#indexing#
if [ -d $PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/INDEX/"$ref2" ]
then
  	echo "Index $ref2 exist."
else
	echo "Index $ref2 does not exist" >&2
	echo "Index will be build."
	cat $ref | sed -E -e 's/>/>chr_/g' >ref.fas
	$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/gmap-2019-03-15/bin/gmap_build -D $PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/INDEX -d $ref2 -C ref.fas
	rm ref.fas
fi

#mapping-RNA#
$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/gmap-2019-03-15/bin/gsnap -D $PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/INDEX -d $ref2 -E 1000 -N1 -Q -O --nofails -A sam --gunzip --split-output=rna_gsnap_$ref2 -t 8 $RNA_fwd $RNA_rev
#mapping-DNA#
$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/gmap-2019-03-15/bin/gsnap -D $PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/INDEX -d $ref2 -E 1000 -Q -O --nofails -A sam --gunzip --split-output=dna_gsnap_$ref2 -t 8 $DNA_fwd $DNA_rev



#####Converting#####
RNA=rna_gsnap_$ref2.concordant_uniq
DNA=dna_gsnap_$ref2.concordant_uniq

$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools view -@ 8 -S -b $RNA > $RNA.bam
$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools view -@ 8 -S -b $DNA > $DNA.bam

$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools sort -@ 8 -o $RNA-sorted $RNA.bam
$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools sort -@ 8 -o $DNA-sorted $DNA.bam

$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools index -b $RNA-sorted
$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools index -b $DNA-sorted

$PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/samtools-1.9/samtools faidx $ref

#####RNA editing calling#####
java -Xmx160g -jar $PC2PFS/hpc-prf-ubo1/ubo10001/Mapping_tools/JACUSA_v1.3.0.jar call-2 -a H:1,M,B,Y -T 2.3 -f V -p 1 -r GSNAP_Jaccusa_DNA-RNA_"$ref2"_"$name"_$date.out $DNA-sorted $RNA-sorted



####clean-up####
#keep files in folder#
mkdir results_"$ref2"_"$name"_$date
mv GSNAP_Jaccusa_DNA-RNA_"$ref2"_"$name"_$date.out results_"$ref2"_"$name"_$date
mv $RNA results_"$ref2"_"$name"_$date
mv $DNA results_"$ref2"_"$name"_$date
mv $ref.fai results_"$ref2"_"$name"_$date
cp $ref results_"$ref2"_"$name"_$date
#remove unused files# (>  "#" in front of rm to keep all files in Editing_detection folder)
#rm -rf ref
#rm *$ref2.nom*
#rm *$ref2.*_*
####show results####
#cp /home/knoop/Editing_detection/results_Picea_glauca_CDS_MT_cor_mature_needle_2018-09-20/$DNA /home/knoop/Editing_detection
