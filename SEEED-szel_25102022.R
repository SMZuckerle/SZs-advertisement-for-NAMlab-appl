### getting started ###

install.packages("dplyr")
install.packages("data.table")
install.packages("tidyverse")

library(dplyr)
library(data.table)
library(tidyverse)
library(BiocManager)
library(seqinr)

setwd("C:/Users/elena/OneDrive/Desktop/Paderborn mirror/New strategie june 2021/Human/AS-change-info/R")

###prior:blast transcripts of interest

Filename="results_3-transcripts.txt" #Blast-Output
Extractor="to-extract.txt" #File with transcript accession and position of snp, comma separated (e.g. NC_12345,988)
Filename2="to-extract-seq.txt" #same as Filename plus nucleotide positions -2 to +2, comma separated (e.g. NC_12345,988,A,C,C,T,G)
Filename3="3transcripts-subset.gff" #gff-subfile of transcripts of interest

fastafile <- read.fasta(file = "GRCh38_latest_protein_singleline.faa", seqtype = "AA",as.string = TRUE, set.attributes = FALSE) #protein fasta file in singleline format

fastafile_transcripts <- read.fasta(file="hitcontigs_Trinity-GG.fasta_11022020_sorted_uniq_12022021_singleline.fas", seqtype="DNA",as.string=TRUE,set.attributes=FALSE) #transcript fasta file in singleline format

###reading in subset gff file & extract information of interest (correlating transcript & protein accession)

a<-
  read.csv(
    Filename3, header = FALSE, stringsAsFactors = FALSE, sep = "\t"
  )

a1 <- subset(a, a$V3=="CDS")
a2 <- a1 %>% select(9)
a3 <- a2[,1]
class(a3)

a4 <- data.frame(do.call(rbind, strsplit(a3, c("-"), fixed=TRUE)))
a5 <- a4 %>% select(2,3)

a6<-cSplit(a5,"X2",sep=";")
a7<-cSplit(a6,"X3",sep=";")
a8<-a7 %>% select(1,3)
colnames(a8)[1]<-"protein"
colnames(a8)[2]<-"transcript"

a9 <- unique(a8)

###importing and merging files to combine SNP-Pos-Info with BLAST output info

x<-
  read.csv(
    Filename, header = FALSE, stringsAsFactors = FALSE, sep = "\t"
  )

x_ex<-  read.csv(
  Extractor, header = FALSE, stringsAsFactors = FALSE, sep = ","
)

x1<- x %>% select(1,2,7,8,9,10)

colnames(x1)[1]<-"transcript"
colnames(x_ex)[1]<-"transcript"
colnames(x_ex)[2]<-"Edit_nt"

xx<- merge.data.frame(x1,x_ex, by= "transcript", all.x = TRUE, all.y = TRUE)

###define blast hits, that cover the snp position

xx$edit1 <- ifelse(xx$Edit_nt >= xx$V7, 1,0)
xx$edit2 <- ifelse(xx$Edit_nt <= xx$V8, 1,0)
xx$edit3 <- ifelse(xx$edit1 + xx$edit2 == 2, 1,0)

xx2 <- subset(xx, xx$edit3==1)

###keep only protein accession, which is related to transcript of interest (GFF-information)

xx3 <- merge.data.frame(xx2,a9, by= "transcript", all.x = TRUE, all.y = TRUE)
s_xx <- subset(xx3, xx3$V2==xx3$protein)

###check for alignment length of blast hits, calculate amino acid position

s_xx$len <- (s_xx$V8-s_xx$V7)
s_xx$AA_pos <-((s_xx$Edit_nt-s_xx$V7)/3+s_xx$V9)
s_xx$AA_pos2 <- floor(s_xx$AA_pos)
colnames(s_xx)[2] <- "protein"

###extract info for AA identity extraction

y <- s_xx %>%
  ungroup() %>%
  select(2,13)

###safe protein RefSeq Fasta File as table, extract sequences of interest

#uses "fastafile" (protein fasta)

subsetlist <- s_xx %>%
  ungroup() %>%
  select(2,14)

names(subsetlist) <- c("ID", "POS")
subset_fasta <- fastafile[names(fastafile) %in% subsetlist$ID]

ttt<- as.data.frame(subset_fasta, header= FALSE)
ttt1<-as.data.frame(t(ttt))
ttt2<-tibble::rownames_to_column(ttt1,"value")

colnames(ttt2)[1]<-"ID"
colnames(y)[1]<-"ID"

subset_fasta_table <- merge.data.frame(y,ttt2, by= "ID", all.x = TRUE, all.y = TRUE)

names(subset_fasta_table) <- c("ID", "POS", "SEQ")

###find out amino acid identity, merge with previous table

subset_fasta_table$AA <- substr(subset_fasta_table$SEQ, subset_fasta_table$POS, subset_fasta_table$POS)
colnames(s_xx)[2] <- "ID"
z <- merge.data.frame(s_xx,subset_fasta_table, by= "ID", all.x = TRUE, all.y = TRUE)

z1 <- z %>%
  select(1,2,13,14,17)

###extracting target sequence -31 to +6

subset_fasta_transcripts <- fastafile_transcripts[names(fastafile_transcripts) %in% x_ex$transcript]

uuu<- as.data.frame(subset_fasta_transcripts, header= FALSE)
uuu1<-as.data.frame(t(uuu))
uuu2<-tibble::rownames_to_column(uuu1,"value")

colnames(uuu2)[1]<-"transcript"
colnames(uuu2)[2]<-"transcript_seq"

subset_fasta_transcript_table <- merge.data.frame(x_ex,uuu2, by= "transcript", all.x = TRUE, all.y = TRUE)

subset_fasta_transcript_table$start <- (subset_fasta_transcript_table$Edit_nt-31)
subset_fasta_transcript_table$end <- (subset_fasta_transcript_table$Edit_nt+6)

subset_fasta_transcript_table$target_seq <- substr(subset_fasta_transcript_table$transcript_seq, subset_fasta_transcript_table$start, subset_fasta_transcript_table$end)

###extracting target sequence -2 to +2

subset_fasta_transcript_table$start_5 <- (subset_fasta_transcript_table$Edit_nt-2)
subset_fasta_transcript_table$end_5 <- (subset_fasta_transcript_table$Edit_nt+2)

subset_fasta_transcript_table$target_seq_sh <- substr(subset_fasta_transcript_table$transcript_seq, subset_fasta_transcript_table$start_5, subset_fasta_transcript_table$end_5)


###find out codon position of snp in the respective protein sequence

z1$AA_pos3 <- floor(z1$AA_pos)
z1$AA_pos4 <- (z1$AA_pos-z1$AA_pos3)
z1$AA_pos4a <- round(z1$AA_pos4, digits=1)
z1$AA_pos5 <- ifelse(z1$AA_pos4a == 0, 1,0)
z1$AA_pos6 <- ifelse(z1$AA_pos4a == 0.3, 2,0)
z1$AA_pos7 <- ifelse(z1$AA_pos4a == 0.7, 3,0)
z1$AA_pos8 <- (z1$AA_pos5+z1$AA_pos6+z1$AA_pos7)

z2 <- z1 %>%
  select(1,2,5,12)

colnames(z2)[3] <- "AA_DNA"
colnames(z2)[4] <- "codonposn."

###find out which AA codon is established by edit

target_5 <- subset_fasta_transcript_table %>%
  select(1,2,9)

target_5_1 <- separate(target_5, target_seq_sh, c("V4", "V5", "V6", "V7", "V8"), sep=c(1:4))
colnames(target_5_1) <- c("transcript", "V2", "V3", "V4", "V5", "V6", "V7")

z3 <- merge.data.frame(z2,target_5_1, by= "transcript", all.x = TRUE, all.y = TRUE)

z3$V5 <- gsub("c", "t", z3$V5)

z3$codon1 <- ifelse(z3$codonposn. == 1, paste(z3$V5, z3$V6, z3$V7, sep=""), NA)
z3$codon2 <- ifelse(z3$codonposn. == 2, paste(z3$V4, z3$V5, z3$V6, sep=""), NA)
z3$codon3 <- ifelse(z3$codonposn. == 3, paste(z3$V3, z3$V4, z3$V5, sep=""), NA)
z3$codon <- paste(z3$codon1, z3$codon2, z3$codon3, sep="")
z3$codon <- gsub('NA','',z3$codon)
z3$codon <- toupper(z3$codon)

#NERDY DADDY FISH credits

z3$C1<-ifelse(z3$codon == "TTT", "P", NA)
z3$C2<-ifelse(z3$codon == "TTC", "P", NA)
z3$C3<-ifelse(z3$codon == "TTA", "L", NA)
z3$C4<-ifelse(z3$codon == "TTG", "L", NA)

z3$C5<-ifelse(z3$codon == "CTT", "L", NA)
z3$C6<-ifelse(z3$codon == "CTC", "L", NA)
z3$C7<-ifelse(z3$codon == "CTA", "L", NA)
z3$C8<-ifelse(z3$codon == "CTG", "L", NA)

z3$C9<-ifelse(z3$codon == "ATT", "I", NA)
z3$C10<-ifelse(z3$codon == "ATC", "I", NA)
z3$C11<-ifelse(z3$codon == "ATA", "I", NA)
z3$C12<-ifelse(z3$codon == "ATG", "M", NA)

z3$C13<-ifelse(z3$codon == "GTT", "V", NA)
z3$C14<-ifelse(z3$codon == "GTC", "V", NA)
z3$C15<-ifelse(z3$codon == "GTA", "V", NA)
z3$C16<-ifelse(z3$codon == "GTG", "V", NA)

z3$C17<-ifelse(z3$codon == "TCT", "S", NA)
z3$C18<-ifelse(z3$codon == "TCC", "S", NA)
z3$C19<-ifelse(z3$codon == "TCA", "S", NA)
z3$C20<-ifelse(z3$codon == "TCG", "S", NA)

z3$C21<-ifelse(z3$codon == "CCT", "P", NA)
z3$C22<-ifelse(z3$codon == "CCC", "P", NA)
z3$C23<-ifelse(z3$codon == "CCA", "P", NA)
z3$C24<-ifelse(z3$codon == "CCG", "P", NA)

z3$C25<-ifelse(z3$codon == "ACT", "Y", NA)
z3$C26<-ifelse(z3$codon == "ACC", "Y", NA)
z3$C27<-ifelse(z3$codon == "ACA",  "*", NA)
z3$C28<-ifelse(z3$codon == "ACG",  "*", NA)

z3$C29<-ifelse(z3$codon == "GCT", "H", NA)
z3$C30<-ifelse(z3$codon == "GCC", "H", NA)
z3$C31<-ifelse(z3$codon == "GCA", "Q", NA)
z3$C32<-ifelse(z3$codon == "GCG", "Q", NA)

z3$C33<-ifelse(z3$codon == "TAT", "N", NA)
z3$C34<-ifelse(z3$codon == "TAC", "N", NA)
z3$C35<-ifelse(z3$codon == "TAA", "K", NA)
z3$C36<-ifelse(z3$codon == "TAG", "K", NA)

z3$C37<-ifelse(z3$codon == "CAT", "D", NA)
z3$C38<-ifelse(z3$codon == "CAC", "D", NA)
z3$C39<-ifelse(z3$codon == "CAA", "E", NA)
z3$C40<-ifelse(z3$codon == "CAG", "E", NA)

z3$C41<-ifelse(z3$codon == "TGT", "C", NA)
z3$C42<-ifelse(z3$codon == "TGC", "C", NA)
z3$C43<-ifelse(z3$codon == "TGA", "*", NA)
z3$C44<-ifelse(z3$codon == "TGG", "W", NA)

z3$C45<-ifelse(z3$codon == "CGT", "R", NA)
z3$C46<-ifelse(z3$codon == "CGC", "R", NA)
z3$C47<-ifelse(z3$codon == "CGA", "R", NA)
z3$C48<-ifelse(z3$codon == "CGG", "R", NA)

z3$AA_edited <- paste(z3$C1,z3$C2,z3$C3,z3$C4,z3$C5,z3$C6,z3$C7,z3$C8,z3$C9,z3$C10,z3$C11,z3$C12,z3$C13,z3$C14,z3$C15,z3$C16,z3$C17,z3$C18,z3$C19,z3$C20,z3$C21,z3$C22,z3$C23,z3$C24,z3$C25,z3$C26,z3$C27,z3$C28,z3$C29,z3$C30,z3$C31,z3$C32,z3$C33,z3$C34,z3$C35,z3$C36,z3$C37,z3$C38,z3$C39,z3$C40,z3$C41,z3$C42,z3$C43,z3$C44,z3$C45,z3$C46,z3$C47,z3$C48, sep="")

z4 <- z3 %>%
  select(1:14,63)

z4$AA_edited <- gsub('NA','',z4$AA_edited)

z5 <- z4 %>%
  select(1:5,14,15)

colnames(z5)[5] <- "edit_pos_transcr"

z5$editname <- str_c(z5$transcript,"eU",z5$edit_pos_transcr,z5$AA_DNA,z5$AA_edited, sep="")

###extract geneID

a<-
  read.csv(
    Filename3, header = FALSE, stringsAsFactors = FALSE, sep = "\t"
  )

b1 <- subset(a, a$V3=="mRNA")
b2 <- b1 %>% select(9)
b3 <- b2[,1]

b4 <- data.frame(do.call(rbind, strsplit(b3, c("-"), fixed=TRUE)))

b5 <- b4 %>% select(2,3)

b6<-cSplit(b5,"X2",sep=";")
b7<-cSplit(b6,"X3",sep=";")
b8<-b7 %>% select(1,3)
colnames(b8)[1]<-"transcript"
colnames(b8)[2]<-"geneID"

b9 <- unique(b8)

b10 <- merge(z5, b8, by="transcript", all.x=TRUE)

b10$AA_DNA[is.na(b10$AA_DNA)] <- "?UTR"

b10$fullname <- paste(b10$geneID, "-", b10$transcript,"eU",b10$edit_pos_transcr,b10$AA_DNA,b10$AA_edited, sep="")

b11 <- b10 %>%
  select(1,10)

###create final information table

final1 <- subset_fasta_transcript_table %>%
  select(1,6)

colheaders <- as.character(c(-31:6))

final2 <- separate(final1, target_seq, colheaders, sep=c(1:37))

final3 <- merge.data.frame(b11, final2, by= "transcript", all.x = TRUE, all.y = TRUE)
