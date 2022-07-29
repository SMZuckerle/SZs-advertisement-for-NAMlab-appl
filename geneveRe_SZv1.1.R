#"C:\Program Files\R\R-4.1.2\Rscript.exe"
#install.packages("dplyr")
#install.packages("data.table")
#install.packages("tidyverse")
#install.packages("seqinr")
library(dplyr)
library(data.table)
library(tidyverse)
library(seqinr)
###################################################################
###R###############
###################################################################
###################################################################
###R###############
###################################################################
###################################################################
###
### Guinevere:
###
### Get Access to all sequences as fasta that have been annotated
### (and not been annotated) in your Genebank file of interest.
###
###
###
###################################################################
#setwd("C:/Users/cyan4/OneDrive/Desktop/R_Home/GeneBank_for_mtg2")
###################################################################
Accession="NC_027000"
Filename="NC_027000.gff3"
Fastafile="NC_027000.fasta"
############Download complete sequence entry from ncbi ############
###################################################################
###################################################################
###################################################################
a<-
  read.csv(
    Filename, header = FALSE, stringsAsFactors = FALSE, sep = "\t"
  )
b<- a %>% select(4,5,9)

b$id <- paste(b$V4,b$V5,sep = "_")
a$id <- paste(a$V4,a$V5,sep = "_")
b1<-b[,c(4,3,1,2)]
#b2<- b1 %>% separate(b2, )
################################
b1$prot1<-str_extract(b1$V9,"protein_id=.*")
b1$prot2<-substr(b1$prot1,1,22)
b1$protein_id <- data.frame(do.call('rbind', strsplit(as.character(b1$prot2),'=',fixed=TRUE)))
################################
#test<-"transcript_id=NP_000066.1"
b1$trans1<-str_extract(b1$V9,"transcript_id=.*")
b1$trans2<-substr(b1$trans1,1,25)
b1$transcript_id <- data.frame(do.call('rbind', strsplit(as.character(b1$trans2),'=',fixed=TRUE)))
################################
b1$info1<-str_extract(b1$V9,"gbkey=.*")
pattern1 <- paste(c(".*;tag",
                    ".*;transcript_id",
                    ".*;protein"))
b1$info2 <- data.frame(do.call('rbind', strsplit(as.character(b1$info1),';',fixed=TRUE)))
################################
b2 <- as.data.table(b1)
b3 <- b2 %>% select(3,4,8,12)
#colnames(b2)
bb1 <- select(b2, c("id","info2.X1","info2.X2","info2.X3"))
b3$id <- paste(b3$V4,b3$V5,sep = "_")
a2 <- a %>% select(1,2,3,10)
c<- merge.data.frame(
  a2,b3, by = "id",
  all.x = TRUE, all.y = TRUE)
c1<- merge.data.frame(
  c,bb1, by = "id",
  all.x = TRUE, all.y = TRUE)
c2 <- unique(c1)

#test <- c[,c('id','V1')] #SELECTION OF COLUMnS BY NAME b3, a2!!!!!!!!!!!!!
#test <- select(c, c("id","V1"))
#USE data.frame "c" to filter and select data based on features, for e.g. fasta extraction
write.table(c, file = paste(Filename,"_GFF-table_format.txt"), quote = FALSE, row.names = FALSE, sep = ";")
###############################################################
fas_a <- read.fasta(
  Fastafile, seqtype = "DNA",as.string = TRUE, set.attributes = FALSE) #protein fasta file in singleline format

ttt<- as.data.frame(fas_a, header= FALSE)
ttt1<-as.data.frame(t(ttt))
ttt2<-tibble::rownames_to_column(ttt1,"value")
colnames(ttt2)[1]<-"Acc"
colnames(c2)[2]<-"Acc"
fas_table <- merge.data.frame(
  c2,
  ttt2,
  by= "Acc", all.x = TRUE, all.y = TRUE)
colnames(fas_table
         )[12]<-"seq"
fas_table$seq1 <- substr(fas_table$seq, fas_table$V4, fas_table$V5
                         )
fas_table1 <- fas_table %>% select(1:11,13
                                   )
colnames(fas_table1
         )[4]<-"type"
fas_table1$id <- gsub("_","..",fas_table1$id
                      )
fas_table1$name <- paste(
  str_extract(fas_table1$info2.X1,"=.*" ),
  str_extract(fas_table1$info2.X2,"=.*"),
  fas_table1$Acc, fas_table1$id,
  str_extract(fas_table1$info2.X3,"=.*"), sep = "_"
  )
fas_table1$name <- gsub("=","",fas_table1$name
                        )
fas_table1$info2.X1 <- gsub("gbkey=","",fas_table1$info2.X1
                            )
colnames(fas_table1)[9]<-"sele"
colnames(fas_table1)[5]<-"sta"
colnames(fas_table1)[6]<-"end"
fas_table2 <- fas_table1 %>% select(1,2,5,6,9,13,12)
###################################################################ALL COOL HERE - NOW GET IGS (GENE) AND INTRONS (CDSs)
gene_table0 <- fas_table1 %>% select(1,2,5,6,9,10)
gene_table <- subset(gene_table0,
                     gene_table0$sele=="Gene")
gene_table <- unique(gene_table)
acc_pos <- gene_table %>% select(1,3,4,6)
colnames(acc_pos)[4]<- "type"
#insert start-end dummies
start_dummi <- unique(
  data.frame(
    acc_pos$Acc,1,1,"start")
)
dummies <- subset(c2,
                     c2$info2.X1=="gbkey=Src")

end_dummi <- unique(
  data.frame(
    acc_pos$Acc,dummies$V5,dummies$V5,"end")
)
names(start_dummi) <- c("Acc","sta","end","type")
names(end_dummi) <- c("Acc","sta","end","type")
acc_pos1 <- rbind(acc_pos,start_dummi,end_dummi)
#colnames(acc_pos)[2]<-"POS"
#colnames(acc_pos)[3]<-"POS2"
#acc_pos1<- sort(acc_pos1$sta, decreasing = FALSE, na.last = TRUE)
acc_pos1 <-acc_pos1[order(acc_pos1$sta),]
acc_pos1$tag <- c(
  1:nrow(acc_pos1)
  )
acc_pos1$tag2 <- c(
  1:nrow(acc_pos1)
)
acc_pos1$tag2 <- acc_pos1$tag+1

acc_pos2 <- acc_pos1 %>% select(1,3,4,6)
colnames(acc_pos2)[4]<-"tag"
#acc_pos2$tag <- c(
#  1:nrow(acc_pos)
#  )
acc_pos3 <- merge.data.frame(
  acc_pos1,acc_pos2, by= "tag", all.x = TRUE, all.y = TRUE)
igs_table <- acc_pos3 %>% select(2,8,3,5,9)
#igs_table0$pair <- ifelse(igs_table0$info2.X2.x==igs_table0$info2.X2.y,1,0)
#igs_table <- subset(igs_table0,igs_table0$pair==0)
colnames(igs_table) [1] <- "Acc"
colnames(igs_table) [2] <- "sta"
colnames(igs_table) [3] <- "end"
colnames(igs_table) [4] <- "gene1"
colnames(igs_table) [5] <- "gene2"
igs_table$gene1 <- gsub("gene=","",igs_table$gene1
)
igs_table$gene2 <- gsub("gene=","",igs_table$gene2
)
igs_table$id <- paste(igs_table$sta, igs_table$end, sep = "..")
igs_table$type <- paste(c("IGS"))
igs_table$name <- paste(
  igs_table$type,
  igs_table$gene1,
  igs_table$gene2,
  igs_table$Acc,
  igs_table$id,
  sep = "_")
igs_table$ori <- ifelse(igs_table$sta<=igs_table$end,1,0)
igs_table02 <- subset(igs_table,
                      igs_table$ori==1)
igs_table03 <- subset(igs_table,
                      igs_table$ori==0)
colnames(igs_table03) [2] <- "end"
colnames(igs_table03) [3] <- "sta"

igs_table01<- rbind(igs_table02,igs_table03)

igs_table_DONT_VIEW<- merge.data.frame(
  igs_table01,ttt2, by= "Acc", all.x = TRUE, all.y = TRUE)
igs_table_DONT_VIEW$seq <- substr(
  igs_table_DONT_VIEW$V1,
  igs_table_DONT_VIEW$sta,
  igs_table_DONT_VIEW$end
)
igs_table1 <- igs_table_DONT_VIEW %>% select(1:9,11
)
###################################################################
####SEQ EXTRACTION################################################

seq00 <- igs_table1 %>% select(8,10)
seq01 <- unique(fas_table2 %>% select(6,7))
colnames(seq01)[2]<- "seq"
seq10 <- rbind(seq00,seq01)
seq10$start <- c(">")
seq10$spacer1 <- c("==")
seq10$spacer2 <- c("==")
seq11 <- seq10 %>% select(4,3,1,5,2)
write.table(seq11, quote = FALSE, col.names = FALSE,row.names = FALSE ,sep = ';', file = "formatA.txt")
za<- "formatA.txt"
za2<-
  readChar("formatA.txt", file.info("formatA.txt")$size)
za3<-
  gsub(";","",za2)


za4<-
  gsub("==","\r",za3)
write(za4, append = FALSE, file = paste(Accession,"_features.fas"))
###################################################################ALL COOL HERE - NOW GET IGS (GENE) AND INTRONS (CDSs)
