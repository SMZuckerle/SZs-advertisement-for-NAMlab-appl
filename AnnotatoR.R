#------This tool filters and formats stand alone BLAST results into geneBank-format ---------
#------The output can be copied and pasted into corresponding geneBank files and ------------
#------visualized with Snapgene--------------------------------------------------------------

setwd(
  "C:/Users/cyan4/OneDrive/Desktop/R_Home"
)

Filename= "results_table_Mpo_cds" #INPUT from BLAST result
GeneBank= "NC_017755"                        #INPUT for desired GB file to be annotated, used as query in blastn

treshold=100                                #set bistscore value for filtering blastn results

#---The BLAST_AnnotatoR.R Script annotates blastn-hits to a given genebank-file that has been used as blastn query---#
#---The outout file will be given the extension "NAME_blast_2_GB.txt---#
#---The content should look like this:  misc_feature (31462..32685)
# /note=BLAST_hit_Amborella_trichopoda_nad5i230g2
#---Copy and paste "NAME_blast_2_GB.txt content into corresponsing genebank-file for ducomentation                ---#

#--------------------------------------------------------- DO NOT CHANGE ANYTHING FROM HERE---

library(dplyr)
library(data.table)
library(tidyverse)

x<-read.csv(
  Filename
  , sep=",", header=FALSE, stringsAsFactors=FALSE)

x<-subset(x,
          V9 >=treshold)

x1<- x %>%
  select (2,1,5,6)

x1$strand <-
  ifelse(x1$V5 > x1$V6, "complement","x")


x1comp <- subset(x1, strand == "complement")
x1comp <- x1comp[,c(1,2,4,3,5)]
colnames(x1comp)[3]<- "V5"
colnames(x1comp)[4]<- "V6"
x1lead <- subset(x1, strand == "x")
x2 <- rbind(x1lead,x1comp)
x2$F5 <-
  ifelse(x1$strand == "complement", "(","x")
x2$F6 <-
  ifelse(x1$strand == "complement", ")","x")

#----------------------------------------------------
x2$F1 <-
  c("WW")
x2$F2 <-
  c("WX")
x2$F3 <-
  c("misc_feature")
x2$F4 <-
  c("#")
x2$F9 <-
  c("XX")
x2$F10 <-
  c("/note=BLAST_hit_")
x2$F13 <-
  c(",")
x2$F14<-
  c("..")

y <- x2[,c(8,9,10,11,5,6,3,15,4,7,8,12,13,2,8,1)]

#-------------------------------------------------------
z<- y[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)]
write.table(z, quote = FALSE, col.names = FALSE,row.names = FALSE ,sep = ';', file = "formatA.txt")

za<- "formatA.txt"
za2<-
  readChar("formatA.txt", file.info("formatA.txt")$size)
za3<-
  gsub(";","",za2)
za4<-
  gsub("WW","\r",za3)
za5<-
  gsub("WX","     ",za4)
za6<-
  gsub("XX","                     ",za5)
za7<-
  gsub("#","    ",za6)
za8<-
  gsub(", ","",za7)
za9<-
  gsub("x","",za8)
za10<-
  gsub(",,",",",za9)
write(za10, append = FALSE, file = paste(GeneBank, "_blast_2_GB.txt"))

#--------------------------------------------------------END-#
