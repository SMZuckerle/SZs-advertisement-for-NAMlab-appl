#MAPPING-EVALUATOR FOR GMAP FILE OUTPUT#----------------------------------#
library(dplyr)
library(data.table)
library(tidyverse)
#---------SETTINGS--------------------------------------------------------#
###########################################################################
#-------------------------------------------------------------------------#

setwd(
  "F:/Germany/Paderbonn-ubo1001/NGS_Dec2019/R script")

#---For RNA evaluation remove line 1-10 from input file-------------------#
#---Set treshold for partial RNA-editing         -------------------------#
#---Script only runs on one fasta as reference correctly for GB annotation#
#---Split gmap .out file into sub-data files i.c. of multiple references--#

RNAediting_cutoff= c(0.029)

#---Enter Mapping reference or comments-----------------------------------#

Comm1= c(",g17")    #"," as seperator; "CP" determine Ref. seq; "g20" for gmap-2020

#---Search and Replace pattern "filename" with new "filename"-------------#

Filename= "GSNAP_Jaccusa_=PPR56_1.2_test.out"

Filterset= "PPR56_interested.txt"



#-------------------------------------------------------------------------#
###########################################################################
#-------------------------------------------------------------------------#

x<-read.csv(
  Filename
  , sep="\t", header=FALSE, stringsAsFactors=FALSE)
z<-read.csv(Filterset,sep = "\r", header=FALSE)

#---------Formatting------------------------------------------------------#
###########################################################################
#-------------------------------------------------------------------------#
x1<- x %>%
  select (1,2,4,5,7,9)
xa<- x %>%
  select (2,10)
xb<- x %>%
  select (2,11)

xa1<-
  xa %>% 
  separate(2, c("DNAtotal","col2"),":")
xb1<-
  xb %>% 
  separate(2, c("RNAtotal","col2"),":")
xa2<-
  xa1 %>% 
  separate(3, c("dA","dC","dG","dT"),",")
xb2<-
  xb1 %>% 
  separate(3, c("A","C","G","U"),",")

x2<-
  merge.data.frame(x1,xa2, by = "V2", all.x = TRUE, all.y = TRUE)
x3<-
  merge.data.frame(x2,xb2, by = "V2", all.x = TRUE, all.y = TRUE)

x3$V2 <-
  as.numeric(as.character(x3$V2))
x3$DNAtotal <-
  as.numeric(as.character(x3$DNAtotal))
x3$RNAtotal <-
  as.numeric(as.character(x3$RNAtotal))
x3$dA <-
  as.numeric(as.character(x3$dA))
x3$A <-
  as.numeric(as.character(x3$A))
x3$dC <-
  as.numeric(as.character(x3$dC))
x3$C <-
  as.numeric(as.character(x3$C))
x3$dG <-
  as.numeric(as.character(x3$dG))
x3$G <-
  as.numeric(as.character(x3$G))
x3$dT <-
  as.numeric(as.character(x3$dT))
x3$U <-
  as.numeric(as.character(x3$U))

#---------Calculations------------------------------------------------------#
x3$Cquality <-
  x3$dC/x3$DNAtotal
x3$Tquality <-
  x3$dT/x3$DNAtotal
x3$Aquality <-
  x3$dA/x3$DNAtotal
x3$Gquality <-
  x3$dG/x3$DNAtotal

x3$Ugain <-
  x3$U/x3$RNAtotal-x3$dT/x3$DNAtotal
x3$Cgain <-
  x3$C/x3$RNAtotal-x3$dC/x3$DNAtotal
x3$Ggain <-
  x3$G/x3$RNAtotal-x3$dG/x3$DNAtotal
x3$Again <-
  x3$A/x3$RNAtotal-x3$dA/x3$DNAtotal

x3$bias1 <-
  x3$Ugain+x3$Cgain
x3$bias2 <-
  x3$Again+x3$Ggain
x3$quick <-
  x3$Cquality*x3$Ugain+x3$Tquality*x3$Cgain+x3$Aquality*x3$Ggain+x3$Gquality*x3$Again

#---------Calculation-FLAGS----------------------------------------------------------

x3$cutoff <-
  ifelse(x3$quick > RNAediting_cutoff,1,0)
x3$eU <-
  ifelse(x3$Cquality+x3$Ugain+x3$Gquality+x3$Again>1,"eU","x")
x3$eC <-
  ifelse(x3$Tquality+x3$Cgain+x3$Aquality+x3$Ggain>1,"eC","x")
x3$partial <-
  ifelse(x3$quick > 0.5, "partial","x")
x3$comp1 <-
  ifelse(x3$Ggain >= 0.1, 1, 0)
x3$comp2 <-
  ifelse(x3$Again >= 0.1, 1, 0)
x3$complement <-
  ifelse(x3$comp1+x3$comp2 >= 1,"complement","x")

x3$commentary <- c(Comm1)

x4 <-
  subset.data.frame(x3, x3$cutoff>=1)

colnames(z)[1] <- "Pos"

colnames(x4)[1] <- "Pos"

test<-
  merge.data.frame(z,x4, by = "Pos", all.x = FALSE, all.y = FALSE)

#---------Clean-up--------------------------------------------------------------------

write.table(x4,file=paste(Filename, "_evaluated_table.txt"), quote = FALSE, sep="\t")

write.table(test,file=paste(Filename, "_filtered_table.txt"), quote = FALSE, sep="\t")

#---------Prepare for gb. Import-------------------------------------------------------
F1 <-
  c("Ö")
F2 <-
  c("Ä")
F3 <-
  c("misc_feature")
F4 <-
  c("#")
F9 <-
  c("Ü")
F10 <-
  c("/note=RNA_editing_")
F13 <-
  c(",")
y <- test %>%
  select (1,29,30,31,34,35)
y$F1 <- F1
y$F2 <- F2
y$F3 <- F3
y$F4 <- F4
y$F8 <- F1
y$F9 <- F9
y$F13 <- F13
y$F16 <- F1
y$F5 <- 
  ifelse(y$complement == "complement", "(","x")
y$F7 <- 
  ifelse(y$complement == "complement", ")","x")
y$F10 <- F10
y2 <- y[,c(7,8,9,10,5,15,1,16,10,11,12,17,2,3,13,4,6,14)]
write.table(y2, quote = FALSE, col.names = FALSE,row.names = FALSE ,sep = ';', file = "format.txt")
z<- "format.txt"
z2<- 
  readChar("format.txt", file.info("format.txt")$size)
z3<- 
  gsub(";","",z2)
z4<- 
  gsub("Ö","\r",z3)
z5<- 
  gsub("Ä","     ",z4)
z6<- 
  gsub("Ü","                     ",z5)
z7<- 
  gsub("#","    ",z6)
z8<- 
  gsub(", ","",z7)
z9<- 
  gsub("x","",z8)
z10<- 
  gsub(",,",",",z9)
write(z10, append = FALSE, file = paste(Filename, "_for_GB.txt"))

#---------DONE-------------------------------------------------------
#####################################################################
#--------------------------------------------------------------------

