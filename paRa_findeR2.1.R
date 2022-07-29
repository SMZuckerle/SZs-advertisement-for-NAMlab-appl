#A#Sort/ Filter Sequences based on Blast results#-----------------
#install.packages("ape")
#install.packages("phytools")
#install.packages("phylotools")
#install.packages("tidyr")
library(tidyr)
library(phytools)
library(phylotools)
library(ape)
library(dplyr)
library(data.table)

#---------------------------##
###PARA_FINDER2.0--------------##
#Para_finder uses blast-output table format to cluster sequences and estimation of their similarity.
# Its use is recommended for reciprocal Blast analyses with multi-query fasta files
# This might be especially useful when fasta sequences from different taxa and/or different loci have been blasted.
# The program will cluster fasta sequences based on the bitscore value, as basis for sequence similarity
# Different increments can be adjusted manually.
# Note that bitscore values heavily rely on sequence lenght.
# Further, it calculates percental likeliness values for each match above the increments, expressing the support for a cluster, based on actual matches/potential maximum matches.
#---------------------------##
##OUTPUT: PARA_FINDER produces four tables by default.
# 1. "_evaluated_table.txt" contains all matches of the blast analyses, respective the chosen increments
# 2. "_pFamilies_table.txt" contains a list of all clusters that have been produced respective the chosen increments
# 3. "_pFam_eval_table.txt" contains same as 2 with likeliness tag
# 4. "_match_count_table.txt" contains number of sequence presence, matches and likeliness-values per
#---------------------------##
#       INPUT             #
##CHOOSE working directory#
setwd(
  "G:/Machina_Eva/FoeRster"
)
##CHOOSE fas file for extraction of correct names. Para_finder can extract fasta headers macthing a group and prepare for sequence extraction with blastcmd#
Filename.fas= "mtg2_collection_structured_20210308.fas"


##CHOOSE .txt resulting from blast. The bistscore value should be placed in column 12. For best use please use standalone blast with the following setting: -outfmt "6"#
# (Later version might include use of distinct headers)
Filename= "results_test_table_mtg2_collection_structured_20210308.txt"

##CHOOSE treshold values. Blast matches above treshold values will be grouped and evaluated.
T0="99" #Treshold: Data below will not appear in later output. Use as filter, only
T1="149" #Set lowest treshold
T2="199" #Set medium treshold
T3="249" #Set highest treshold
#---------------------------##

#do not change from here on
################################################################################################
#GETTING-FASTA-NAMES: Preparation of list for later Extraction based on produced clusters-----##
fasta<-read.FASTA(Filename.fas)
names<-names(fasta)
name_frame_fas<-data.frame(names)

#write(name_frame_fas, append = FALSE, file = paste(name_frame_fas, "_all_names.txt"))

name_frame_fas$A = as.character(
  lapply(
    strsplit(
      as.character(name_frame_fas$names),
      split="_"),
    "[", 1)
)

name_frame_fas$B = as.character(
  lapply(
    strsplit(
      as.character(name_frame_fas$names),
      split="_"),
    "[", 2)
)

name_frame_fas$C = as.character(
  lapply(
    strsplit(
      as.character(name_frame_fas$names),
      split="_"),
    "[", 3)
)

#Filter by values after BLAST analsis#
#--------------------------------------------#----------------READING INPUT TABLE######
#Name results_table_DATA_blast = x
x<-read.csv(
  Filename, sep="\t", header=FALSE, stringsAsFactors=FALSE)

x$V12 <- as.numeric(as.character(x$V12))

#x <- x %>% select(1,2,9)

x<-as.data.table(
  x)[,sum(V12),by = .(V1,V2)]
names(x)<-make.unique(
  names(x))
write.csv(x,file=paste(Filename, "_aggregated.txt"), row.names = FALSE)
#V12 is filtering for bitscore-values#

#-------------------------------#

x$q1 = as.character(
  lapply(
    strsplit(
      as.character(x$V1),
      split="_"),
    "[", 1)
)
x$q2 = as.character(
  lapply(
    strsplit(
      as.character(x$V1),
      split="_"),
    "[", 2)
)
x$q3 = as.character(
  lapply(
    strsplit(
      as.character(x$V1),
      split="_"),
    "[", 3)
)


x$s1 = as.character(
  lapply(
    strsplit(
      as.character(x$V2),
      split="_"),
    "[", 1)
)
x$s2 = as.character(
  lapply(
    strsplit(
      as.character(x$V2),
      split="_"),
    "[", 2)
)
x$s3 = as.character(
  lapply(
    strsplit(
      as.character(x$V2),
      split="_"),
    "[", 3)
)
#------------COUNTING AND DATA PREPARATION- INDEXING -------------------#
split_list1<-split.data.frame(
  x,x$q1)

Qx <- as.data.frame(
  names(split_list1)
)
x_names <- names(split_list1)

Nx <- nrow(Qx)

Qx$Ny <- 1:nrow(Qx)

#colnames(Qx)[1] <- "q1"

#x<- merge.data.frame(x,Qx, by = "q1", all.x = TRUE, all.y = TRUE)

#colnames(Qx)[1] <- "s1"
#colnames(Qx)[2] <- "subjectID"

#x<- merge.data.frame(x,Qx, by = "s1", all.x = TRUE, all.y = TRUE)

colnames(x)[1] <- "Qcomplete"
colnames(x)[2] <- "Scomplete"
colnames(x)[3] <- "Score"
colnames(x)[4] <- "Qloc"
colnames(x)[5] <- "Qtax1"
colnames(x)[6] <- "Qtax2"
colnames(x)[7] <- "Sloc"
colnames(x)[8] <- "Stax1"
colnames(x)[9] <- "Stax2"

x$Score <- as.numeric(as.character(x$Score))

x$TT0<- ifelse(x$Score >T0,1,0)
x$TT1<- ifelse(x$Score >T1,1,0)
x$TT2<- ifelse(x$Score >T2,1,0)
x$TT3<- ifelse(x$Score >T3,1,0)

x$identLoc <- ifelse(x$Qloc==x$Sloc,1,0)
x$identTax <- ifelse(x$Qtax1==x$Stax1,1,0)


#Min-Max----------------------------------------------------------CLACUTATION OF MININUM AND MAXIMUM SCORES FOR EACH CLUSTER------------------#

XTorto<- subset(x,
                x$identLoc ==1)
#XTorto<- subset(XLorto,
#  XTorto$identTax ==0)


XTorto$Task<-paste(XTorto$Qloc, XTorto$Sloc, sep="_")

Task_ort <- XTorto %>%
  select(Task,Score)

Task_ort1 <- Task_ort %>%
  group_by(Task) %>%
  mutate(
    MaxScoreByTask = max(Score, na.rm = T),
    MinScoreByTask = min(Score, na.rm = T)
  ) %>%
  arrange(Task)
XLanal<- subset(x,
                x$identLoc ==0)

XLanal$Task<-paste(XLanal$Qloc, XLanal$Sloc, sep="_")

Task_Ana <- XLanal %>% select(Task,Score)

Task_Ana1 <- Task_Ana %>%
  group_by(Task) %>%
  mutate(
    MaxScoreByTask = max(Score, na.rm = T),
    MinScoreByTask = min(Score, na.rm = T)
  ) %>%
  arrange(Task)

y<- Task_ort1

y$Query = as.character(
  lapply(
    strsplit(
      as.character(y$Task),
      split="_"),
    "[", 1)
)

yX <- y %>% select(MinScoreByTask,MaxScoreByTask,Query,Task)
yX1 <- unique(yX)
colnames(yX1)[3] <-"Qloc"

yy<- merge.data.frame(x,yX1, by = "Qloc", all.x = TRUE, all.y = TRUE)


yy$Score <- as.numeric(as.character(yy$Score))
yy$MinScoreByTask <- as.numeric(as.character(yy$MinScoreByTask))
yy$MaxScoreByTask <- as.numeric(as.character(yy$MaxScoreByTask))

yy$T_MDO<- ifelse(yy$Score >=yy$MinScoreByTask,1,0)


XY<- yy %>%
  select (1,5,7,8,4,14,15,19,10,11,12,13,16,17,2,3)

sk<-Task_ort1 %>%
  select (1,3,4)

sk<-unique(sk)

#write.table(sk,file=paste(Filename, "_task_min_max.txt"), quote = FALSE, sep="\t", row.names = FALSE)

write.table(XY,file=paste(Filename, "_evaluated_table.txt"), quote = FALSE, sep="\t", row.names = FALSE)
#-----------------------------------------------------------------#--------------------------------------GENERATION OF CLUSTERS FOR OUTPUT----#

z<- XY %>% select(1:12)

#zP<- subset(z,
#             z$identLoc ==0)

zPT3<- z %>% select(1,3,12)
zPT3 <- unique(zPT3)

zPT3.1 <-as.data.table(
  zPT3)[,sum(TT3),by = .(Qloc,Sloc)]

zPT3.1$V1[zPT3.1$V1==0] <- NA

MT3<- zPT3.1 %>% spread(Qloc,V1)
MT3<- as.data.frame(MT3)
#MT3.1<-transform(MT3, char = as.numeric(char))
#MT3 <- as.numeric(as.character(MT3[2:162]))
#MT3 <- as.data.frame(MT3)
MT3.1<- ifelse(MT3[1:162] != 0, MT3$Sloc, NA)
MT3.2<- as.data.frame(MT3.1)
MT3.2<- MT3.2 %>% select(2:162)
MT3.2<- t(MT3.2)
MT3.2<- as.data.frame(MT3.2)
pFam_TT3<- unite(MT3.2, pFam_TT3, sep = ";", remove = TRUE, na.rm = TRUE)
##--------------------------------------------------------------------------##
zPT2<- z %>% select(1,3,11)
zPT2 <- unique(zPT2)

zPT2.1 <-as.data.table(
  zPT2)[,sum(TT2),by = .(Qloc,Sloc)]

zPT2.1$V1[zPT2.1$V1==0] <- NA

MT2<- zPT2.1 %>% spread(Qloc,V1)
MT2<- as.data.frame(MT2)
#MT3.1<-transform(MT3, char = as.numeric(char))
#MT3 <- as.numeric(as.character(MT3[2:162]))
#MT3 <- as.data.frame(MT3)
MT2.1<- ifelse(MT2[1:162] != 0, MT2$Sloc, NA)
MT2.2<- as.data.frame(MT2.1)
MT2.2<- MT2.2 %>% select(2:162)
MT2.2<- t(MT2.2)
MT2.2<- as.data.frame(MT2.2)
pFam_TT2<- unite(MT2.2, pFam_TT2, sep = ";", remove = TRUE, na.rm = TRUE)
##--------------------------------------------------------------------------##
zPT1<- z %>% select(1,3,10)
zPT1 <- unique(zPT1)

zPT1.1 <-as.data.table(
  zPT1)[,sum(TT1),by = .(Qloc,Sloc)]

zPT1.1$V1[zPT1.1$V1==0] <- NA

MT1<- zPT1.1 %>% spread(Qloc,V1)
MT1<- as.data.frame(MT1)
#MT3.1<-transform(MT3, char = as.numeric(char))
#MT3 <- as.numeric(as.character(MT3[2:162]))
#MT3 <- as.data.frame(MT3)
MT1.1<- ifelse(MT1[1:162] != 0, MT1$Sloc, NA)
MT1.2<- as.data.frame(MT1.1)
MT1.2<- MT1.2 %>% select(2:162)
MT1.2<- t(MT1.2)
MT1.2<- as.data.frame(MT1.2)
pFam_TT1<- unite(MT1.2, pFam_TT1, sep = ";", remove = TRUE, na.rm = TRUE)
##--------------------------------------------------------------------------##
zPMDO<- z %>% select(1,3,8)
zPMDO <- unique(zPMDO)

zPMDO.1 <-as.data.table(
  zPMDO)[,sum(T_MDO),by = .(Qloc,Sloc)]

zPMDO.1$V1[zPMDO.1$V1==0] <- NA

MMDO<- zPMDO.1 %>% spread(Qloc,V1)
MMDO<- as.data.frame(MMDO)
#MT3.1<-transform(MT3, char = as.numeric(char))
#MT3 <- as.numeric(as.character(MT3[2:162]))
#MT3 <- as.data.frame(MT3)
MMDO.1<- ifelse(MMDO[1:162] != 0, MT1$Sloc, NA)
MMDO.2<- as.data.frame(MMDO.1)
MMDO.2<- MMDO.2 %>% select(2:162)
MMDO.2<- t(MMDO.2)
MMDO.2<- as.data.frame(MMDO.2)
pFam_TMDO<- unite(MMDO.2, pFam_TMDO, sep = ";", remove = TRUE, na.rm = TRUE)
#--------------------------------------------------------------#
#pFamilies<-gsub("\n", "", pFamilies)
#pFamilies1<- pFamilies %>% select(1)
pFamilies1 <- merge(pFam_TT3, pFam_TT2, by=0, all=TRUE)
pFamilies2 <- merge(pFam_TT1, pFam_TMDO, by=0, all=TRUE)
pFamilies3 <- merge(pFamilies1, pFamilies2, by="row.names", all=TRUE)

pFamilies<- pFamilies3 %>% select(1,2,3,4,6,7)
colnames(pFamilies)[1] <-"ID"
colnames(pFamilies)[2] <-"Qloc"
write.table(pFamilies,file=paste(Filename, "_pFamilies_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

#-----------#CALCULATING Paralog_Scores Score how often two sequences match above treshold, divided by the number of max-number of matches sequences per blast#
c1<- XY %>% select(1,3,6,8,9,10,11,12)
c2<- XY %>% select(3,1,6,8,9,10,11,12)

c12 <- rbind(c1, c2)
c12<- subset(c12,
             c12$identLoc ==0)

c12mdo<- c12 %>% select(1,2,4)
Sc12mdo<-as.data.table(
  c12mdo)[,sum(T_MDO),by = .(Qloc,Sloc)]
Sc12mdo$tag<-paste(Sc12mdo$Qloc,Sc12mdo$Sloc, sep = "_")
colnames(Sc12mdo)[3] <-"T_MDO"

c12TT0<- c12 %>% select(1,2,5)
Sc12TT0<-as.data.table(
  c12TT0)[,sum(TT0),by = .(Qloc,Sloc)]
Sc12TT0$tag<-paste(Sc12TT0$Qloc,Sc12TT0$Sloc, sep = "_")
colnames(Sc12TT0)[3] <-"TT0"

c12TT1<- c12 %>% select(1,2,6)
Sc12TT1<-as.data.table(
  c12TT1)[,sum(TT1),by = .(Qloc,Sloc)]
Sc12TT1$tag<-paste(Sc12TT1$Qloc,Sc12TT1$Sloc, sep = "_")
colnames(Sc12TT1)[3] <-"TT1"

c12TT2<- c12 %>% select(1,2,7)
Sc12TT2<-as.data.table(
  c12TT2)[,sum(TT2),by = .(Qloc,Sloc)]
Sc12TT2$tag<-paste(Sc12TT2$Qloc,Sc12TT2$Sloc, sep = "_")
colnames(Sc12TT2)[3] <-"TT2"

c12TT3<- c12 %>% select(1,2,8)
Sc12TT3<-as.data.table(
  c12TT3)[,sum(TT3),by = .(Qloc,Sloc)]
Sc12TT3$tag<-paste(Sc12TT3$Qloc,Sc12TT3$Sloc, sep = "_")
colnames(Sc12TT3)[3] <-"TT3"

#c12$count <-c(1)
#c12count<- c12 %>% select(1,2,9)
#Sc12count<-as.data.table(
#  c12count)[,sum(count),by = .(Qloc,Sloc)]
#Sc12count$tag<-paste(Sc12count$Qloc,Sc12count$Sloc, sep = "_")
#colnames(Sc12count)[3] <-"C"
name_count<- table(name_frame_fas$A)
name_count.1<- as.data.frame(name_count)
colnames(name_count.1)[1] <-"Qloc"
colnames(name_count.1)[2] <-"Qcount"
c12$ID<-c(1:nrow(c12))
c12.1 <- merge(c12,name_count.1, by= "Qloc")
colnames(name_count.1)[1] <-"Sloc"
colnames(name_count.1)[2] <-"Scount"
c12.2 <- merge(c12.1,name_count.1, by= "Sloc")
c12.3<- c12.2 %>%select(2,1,3:ncol(c12.2))
c12.3$seq_count<- paste(c12.3$Qcount*c12.3$Scount)
c12.3$tag<-paste(c12.3$Qloc,c12.3$Sloc, sep = "_")


rSc12count <- c12.3 %>% select(13,12)
rSc12mdo<- Sc12mdo %>% select(3,4)
rSc12TT0<- Sc12TT0 %>% select(3,4)
rSc12TT1<- Sc12TT1 %>% select(3,4)
rSc12TT2<- Sc12TT2 %>% select(3,4)
rSc12TT3<- Sc12TT3 %>% select(3,4)

ScM <- merge(
  merge(
    merge(
      merge(
        merge(
          rSc12count,rSc12mdo,by="tag"),
        rSc12TT0, by = "tag"),
      rSc12TT1, by = "tag"),
    rSc12TT2, by = "tag"),
  rSc12TT3, by = "tag")

ScM2<-ScM[, 3:7]/2

ScM2$ID<-c(1:nrow(ScM))
ScM$ID<-c(1:nrow(ScM))
ScM3<- ScM %>% select(1,2,8)
ScM1.1<- merge(ScM3,ScM2, by = "ID")
ScM1.2<- ScM1.1 %>% select(2:8)
ScM1.2 <- unique(ScM1.2)
ScM1.2[, 2:7] <- sapply(ScM1.2[, 2:7], as.numeric)
ScM1.2$MDO_Pscore<-c(ScM1.2$T_MDO/ScM1.2$seq_count)
ScM1.2$TT0_Pscore<-c(ScM1.2$TT0/ScM1.2$seq_count)
ScM1.2$TT1_Pscore<-c(ScM1.2$TT1/ScM1.2$seq_count)
ScM1.2$TT2_Pscore<-c(ScM1.2$TT2/ScM1.2$seq_count)
ScM1.2$TT3_Pscore<-c(ScM1.2$TT3/ScM1.2$seq_count)

#C should be larger than one, else the match only appears once!
#COnsequently, total number of Qloci needs to counted

ScM2.1<- subset(ScM1.2, ScM1.2$seq_count>2)
ScM2.2<- subset(ScM2.1, ScM2.1$TT0!=0)
is.num <- sapply(ScM2.2, is.numeric)
ScM2.2[is.num] <- lapply(ScM2.2[is.num], round, 4)

ScM2.2$Qloc = as.character(
  lapply(
    strsplit(
      as.character(ScM2.2$tag),
      split="_"),
    "[", 1)
)

ScM2.2$Sloc = as.character(
  lapply(
    strsplit(
      as.character(ScM2.2$tag),
      split="_"),
    "[", 2)
)
#merge and write with evaluated table for manual in detail screening
#----------------------------------------------------------------------------------------------#-----INCLUDE CALCULATION INTO CLUSERS--------------#
ScM3.TT3<- ScM2.2 %>% select(13,14,12)
ScM3.TT3$Sloc_TT3<-paste(ScM3.TT3$Sloc,ScM3.TT3$TT3_Pscore, sep="_")
ScM3.TT3$TT3_Pscore[ScM3.TT3$TT3_Pscore==0] <- NA
ScM3.TT3.1<- ScM3.TT3 %>% select(1,4,3)
MScM3.TT3.1<- ScM3.TT3.1 %>% spread(Qloc,TT3_Pscore)
MScM3.TT3.1<- as.data.frame(MScM3.TT3.1)

MScM3.TT3.1<- ifelse(MScM3.TT3.1[1:ncol(MScM3.TT3.1)] != 0, MScM3.TT3.1$Sloc_TT3, NA)
MScM3.TT3.2<- as.data.frame(MScM3.TT3.1)
MScM3.TT3.2<- MScM3.TT3.2 %>% select(2:ncol(MScM3.TT3.2))
MScM3.TT3.2<- t(MScM3.TT3.2)
MScM3.TT3.2<- as.data.frame(MScM3.TT3.2)
pFam_MScM3.TT3.2<- unite(MScM3.TT3.2, pFam_TT3, sep = ";", remove = TRUE, na.rm = TRUE)

#------------------------------------------------------------------------------------------------#
ScM3.TT2<- ScM2.2 %>% select(13,14,11)
ScM3.TT2$Sloc_TT2<-paste(ScM3.TT2$Sloc,ScM3.TT2$TT2_Pscore, sep="_")
ScM3.TT2$TT2_Pscore[ScM3.TT2$TT2_Pscore==0] <- NA
ScM3.TT2.1<- ScM3.TT2 %>% select(1,4,3)
MScM3.TT2.1<- ScM3.TT2.1 %>% spread(Qloc,TT2_Pscore)
MScM3.TT2.1<- as.data.frame(MScM3.TT2.1)

MScM3.TT2.1<- ifelse(MScM3.TT2.1[1:ncol(MScM3.TT2.1)] != 0, MScM3.TT2.1$Sloc_TT2, NA)
MScM3.TT2.2<- as.data.frame(MScM3.TT2.1)
MScM3.TT2.2<- MScM3.TT2.2 %>% select(2:ncol(MScM3.TT2.2))
MScM3.TT2.2<- t(MScM3.TT2.2)
MScM3.TT2.2<- as.data.frame(MScM3.TT2.2)
pFam_MScM3.TT2.2<- unite(MScM3.TT2.2, pFam_TT2, sep = ";", remove = TRUE, na.rm = TRUE)
#------------------------------------------------------------------------------------------------#

ScM3.TT1<- ScM2.2 %>% select(13,14,10)
ScM3.TT1$Sloc_TT1<-paste(ScM3.TT1$Sloc,ScM3.TT1$TT1_Pscore, sep="_")
ScM3.TT1$TT1_Pscore[ScM3.TT1$TT1_Pscore==0] <- NA
ScM3.TT1.1<- ScM3.TT1 %>% select(1,4,3)
MScM3.TT1.1<- ScM3.TT1.1 %>% spread(Qloc,TT1_Pscore)
MScM3.TT1.1<- as.data.frame(MScM3.TT1.1)

MScM3.TT1.1<- ifelse(MScM3.TT1.1[1:ncol(MScM3.TT1.1)] != 0, MScM3.TT1.1$Sloc_TT1, NA)
MScM3.TT1.2<- as.data.frame(MScM3.TT1.1)
MScM3.TT1.2<- MScM3.TT1.2 %>% select(2:ncol(MScM3.TT1.2))
MScM3.TT1.2<- t(MScM3.TT1.2)
MScM3.TT1.2<- as.data.frame(MScM3.TT1.2)
pFam_MScM3.TT1.2<- unite(MScM3.TT1.2, pFam_TT1, sep = ";", remove = TRUE, na.rm = TRUE)

#------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------#

ScM3.T_MDO<- ScM2.2 %>% select(13,14,8)
ScM3.T_MDO$Sloc_T_MDO<-paste(ScM3.T_MDO$Sloc,ScM3.T_MDO$MDO_Pscore, sep="_")
ScM3.T_MDO$T_MDO_Pscore[ScM3.T_MDO$MDO_Pscore==0] <- NA
ScM3.T_MDO.1<- ScM3.T_MDO %>% select(1,4,3)
MScM3.T_MDO.1<- ScM3.T_MDO.1 %>% spread(Qloc,MDO_Pscore)
MScM3.T_MDO.1<- as.data.frame(MScM3.T_MDO.1)

MScM3.T_MDO.1<- ifelse(MScM3.T_MDO.1[1:ncol(MScM3.T_MDO.1)] != 0, MScM3.T_MDO.1$Sloc_T_MDO, NA)
MScM3.T_MDO.2<- as.data.frame(MScM3.T_MDO.1)
MScM3.T_MDO.2<- MScM3.T_MDO.2 %>% select(2:ncol(MScM3.T_MDO.2))
MScM3.T_MDO.2<- t(MScM3.T_MDO.2)
MScM3.T_MDO.2<- as.data.frame(MScM3.T_MDO.2)
pFam_MScM3.T_MDO.2<- unite(MScM3.T_MDO.2, pFam_T_MDO, sep = ";", remove = TRUE, na.rm = TRUE)

#------------------------------------------------------------------------------------------------#
write.table(ScM2.2,file=paste(Filename, "_match-count_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)
#write.table(pFam_MScM3.T_MDO.2,file=paste(Filename, "_MDO-fam_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
#write.table(pFam_MScM3.TT1.2,file=paste(Filename, "_TT1-fam_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
#write.table(pFam_MScM3.TT2.2,file=paste(Filename, "_TT2-fam_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
#write.table(pFam_MScM3.TT3.2,file=paste(Filename, "_TT3-fam_table.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = TRUE)
#------------------------------------------------------------------------------------------------#

pFam_MScM3.T_MDO.2$Fam_set<-c("TMDO")
pFam_MScM3.T_MDO.2$Qloc<-c(row.names(pFam_MScM3.T_MDO.2))
pFam_MScM3.T_MDO.3<- pFam_MScM3.T_MDO.2 %>% select(2,3,1)
colnames(pFam_MScM3.T_MDO.3)[3] <-"pPara"
pFam_MScM3.TT1.2$Fam_set<-c("TT1")
pFam_MScM3.TT1.2$Qloc<-c(row.names(pFam_MScM3.TT1.2))
pFam_MScM3.TT1.3<- pFam_MScM3.TT1.2 %>% select(2,3,1)
colnames(pFam_MScM3.TT1.3)[3] <-"pPara"
pFam_MScM3.TT2.2$Fam_set<-c("TT2")
pFam_MScM3.TT2.2$Qloc<-c(row.names(pFam_MScM3.TT2.2))
pFam_MScM3.TT2.3<- pFam_MScM3.TT2.2 %>% select(2,3,1)
colnames(pFam_MScM3.TT2.3)[3] <-"pPara"
pFam_MScM3.TT3.2$Fam_set<-c("TT3")
pFam_MScM3.TT3.2$Qloc<-c(row.names(pFam_MScM3.TT3.2))
pFam_MScM3.TT3.3<- pFam_MScM3.TT3.2 %>% select(2,3,1)
colnames(pFam_MScM3.TT3.3)[3] <-"pPara"
pFam_join1 <- rbind(pFam_MScM3.T_MDO.3, pFam_MScM3.TT1.3, pFam_MScM3.TT2.3, pFam_MScM3.TT3.3, row.names=FALSE)


pFam_join1.1<-pFam_join1 %>% mutate_all(na_if,"")
write.table(pFam_join1.1,file=paste(Filename, "_pFam_eval_table2.txt"), quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)


############################DONE###########################################################

#       USAGE FROM HERE ON IS OPTIONAL
#       IDENTIFIED CLUSTERS ABOVE CAN BE WRITTEN INTO INDEXES FOR BLASTCMD EXTRACTION OF FASTA-SEQUENCES

#Get groups of names and values#----------------------------------#


names1 <- names(split_list1)

for(
  i in seq_along(names1)
){
  assign(
    names1[i], split_list1[[i]])
}

##############MAGIC: CHOSE THE LOCUS CLUSTER

Cluster="trnSi39g2"

Cluster_names <- data.frame(Cluster$s3)
colnames(Cluster_names)[1] <-"s3"
colnames(name_frame_fas)[4] <-"s3"

desired_fastas<-semi_join(
  name_frame_fas,Cluster_names, (by="s3")
)

desired_fastas_alone <-desired_fastas[1]
list_test <-list(desired_fastas_alone)
list(desired_fastas_alone)
write(desired_fastas_alone$names, file= paste(Cluster, "_cluster.txt"))

#GO TO G:\Machina_Eva\BLASTplus\ncbi-blast-2.10.1+\bin#

#blastdbcmd -db mtg2_collection_20210129.fas -dbtype nucl -entry_batch cluster.txt -outfmt "%f" -out cluster.fas


#B#Extract subsets of sequences based on Filter step#-------------

#load fasta#

#get fasta names#

#filter for fasta names meeting criteria of step A#

#C#Write Fasta - Align Fasta with mafft in console - Generate phyloTree in IQtree###
