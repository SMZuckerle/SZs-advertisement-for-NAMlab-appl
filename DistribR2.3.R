#--Analyse and structure binary character data (e.g. presence/absence)--#
#--INPUT data should be structured as table with characters as rows ----#
#--and taxa as columns--------------------------------------------------#
#--Script will produce condensed subsets for overview, usage in Venn ---#
#--analyses and plots --------------------------------------------------#

setwd(
  ""
)

#
#----------------------------------------------------------------------#
Filename= "mtg2_overview2021_4.csv"
F_order="mtg2_overview2021_F_order2.csv"
#----------------------------------------------------------------------#
x<-
  read.csv(
    Filename, header = TRUE, stringsAsFactors = FALSE, sep = ";"
  )

Forder<-
  read.csv(
    F_order, header = TRUE, stringsAsFactors = FALSE, sep = ","
  )
colnames(Forder)[1]<-"intron"
colnames(Forder)[2]<-"id"

#----------------------------------------------------------------------#
install.packages(eulerr)
install.packages(tidyverse)
install.packages(dplyr)
library(dplyr)
library(tidyverse)
library(eulerr)
#library(data.table)
#library(ggplot2)
#library(reshape2)


#----------------------------------------------------------------------#


x1 <- x %>%
  select (1,2,4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188, 192)


x2 <- x %>%
  select (1,4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188, 192)
colnames(x2)[1]<-"locus"
write.table(x1, file = paste(Filename, "_presence.txt"), quote = FALSE, row.names = FALSE, sep = ",")


x3 <- x %>%
  select(1,2)
xT1 <- x %>%
  select(1,2)
xC <- x %>%
  select(1,2)
#AAL-ancestral streptophyte algae, polyphyletic; Mesostigma, Chlorokybus, Klebmordium
x3$AAL<-
  ifelse(x1[3]+x1[4]+x1[5]==0,0,1)
xT1$AAL<-
  (x1[3]+x1[4]+x1[5])/3
xC$AAL<-
  (x1[3]+x1[4]+x1[5])
#CHA - Charophytes, monophyletic; Chara, Nitella
x3$CHA<-
  ifelse(x1[6]+x1[7]==0,0,1)
xT1$CHA<-
  (x1[6]+x1[7])/2
xC$CHA<-
  (x1[6]+x1[7])
#CZA;COL-ZYG - sister group to EMB, appears monophyletic in cpDNA based data,
x3$CZA<- 
  ifelse(x1[8]+x1[9]+x1[10]+x1[11]+x1[12]+x1[13]+x1[14]==0,0,1)
xT1$CZA<- 
  (x1[8]+x1[9]+x1[10]+x1[11]+x1[12]+x1[13]+x1[14])/7
xC$CZA<- 
  (x1[8]+x1[9]+x1[10]+x1[11]+x1[12]+x1[13]+x1[14])

x3$COL<-
  ifelse(x1[8]+x1[9]==0,0,1)
xT1$COL<- 
  (x1[8]+x1[9])/2
xC$COL<- 
  (x1[8]+x1[9])

x3$ZYG<- 
  ifelse(x1[10]+x1[11]+x1[12]+x1[13]+x1[14]==0,0,1)
xT1$ZYG<- 
  (x1[10]+x1[11]+x1[12]+x1[13]+x1[14])/5
xC$ZYG<- 
  (x1[10]+x1[11]+x1[12]+x1[13]+x1[14])

ALG<- x1 %>%
  select(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
CHACZA<- x1 %>%
  select(6,7,8,9,10,11,12,13,14)

x3$LIV<- 
  ifelse(x1[15]+x1[16]+x1[17]+x1[18]+x1[19]+x1[20]==0,0,1)
xT1$LIV<- 
  (x1[15]+x1[16]+x1[17]+x1[18]+x1[19]+x1[20])/6
xC$LIV<- 
  (x1[15]+x1[16]+x1[17]+x1[18]+x1[19]+x1[20])
LIV<- x1 %>%
  select(15,16,17,18,19,20)

x3$MOS<- 
  ifelse(x1[21]+x1[22]+x1[23]+x1[24]+x1[25]+x1[26]==0,0,1)
xT1$MOS<- 
  (x1[21]+x1[22]+x1[23]+x1[24]+x1[25]+x1[26])/6
xC$MOS<- 
  (x1[21]+x1[22]+x1[23]+x1[24]+x1[25]+x1[26])
MOS<- x1 %>%
  select(21,22,23,24,25,26)

x3$HOR<- 
  ifelse(x1[27]+x1[28]+x1[29]+x1[30]+x1[31]==0,0,1)
xT1$HOR<- 
  (x1[27]+x1[28]+x1[29]+x1[30]+x1[31])/5
xC$HOR<- 
  (x1[27]+x1[28]+x1[29]+x1[30]+x1[31])
HOR<- x1 %>%
  select(27,28,29,30,31)

x3$LYC<- 
  ifelse(x1[32]+x1[33]+x1[34]==0,0,1)
xT1$LYC<- 
  (x1[32]+x1[33]+x1[34])/3
xC$LYC<- 
  (x1[32]+x1[33]+x1[34])
LYC<- x1 %>%
  select(32,33,34)

x3$MON<- 
  ifelse(x1[35]+x1[36]+x1[37]+x1[38]==0,0,1)
xT1$MON<- 
  (x1[35]+x1[36]+x1[37]+x1[38])/4
xC$MON<- 
  (x1[35]+x1[36]+x1[37]+x1[38])
x3$GYM<- 
  ifelse(x1[39]+x1[40]+x1[41]+x1[42]+x1[43]==0,0,1)
xT1$GYM<- 
  (x1[39]+x1[40]+x1[41]+x1[42]+x1[43])/5
xC$GYM<- 
  (x1[39]+x1[40]+x1[41]+x1[42]+x1[43])
x3$ANG<- 
  ifelse(x1[44]+x1[45]+x1[46]+x1[47]+x1[48]+x1[49]+x1[50]==0,0,1)
xT1$ANG<- 
  (x1[44]+x1[45]+x1[46]+x1[47]+x1[48]+x1[49]+x1[50])/7
xC$ANG<- 
  (x1[44]+x1[45]+x1[46]+x1[47]+x1[48]+x1[49]+x1[50])

TRA<- x1 %>%
  select(32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)
EUP<- x1 %>%
  select(35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)

x3$ALG<-
  ifelse(x3$AAL+x3$CHA+x3$COL+x3$ZYG==0,0,1)
x3$BRY<-
  ifelse(x3$LIV+x3$MOS+x3$HOR==0,0,1)
x3$TRA<-
  ifelse(x3$LYC+x3$MON+x3$GYM+x3$ANG==0,0,1)
x3$EUP<-
  ifelse(x3$MON+x3$GYM+x3$ANG==0,0,1)
x3$EMB<-
  ifelse(x3$MON+x3$GYM+x3$ANG==0,0,1)

xT2<- xT1 %>%
  select(1,2,3,4,6,7,8,9,10,11,12,13,14)

xC2<- xC %>%
  select(1,2,3,4,6,7,8,9,10,11,12,13,14)

xC2$sm<- 
  (xC2[3]+xC2[4]+xC2[5]+xC2[6]+xC2[7]+xC2[8]+xC2[9]+xC2[10]+xC2[11]+xC2[12]+xC2[13])

xC2<- xC2 %>%
  select(1,2,14,3,4,5,6,7,8,9,10,11,12,13)

xCsm<- xC2 %>%
  select(1,3)
colnames(xT2)[1]<- "intron"
colnames(xCsm)[1]<- "intron"
xTC1<- merge.data.frame(xT2,xCsm, by = "intron", all.x = TRUE, all.y = TRUE)

xTC2<-as.data.table((xTC1))
write.table(xTC2, file = paste(Filename, "_sm_gainloss_percent.txt"), quote = FALSE, row.names = FALSE, sep = ",")
xC3<-as.data.table((xC2))
xT2<-as.data.table((xT2))
write.table(xT2, file = paste(Filename, "_gainloss_percent.txt"), quote = FALSE, row.names = FALSE, sep = ",")


#-----WE ARE COOL UNTIL HERE----------------------------####
LIV2<-as.data.frame(t(LIV))
MOS2<-as.data.frame(t(MOS))
HOR2<-as.data.frame(t(HOR))
TRA2<-as.data.frame(t(TRA))
LYC2<-as.data.frame(t(LYC))
EUP2<-as.data.frame(t(EUP))
ALG2<-as.data.frame(t(ALG))
CHACZA2<-as.data.frame(t(CHACZA))

write.table(LIV2, file = "LIV2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(MOS2, file = "MOS2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(HOR2, file = "HOR2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(TRA2, file = "TRA2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(LYC2, file = "LYC2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(EUP2, file = "EUP2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(ALG2, file = "ALG2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")
write.table(CHACZA2, file = "CHACZA2_distrib_202110.txt", append = FALSE, quote = FALSE, row.names = TRUE, col.names = FALSE ,sep = " ")



colnames(x3)[1]<-"intron"
colnames(x3)[2]<-"family"
x3[is.na(x3)]<-0
x3$intron<- as.character(x3$intron)
x3$family<- as.character(x3$family)
x3$intron<-factor(x3$intron,levels = x3$intron)
x3 %>% map_df(rev)

write.table(x3, file = paste(Filename, "_lineages.txt"), quote = FALSE, row.names = FALSE, sep = ",")

#----------------------------------------------------------------------#
x4<- x3 %>%
  select(1,2,3,4,6,7,8,9,10,11,12,13,14,15)
xF<- 
  subset.data.frame(x4,x4$family!="S")



#xF2<- xFo %>%
#  select(11,2,3,4,5,6,7,8,9,10)

xS<- 
  subset.data.frame(x4,x4$family=="S")

#x5<- x3 %>%
#  select(1,2,3,4,6,7,8,9,10,11,12,13,14)

#x5$ALG<-
#  ifelse(x5[3]+x5[4]+x5[5]==0,0,1)
#x5$BRY<-
#  ifelse(x5[6]+x5[7]+x5[8]==0,0,1)
#x5$TRA<-
#  ifelse(x5[9]+x5[10]+x5[11]+x5[12]==0,0,1)
#x5$EMB<-
#  ifelse(x5[6]+x5[7]+x5[8]+x5[9]+x5[10]+x5[11]+x5[12]==0,0,1)

xFagg<- as.data.frame(xF %>%
  select(2:13))
#second = function(x) {
#  if (length(x) == 1)
#    return(x)
#  return(sort(x, decreasing = TRUE)[2])}
#xFagg1<- aggregate(xFagg, by = list (xFagg$family), FUN = sum)
xFagg1<-xFagg %>%
  group_by(family) %>%
  summarise_at(vars(AAL,CHA,COL,ZYG,LIV,MOS,HOR,LYC,MON,GYM,ANG),
               list(name = sum))

#xFagg1<- xFagg1 %>%
#  select(2:13)
#xFagg1$F_<-c("F")
#xFagg1$familyF<-paste(xFagg1$F_,xFagg1$family, sep = "")
#xFagg1<- xFagg1 %>%
#  select(14,2:12)

#xFagg1<-as.data.table(
#  xFagg)[,sum(xFagg[2:12]),by = .(xFagg[1])]


#------TRANSFORM DISTRIBUTION INTO VENN----------------------------------------------------------------#

y1<- x3

y1$AAL<- ifelse(y1$AAL==1,y1$intron,NA)
y1$CHA<- ifelse(y1$CHA==1,y1$intron,NA)
y1$CZA<- ifelse(y1$CZA==1,y1$intron,NA)
y1$COL<- ifelse(y1$COL==1,y1$intron,NA)
y1$ZYG<- ifelse(y1$ZYG==1,y1$intron,NA)
y1$LIV<- ifelse(y1$LIV==1,y1$intron,NA)
y1$MOS<- ifelse(y1$MOS==1,y1$intron,NA)
y1$HOR<- ifelse(y1$HOR==1,y1$intron,NA)
y1$LYC<- ifelse(y1$LYC==1,y1$intron,NA)
y1$MON<- ifelse(y1$MON==1,y1$intron,NA)
y1$GYM<- ifelse(y1$GYM==1,y1$intron,NA)
y1$ANG<- ifelse(y1$ANG==1,y1$intron,NA)
y1$ALG<- ifelse(y1$ALG==1,y1$intron,NA)
y1$BRY<- ifelse(y1$BRY==1,y1$intron,NA)
y1$TRA<- ifelse(y1$TRA==1,y1$intron,NA)
y1$EMB<- ifelse(y1$EMB==1,y1$intron,NA)
y1$EUP<- ifelse(y1$EUP==1,y1$intron,NA)
x3$id  <- 1:nrow(x3)

ids <- x3 %>%
  select(1,20)

ids$intron<- as.character(ids$intron)

ids$AAL<- ifelse(y1$AAL==ids$id,ids$intron,NA)
ids$CHA<- ifelse(y1$CHA==ids$id,ids$intron,NA)
ids$CZA<- ifelse(y1$CZA==ids$id,ids$intron,NA)
ids$COL<- ifelse(y1$COL==ids$id,ids$intron,NA)
ids$ZYG<- ifelse(y1$ZYG==ids$id,ids$intron,NA)
ids$LIV<- ifelse(y1$LIV==ids$id,ids$intron,NA)
ids$MOS<- ifelse(y1$MOS==ids$id,ids$intron,NA)
ids$HOR<- ifelse(y1$HOR==ids$id,ids$intron,NA)
ids$LYC<- ifelse(y1$LYC==ids$id,ids$intron,NA)
ids$MON<- ifelse(y1$MON==ids$id,ids$intron,NA)
ids$GYM<- ifelse(y1$GYM==ids$id,ids$intron,NA)
ids$ANG<- ifelse(y1$ANG==ids$id,ids$intron,NA)

ids$ALG<- ifelse(y1$ALG==ids$id,ids$intron,NA)
ids$BRY<- ifelse(y1$BRY==ids$id,ids$intron,NA)
ids$TRA<- ifelse(y1$TRA==ids$id,ids$intron,NA)
ids$EMB<- ifelse(y1$EMB==ids$id,ids$intron,NA)

write.table(ids, file = paste(Filename, "_for_VENN.txt"), quote = FALSE, row.names = FALSE, sep = ",")


#----------------------------------------------------------------#
#y2<- y1 %>%
#  select(2,3,4,5,6,7,8,9,10,11,12,13)


y2<-lapply(y1,as.numeric)

##continue here
AAL<- paste0(y2$AAL[!is.na(y2$AAL)], collapse = " ")
CHA<- paste0(y2$CHA[!is.na(y2$CHA)], collapse = " ")
CZA<- paste0(y2$CZA[!is.na(y2$CZA)], collapse = " ")
#ALG<- paste0(y2$ALG[!is.na(y2$ALG)], collapse = " ")
LIV<- paste0(y2$LIV[!is.na(y2$LIV)], collapse = " ")
MOS<- paste0(y2$MOS[!is.na(y2$MOS)], collapse = " ")
HOR<- paste0(y2$HOR[!is.na(y2$HOR)], collapse = " ")
LYC<- paste0(y2$LYC[!is.na(y2$LYC)], collapse = " ")
MON<- paste0(y2$MON[!is.na(y2$MON)], collapse = " ")
GYM<- paste0(y2$GYM[!is.na(y2$GYM)], collapse = " ")
ANG<- paste0(y2$ANG[!is.na(y2$ANG)], collapse = " ")

ALG<- paste0(y2$ALG[!is.na(y2$ALG)], collapse = " ")
TRA<- paste0(y2$TRA[!is.na(y2$TRA)], collapse = " ")
BRY<- paste0(y2$BRY[!is.na(y2$BRY)], collapse = " ")
EMB<- paste0(y2$EMB[!is.na(y2$EMB)], collapse = " ")

AAL<-as.numeric(strsplit(AAL, " ")[[1]])
CHA<-as.numeric(strsplit(CHA, " ")[[1]])
CZA<-as.numeric(strsplit(CZA, " ")[[1]])
#ALG<-as.numeric(strsplit(ALG, " ")[[1]])
LIV<-as.numeric(strsplit(LIV, " ")[[1]])
MOS<-as.numeric(strsplit(MOS, " ")[[1]])
HOR<-as.numeric(strsplit(HOR, " ")[[1]])
LYC<-as.numeric(strsplit(LYC, " ")[[1]])
MON<-as.numeric(strsplit(MON, " ")[[1]])
GYM<-as.numeric(strsplit(GYM, " ")[[1]])
ANG<-as.numeric(strsplit(ANG, " ")[[1]])

ALG<-as.numeric(strsplit(ALG, " ")[[1]])
TRA<-as.numeric(strsplit(TRA, " ")[[1]])
BRY<-as.numeric(strsplit(BRY, " ")[[1]])
EMB<-as.numeric(strsplit(EMB, " ")[[1]])
#-----PLOTTING VENN DIAGRAMS AND EULER-----------------------------------------------------------------#
FAG<- xFagg1

FAG$AAL<- ifelse(FAG$AAL_name!=0,FAG$family,NA)
FAG$CHA<- ifelse(FAG$CHA_name!=0,FAG$family,NA)
FAG$COL<- ifelse(FAG$COL_name!=0,FAG$family,NA)
FAG$ZYG<- ifelse(FAG$ZYG_name!=0,FAG$family,NA)
FAG$LIV<- ifelse(FAG$LIV_name!=0,FAG$family,NA)
FAG$MOS<- ifelse(FAG$MOS_name!=0,FAG$family,NA)
FAG$HOR<- ifelse(FAG$HOR_name!=0,FAG$family,NA)
FAG$LYC<- ifelse(FAG$LYC_name!=0,FAG$family,NA)
FAG$MON<- ifelse(FAG$MON_name!=0,FAG$family,NA)
FAG$GYM<- ifelse(FAG$GYM_name!=0,FAG$family,NA)
FAG$ANG<- ifelse(FAG$ANG_name!=0,FAG$family,NA)

FAG$ALG<- 
  ifelse(FAG[2]+FAG[3]+FAG[4]+FAG[5]==0,NA,FAG$family)
FAG$BRY<- 
  ifelse(FAG[6]+FAG[7]+FAG[8]==0,NA,FAG$family)
FAG$TRA<- 
  ifelse(FAG[9]+FAG[10]+FAG[11]+FAG[12]==0,NA,FAG$family)
FAG$EMB<- 
  ifelse(FAG[6]+FAG[7]+FAG[8]+FAG[9]+FAG[10]+FAG[11]+FAG[12]==0,NA,FAG$family)

#FAG$ALG<- ifelse(FAG$ALG_name!=0,FAG$family,NA)
#FAG$BRY<- ifelse(FAG$BRY_name!=0,FAG$family,NA)
#FAG$TRA<- ifelse(FAG$TRA_name!=0,FAG$family,NA)
#FAG$EMB<- ifelse(FAG$EMB_name!=0,FAG$family,NA)
FAG2<- FAG %>%
  select(1,13:23)
FAG2<-lapply(FAG,as.numeric)

##continue here
fAAL<- paste0(FAG2$AAL[!is.na(FAG2$AAL)], collapse = " ")
fCHA<- paste0(FAG2$CHA[!is.na(FAG2$CHA)], collapse = " ")
fCOL<- paste0(FAG2$COL[!is.na(FAG2$COL)], collapse = " ")
fZYG<- paste0(FAG2$ZYG[!is.na(FAG2$ZYG)], collapse = " ")
#ALG<- paste0(FAG2$ALG[!is.na(FAG2$ALG)], collapse = " ")
fLIV<- paste0(FAG2$LIV[!is.na(FAG2$LIV)], collapse = " ")
fMOS<- paste0(FAG2$MOS[!is.na(FAG2$MOS)], collapse = " ")
fHOR<- paste0(FAG2$HOR[!is.na(FAG2$HOR)], collapse = " ")
fLYC<- paste0(FAG2$LYC[!is.na(FAG2$LYC)], collapse = " ")
fMON<- paste0(FAG2$MON[!is.na(FAG2$MON)], collapse = " ")
fGYM<- paste0(FAG2$GYM[!is.na(FAG2$GYM)], collapse = " ")
fANG<- paste0(FAG2$ANG[!is.na(FAG2$ANG)], collapse = " ")

fALG<- paste0(FAG2$ALG[!is.na(FAG2$ALG)], collapse = " ")
fTRA<- paste0(FAG2$TRA[!is.na(FAG2$TRA)], collapse = " ")
fBRY<- paste0(FAG2$BRY[!is.na(FAG2$BRY)], collapse = " ")
fEMB<- paste0(FAG2$EMB[!is.na(FAG2$EMB)], collapse = " ")

fAAL<-as.numeric(strsplit(fAAL, " ")[[1]])
fCHA<-as.numeric(strsplit(fCHA, " ")[[1]])
fCOL<-as.numeric(strsplit(fCOL, " ")[[1]])
fZYG<-as.numeric(strsplit(fZYG, " ")[[1]])
#ALG<-as.numeric(strsplit(ALG, " ")[[1]])
fLIV<-as.numeric(strsplit(fLIV, " ")[[1]])
fMOS<-as.numeric(strsplit(fMOS, " ")[[1]])
fHOR<-as.numeric(strsplit(fHOR, " ")[[1]])
fLYC<-as.numeric(strsplit(fLYC, " ")[[1]])
fMON<-as.numeric(strsplit(fMON, " ")[[1]])
fGYM<-as.numeric(strsplit(fGYM, " ")[[1]])
fANG<-as.numeric(strsplit(fANG, " ")[[1]])

fALG<-as.numeric(strsplit(fALG, " ")[[1]])
fTRA<-as.numeric(strsplit(fTRA, " ")[[1]])
fBRY<-as.numeric(strsplit(fBRY, " ")[[1]])
fEMB<-as.numeric(strsplit(fEMB, " ")[[1]])
#-----PLOTTING VENN DIAGRAMS AND EULER-----------------------------------------------------------------#
yS1<- subset.data.frame(y1,y1$family=="S")

YS2<-lapply(yS1,as.numeric)

##continue here
sAAL<- paste0(YS2$AAL[!is.na(YS2$AAL)], collapse = " ")
sCHA<- paste0(YS2$CHA[!is.na(YS2$CHA)], collapse = " ")
sCOL<- paste0(YS2$COL[!is.na(YS2$COL)], collapse = " ")
sZYG<- paste0(YS2$ZYG[!is.na(YS2$ZYG)], collapse = " ")
#ALG<- paste0(YS2$ALG[!is.na(YS2$ALG)], collapse = " ")
sLIV<- paste0(YS2$LIV[!is.na(YS2$LIV)], collapse = " ")
sMOS<- paste0(YS2$MOS[!is.na(YS2$MOS)], collapse = " ")
sHOR<- paste0(YS2$HOR[!is.na(YS2$HOR)], collapse = " ")
sLYC<- paste0(YS2$LYC[!is.na(YS2$LYC)], collapse = " ")
sMON<- paste0(YS2$MON[!is.na(YS2$MON)], collapse = " ")
sGYM<- paste0(YS2$GYM[!is.na(YS2$GYM)], collapse = " ")
sANG<- paste0(YS2$ANG[!is.na(YS2$ANG)], collapse = " ")

sALG<- paste0(YS2$ALG[!is.na(YS2$ALG)], collapse = " ")
sTRA<- paste0(YS2$TRA[!is.na(YS2$TRA)], collapse = " ")
sBRY<- paste0(YS2$BRY[!is.na(YS2$BRY)], collapse = " ")
sEMB<- paste0(YS2$EMB[!is.na(YS2$EMB)], collapse = " ")

sAAL<-as.numeric(strsplit(sAAL, " ")[[1]])
sCHA<-as.numeric(strsplit(sCHA, " ")[[1]])
sCOL<-as.numeric(strsplit(sCOL, " ")[[1]])
sZYG<-as.numeric(strsplit(sZYG, " ")[[1]])
#ALG<-as.numeric(strsplit(ALG, " ")[[1]])
sLIV<-as.numeric(strsplit(sLIV, " ")[[1]])
sMOS<-as.numeric(strsplit(sMOS, " ")[[1]])
sHOR<-as.numeric(strsplit(sHOR, " ")[[1]])
sLYC<-as.numeric(strsplit(sLYC, " ")[[1]])
sMON<-as.numeric(strsplit(sMON, " ")[[1]])
sGYM<-as.numeric(strsplit(sGYM, " ")[[1]])
sANG<-as.numeric(strsplit(sANG, " ")[[1]])

sALG<-as.numeric(strsplit(sALG, " ")[[1]])
sTRA<-as.numeric(strsplit(sTRA, " ")[[1]])
sBRY<-as.numeric(strsplit(sBRY, " ")[[1]])
sEMB<-as.numeric(strsplit(sEMB, " ")[[1]])


#NOT WORKING YET DUE TO UNDEFINED GROUPS! 

##AAL has no introns?!


set1 <- list(ALG,BRY,TRA)
names(set1)<-c("ALG","BRY","TRA")

set2 <- list(ALG,LIV,MOS,HOR)
names(set2)<-c("ALG","LIV","MOS","HOR")

set3 <- list(LIV,MOS,HOR,TRA)
names(set3)<-c("LIV","MOS","HOR","TRA")

set4 <- list(AAL,CHA,CZA,EMB)
names(set4)<-c("AAL","CHA","CZA","EMB")

set5 <- list(ALG,EMB)
names(set5)<-c("ALG","EMB")

fset1 <- list(fLIV,fMOS,fHOR,fLYC)
names(fset1)<-c("fLIV","fMOS","fHOR","fLYC")

fset2 <- list(fALG,fBRY,fTRA)
names(fset2)<-c("fALG","fBRY","fTRA")

sSet1 <- list(sLIV,sMOS,sHOR,sLYC)
names(sSet1)<-c("fLIV","fMOS","fHOR","fLYC")

sSet2 <- list(sALG,sBRY,sTRA)
names(sSet2)<-c("sALG","sBRY","sTRA")

#---------------------------------------------------------------------------------------#


plot(euler(set1, shape = "ellipse"), quantities = TRUE)
plot(euler(fset1, shape = "ellipse"), quantities = TRUE)
plot(euler(fset2, shape = "ellipse"), quantities = TRUE)
plot(euler(sSet1, shape = "ellipse"), quantities = TRUE)
plot(euler(sSet2, shape = "ellipse"), quantities = TRUE)
plot(euler(set2, shape = "ellipse"), quantities = TRUE)
plot(euler(set3, shape = "ellipse"), quantities = TRUE)
plot(euler(set4, shape = "ellipse"), quantities = TRUE)
plot(euler(set5, shape = "ellipse"), quantities = TRUE)

plot(venn(set1))
plot(venn(fset1))
plot(venn(fset2))
plot(venn(sSet2))
plot(venn(set2))
plot(venn(set3))
plot(venn(set4))
#----SELECT DISTINCT FAMILIES OR INTRONS------------------------------------------------------------------#



F2 <- subset.data.frame(xF,xF$family==2)
FT2 <- subset.data.frame(xT2,xT2$family!="S")
ST2 <- subset.data.frame(xT2,xT2$family=="S")
S1 <- subset.data.frame(x3,subset = intron %in% c("cox1i732g2","nad2i709g2","nad7i140g2","trnSi43g2"))
M1 <- subset.data.frame(x3,subset = intron %in% c("cox1i1149g2","nad9i246g2"))

write.table(FT2, file = paste(Filename, "_families_shaded_order.txt"), quote = FALSE, row.names = FALSE, sep = ",")
write.table(xF, file = paste(Filename, "_families_binary.txt"), quote = FALSE, row.names = FALSE, sep = ",")
write.table(FAG, file = paste(Filename, "_families.txt"), quote = FALSE, row.names = FALSE, sep = ",")

xFo<- merge.data.frame(FT2,Forder, by = "intron", all.x = TRUE, all.y = TRUE)
Fo<-xFo[order(xFo$id), ]

#----PLOTTING OPTIONS------------------------------------------------------------------#
medium<-0.5

ggplot(melt(xT2), aes(variable, intron, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  scale_color_gradient2(low = "blue", mid = "white",high = "grey", midpoint = medium) +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(melt(Fo), aes(variable, id, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  scale_color_gradient2(low = "blue", mid = "white",high = "grey", midpoint = medium) +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(x4), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(xF), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(xF2), aes(variable, id, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(xS), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  scale_fill_grey(start = 0, end = .1) +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(F2), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(S1), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(melt(M1), aes(variable, intron, fill = family, alpha = value)) + 
  geom_tile(colour = "white") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  coord_trans(y="reverse") +
  geom_density() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))


