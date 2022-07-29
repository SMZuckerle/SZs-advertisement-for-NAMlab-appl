#---JCSC-------------------------------------------------------------------#
#---JACUSA CROSS VALIDATION AND SUBSETTING Tool for R----------------------#
#---Simon Zumkeller 2020---------------------------------------------------#
#--------------------------------------------------------------------------#

#---------SETTINGS---------------------------------------------------------#
############################################################################
#---Set your working directory where your INPUT and OUTPUT files are/------#
#---will be located--------------------------------------------------------#

setwd(
  "G:/Machina_Eva/R_home/Cross_ValidatoR")

#---Enter the INPUT filename into the space between hyphens----------------#
#--------------------------------------------------------------------------#

Filename= "Ecoli-offtargets-19-datasets-VK20201209.csv"

#---OUTPUT----------------------------------------------------------------#
#---One table should be generated-----------------------------------------#
#-------------------------------------------------------------------------#
###########################################################################
#-------------------------------------------------------------------------#

library(dplyr)
library(data.table)
library(tidyverse)

x<-read.csv(
  Filename
  , sep=";", header=TRUE, stringsAsFactors=TRUE)
x$strand <-
  ifelse(x$REF == "A"| x$REF =="G", "1","0")
x$DATASET <-
  gsub('_','-',x$DATASET)
colnames(x)[1] <- "POS"

y<- x %>%
  select (1,2)

z<- x %>%
  select (1,3,16)
z<-z[!duplicated(z),]

Read_aggr_DNA <- as.data.table(
  x)[, sum(DNA), by = .(POS)]
colnames(
  Read_aggr_DNA)[2] <- "DNA"
Read_aggr_A_DNA <- as.data.table(
  x)[, sum(A_DNA), by = .(POS)]
colnames(
  Read_aggr_A_DNA)[2] <- "A_DNA"
Read_aggr_C_DNA <- as.data.table(
  x)[, sum(C_DNA), by = .(POS)]
colnames(
  Read_aggr_C_DNA)[2] <- "C_DNA"
Read_aggr_G_DNA <- as.data.table(
  x)[, sum(G_DNA), by = .(POS)]
colnames(
  Read_aggr_G_DNA)[2] <- "G_DNA"
Read_aggr_T_DNA <- as.data.table(
  x)[, sum(T_DNA), by = .(POS)]
colnames(
  Read_aggr_T_DNA)[2] <- "T_DNA"
Read_aggr_RNA <- as.data.table(
  x)[, sum(RNA), by = .(POS)]
colnames(
  Read_aggr_RNA)[2] <- "RNA"
Read_aggr_A_RNA <- as.data.table(
  x)[, sum(A_RNA), by = .(POS)]
colnames(
  Read_aggr_A_RNA)[2] <- "A_RNA"
Read_aggr_C_RNA <- as.data.table(
  x)[, sum(C_RNA), by = .(POS)]
colnames(
  Read_aggr_C_RNA)[2] <- "C_RNA"
Read_aggr_G_RNA <- as.data.table(
  x)[, sum(G_RNA), by = .(POS)]
colnames(
  Read_aggr_G_RNA)[2] <- "G_RNA"
Read_aggr_T_RNA <- as.data.table(
  x)[, sum(T_RNA), by = .(POS)]
colnames(
  Read_aggr_T_RNA)[2] <- "T_RNA"

Read_merge<- merge(
  merge(
    merge(
      merge(
        merge(
          merge(
            merge(
              merge(
                merge(
    Read_aggr_DNA,Read_aggr_A_DNA, by = "POS"),
  Read_aggr_C_DNA, by = "POS"),
  Read_aggr_G_DNA, by = "POS"),
  Read_aggr_T_DNA, by = "POS"),
  Read_aggr_RNA, by = "POS"),
  Read_aggr_A_RNA, by = "POS"),
  Read_aggr_C_RNA, by = "POS"),
  Read_aggr_G_RNA, by = "POS"),
  Read_aggr_T_RNA, by = "POS")

y <- y[,c(2,1)]

colnames(
  y)[1] <- "frame"
colnames(
  y)[2] <- "POS"

split_list_y<-split.data.frame(
  y,y$frame)

split_list_y[1]
split_list_y[2]

y_names <- names(split_list_y)
y_names2 <- as.data.frame(
  names(split_list_y))

for(
  i in seq_along(y_names)
){
  assign(
    y_names[i], split_list_y[[i]])
}

common_df_all <- merge(
  merge(
    merge(
      merge(
        merge(
          merge(
            merge(
              merge(
                merge(
                  merge(
                    merge(
                      merge(
                        merge(
                          merge(
                            merge(
                              merge(
                                merge(
                                  merge(
                                    merge(
  `0-WT`,
  `1-PPR65-1`,by = "POS", all=TRUE),
  `2-PPR65-2`,by = "POS", all=TRUE),
  `3-PPR56-3`,by = "POS", all=TRUE),
  `4-PPR56-4`,by = "POS", all=TRUE),
  `5-PPR56-notarget`,by = "POS", all=TRUE),
  `6-PPR56-nad3`,by = "POS", all=TRUE),
  `7-PPR56-nad3-nad4`,by = "POS", all=TRUE),
  `8-PPR65-S4NN`,by = "POS", all=TRUE),
  `9-PPR65-P3TN`,by = "POS", all=TRUE),
  `10-PPR78-N54P-C4`,by = "POS", all=TRUE),
  `11-PPR78-N54P-R200`,by = "POS", all=TRUE),
  `12-PPR98`,by = "POS", all=TRUE),
  `13-NODE1476`,by = "POS", all=TRUE),
  `14-PPR56-NOD2497`,by = "POS", all=TRUE),
  `15-PPR65-P9TD`,by = "POS", all=TRUE),
  `16-PPR56-PPR45`,by = "POS", all=TRUE),
  `17-PPR65-P6`,by = "POS", all=TRUE),
  `18-NODE2497-atpA`,by = "POS", all=TRUE),
  `19-PPR65-P3-P6`,by = "POS", all=TRUE)

colnames(
  common_df_all)[2] <- "0_WT"
colnames(
  common_df_all)[3] <- "1_PPR65_1"
colnames(
  common_df_all)[4] <- "2_PPR65_2"
colnames(
  common_df_all)[5] <- "3_PPR56_3"
colnames(
  common_df_all)[6] <- "4_PPR56_4"
colnames(
  common_df_all)[7] <- "5_PPR56_notarget"
colnames(
  common_df_all)[8] <- "6_PPR56_nad3"
colnames(
  common_df_all)[9] <- "7_PPR56_nad3_nad4"
colnames(
  common_df_all)[10] <- "8_PPR65_S4NN"
colnames(
  common_df_all)[11] <- "9_PPR65_P3TN"
colnames(
  common_df_all)[12] <- "10_PPR78_N54P_C4"
colnames(
  common_df_all)[13] <- "11_PPR78_N54P_R200"
colnames(
  common_df_all)[14] <- "12_PPR98"
colnames(
  common_df_all)[15] <- "13_NODE1476"
colnames(
  common_df_all)[16] <- "14_PPR56-NOD2497"
colnames(
  common_df_all)[17] <- "15_PPR65_P9TD"
colnames(
  common_df_all)[18] <- "16_PPR56-PPR45"
colnames(
  common_df_all)[19] <- "17_PPR65_P6"
colnames(
  common_df_all)[20] <- "18_NODE2497_atpA"
colnames(
  common_df_all)[21] <- "19_PPR65_P3_P6"

common_df_all<-
  merge.data.frame(
    common_df_all,z, by = "POS", all.x = TRUE, all.y = TRUE)
common_df_all2<-
  merge.data.frame(
    common_df_all,Read_merge, by = "POS", all.x = TRUE, all.y = TRUE)

write.csv(
  common_df_all2,file=paste(
    Filename, "_JCSC.csv"), quote = FALSE)
####################################################################

####################################################################



##############################################################DONE##
