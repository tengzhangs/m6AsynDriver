# f1<-"/home/disk1/zhangteng/Reasearch/big_data/synMALL/tabix_chr1.tsv"
# library(data.table)
# chr1_tabix <- fread(f1)
# colnames(chr1_tabix)[1] <- "CHROM"
# select_chr1_tabix <- chr1_tabix[,22:ncol(chr1_tabix)]
# select_chr1_tabix2 <- chr1_tabix[,1:4]
# new_chr1_tabix <- cbind(select_chr1_tabix2,select_chr1_tabix)
library(data.table)
library(stringr)
#mDM_synom$label <- c(rep(1,1444),rep(0,nrow(mDM_synom)-1444))
fa <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_syn_drivermutations.bed"
driver_synom_mDM_sites <- read.delim2(fa,header = F)
fb <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\_mDM_syn_passenger_mutations_right.bed"
passenger_synom_mDM_sites <- read.delim2(fb,header = F)
consist_pos <- intersect(driver_synom_mDM_sites$V2,passenger_synom_mDM_sites$V2)
driver_synom_mDM <- driver_synom_mDM_sites[which(is.na(match(driver_synom_mDM_sites$V2,consist_pos))),]
passenger_synom_mDM <- passenger_synom_mDM_sites[which(is.na(match(passenger_synom_mDM_sites$V2,consist_pos))),]
mDM_synom<-rbind(driver_synom_mDM,passenger_synom_mDM)
pos_index <- match(driver_synom_mDM$V2,mDM_synom$V2)
neg_index <- match(passenger_synom_mDM$V2,mDM_synom$V2)
labels_value <- vector(length = nrow(mDM_synom))
labels_value[pos_index] <- 1
labels_value[neg_index] <- 0
mDM_synom$label <- labels_value

colnames(mDM_synom)[1:5] <- c("CHROM","POS","ID","REF","ALT")
synMALL_files <- paste0("tabix_",c( paste0("chr",1:22),"chrX"),".tsv")
file_dir <- "/home/zhangteng/disk1/zhangteng/Reasearch/big_data/synMALL/"
select_synMALL_feat <- data.frame()
for (i in 1:length(synMALL_files)) {
  file_one <- paste0(file_dir,synMALL_files[i])
  one_tabix <- fread(file_one)
  colnames(one_tabix)[1] <- "CHROM"
  select_one_tabix <- one_tabix[,22:ncol(one_tabix)]
  select_one_tabix2 <- one_tabix[,1:4]
  new_one_tabix <- cbind(select_one_tabix2,select_one_tabix)
  one_select_mDM_synom <- mDM_synom[mDM_synom$CHROM==str_extract(synMALL_files[i], "chr[0-9X]+"),]
  common_pos <- intersect(new_one_tabix$POS,one_select_mDM_synom$POS)
  one_mDM_synom <- one_select_mDM_synom[which(!is.na(match(one_select_mDM_synom$POS,common_pos))),]
  select_new_one_tabix <- new_one_tabix[which(!is.na(match(new_one_tabix$POS,common_pos))),]
  last_select_one <- data.frame()
  last_mDM_synom <- data.frame()
  for (j in 1:length(common_pos)) {
    one_pos_mDM <- one_mDM_synom[one_mDM_synom$POS==common_pos[j],]
    last_mDM_synom <- rbind(last_mDM_synom,one_pos_mDM)
    one_pos_tabix <- select_new_one_tabix[select_new_one_tabix$POS==common_pos[j],]
    if(nrow(one_pos_tabix)==1){
      last_select_one <- rbind(last_select_one,one_pos_tabix)
    }
    if(nrow(one_pos_tabix)>1){
      onepos_tabix <- one_pos_tabix[1,]
      last_select_one <- rbind(last_select_one,onepos_tabix)
    }
    
  }
  #select_new_one <- select_new_one_tabix[which(!is.na(select_new_one_tabix$POS)),]
  last_select_one$label <- as.numeric(as.character(last_mDM_synom$label))
  select_synMALL_feat <- rbind(select_synMALL_feat,last_select_one)
}
save(select_synMALL_feat,file = "/home/zhangteng/Research/m6A_mutation/data/driver_passenger/synMALL_data/synMALL_feat.Rdata")
