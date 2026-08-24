fa <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\processed\\mutant_251sequences.csv"
muta_MEF <- read.csv(fa)
###
fb <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\processed\\wildtype_251sequences.csv"
wt_MEF <- read.csv(fb)

diff_MEF <- ( muta_MEF$MFE)-(wt_MEF$MFE)
label <- c(rep(1,1439),rep(0,(nrow(muta_MEF)-1439)))
diff_MEF_infor <- data.frame(diff_MEF=diff_MEF,label=label)

fa <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_syn_driver_mutations.bed"
driver_synom_mDM_sites <- read.delim2(fa,header = F)
fb <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_syn_passenger_mutations_right.bed"
passenger_synom_mDM_sites <- read.delim2(fb,header = F)
consist_pos <- intersect(driver_synom_mDM_sites$V2,passenger_synom_mDM_sites$V2)
driver_synom_mDM <- driver_synom_mDM_sites[which(is.na(match(driver_synom_mDM_sites$V2,consist_pos))),]
passenger_synom_mDM <- passenger_synom_mDM_sites[which(is.na(match(passenger_synom_mDM_sites$V2,consist_pos))),]
mDM_synom <- rbind(driver_synom_mDM,passenger_synom_mDM)
colnames(mDM_synom) <- c("CHROM","POS","ID","REF","ALT")
pos_index <- match(driver_synom_mDM$V2,mDM_synom$POS)
neg_index <- match(passenger_synom_mDM$V2,mDM_synom$POS)
labels_value <- vector(length = nrow(mDM_synom))
labels_value[pos_index] <- 1
labels_value[neg_index] <- 0
mDM_synom$label <- labels_value
RNA_fold_feat <- data.frame(mDM_synom,diff_MEF=diff_MEF)

write.csv(RNA_fold_feat,file ="D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\processed\\diff_MEF.csv",row.names = F )
