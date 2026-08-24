load("D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\synMALL_feat.Rdata")
synMALL_pos <- as.numeric(as.character(select_synMALL_feat$POS))
fa <-  "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\pos_WER_encode.csv"
fb <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\neg_WER_encode.csv"

pos_WERencodes <- read.csv(fa)
new_WERencodes <- read.csv(fb)
WERencodes <- rbind(pos_WERencodes,new_WERencodes)
WER_sum=rowSums(WERencodes[,-c(1:2)])
max_WER <- max(WER_sum)
min_WER <- min(WER_sum)
WERsum_norm <- (WER_sum-min_WER)/(max_WER-min_WER)
WERencodes_sum <- data.frame(WERencodes[,1:2],WER_sum=WERsum_norm)

# f1 <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\processed\\RNAfold_101_feat.csv"
# RNA_fold_feat <- read.csv(f1)

# new_addfeat <- data.frame()
# for (i in 1:nrow(WERencodes)) {
#   one_feat <- cbind(select_synMALL_feat[select_synMALL_feat$POS==synMALL_pos[i],],WERencodes[WERencodes$POS==synMALL_pos[i],-c(1:2)])
#   new_addfeat <- rbind(new_addfeat,one_feat)
# }

#feat_datas <- new_addfeat[,5:ncol(new_addfeat)]
feat_datas <- select_synMALL_feat[,5:ncol(select_synMALL_feat)]
library(dplyr)
library(mice)
char_cols <- feat_datas %>% 
  select(where(is.character)) %>% 
  names()

select_colname <- colnames(feat_datas)[which(is.na(match(colnames(feat_datas),char_cols)))]
feat_data <- feat_datas[,colnames(feat_datas)%in%select_colname]

rm_col <- vector(length = ncol(feat_data))

for (i in 1:ncol(feat_data)) {
  one_feature <- feat_data[,i]
  if(length(which(is.na(one_feature)))/length(one_feature)>=0.7){
    rm_col[i] <- i
  }
  if(length(which(is.na(one_feature)))/length(one_feature)<0.7){
    rm_col[i] <- NA
  }
}
new_feat <- feat_data[,-which(!is.na(rm_col))]

select_col <- vector(length = ncol(new_feat))
for (i in 1:ncol(new_feat)) {
  one_feat <- new_feat[,i]
  if(length(which(is.na(one_feat)))==0){
    select_col[i] <- TRUE
  }
  if(length(which(is.na(one_feat)))!=0){
    select_col[i] <- FALSE
  }
}
full_feat <- new_feat[,select_col]
full_feat <- as.data.frame(apply(full_feat,2,as.numeric))
full_feat_standardized <- full_feat %>%
  mutate(across(where(is.numeric), ~ {
    (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
  }))

load("D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\imputed_NA_value.Rdata")
complete_data <- complete(imputed, 1)
new_complete <- complete_data[,-grep("macie_anyclass",colnames(complete_data))]
filldata_standardized <- new_complete %>%
  mutate(across(where(is.numeric), ~ {
    (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
  }))
##synMALL feat
synmall_feat <- cbind(filldata_standardized,full_feat_standardized)
synmall_feat <- synmall_feat[,-c(grep("m6a_peak",colnames(synmall_feat)))]
#synmall_feat_select <- synmall_feat[,c(1:106, grep("max",colnames(synmall_feat)),133:ncol(synmall_feat))]
synmall_feat_pos <- data.frame(POS=synMALL_pos,synmall_feat)


# new_addfeat <- data.frame()
# for (i in 1:nrow(WERencodes)) {
#   one_feat <- cbind(WERencodes[WERencodes$POS==synMALL_pos[i],-c(1:2)],synmall_feat_pos[synmall_feat_pos$POS==synMALL_pos[i],-1])
#   new_addfeat <- rbind(new_addfeat,one_feat)
# }
# f1 <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\processed\\RNAfold_101_feat.csv"
# RNA_fold_feat <- read.csv(f1)
f1 <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\processed\\diff_MEF.csv"
RNA_fold_feat <- read.csv(f1)
f6 <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\GC_content\\GC.csv"
GC_content <- read.csv(f6)
GC_contents <- GC_content[,c(1:4,grep("GCContent",colnames(GC_content)))]
#GC_contents_rowMean <- data.frame(POS=as.numeric(GC_contents$pos),GC_count=rowMeans(GC_contents[,-c(1:4)]))
GC_contents_rowMean <- data.frame(POS=as.numeric(GC_contents$pos),GC_count=GC_contents[,grep("X500GCContent",colnames(GC_contents))])
###Number of m6A peaks
# fa <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\m6A_Data_REPIC_m6Aatlas_GEO_hg38.bed"
# m6A_peak <- read.delim2(fa,header = F)
# library(GenomicFeatures)
# m6A_GR <- GRanges(
#   seqnames = as.character(m6A_peak$V1),
#   ranges = IRanges(start =as.numeric(as.character(m6A_peak$V2)), 
#                    end = as.numeric(as.character(m6A_peak$V3))))
# 
# synsites_GR <-GRanges(seqnames = paste0("chr",as.character(select_synMALL_feat$CHROM)),
#                     IRanges(start = as.numeric(as.character(select_synMALL_feat$POS)),
#                             width = 1))
# overlap_m6A_num <- vector()
# for (i in 1:length(synsites_GR)) {
#   one_overlap <- findOverlaps(synsites_GR[i],m6A_GR)
#   one_overlap_num <- length(one_overlap@to)
#   overlap_m6A_num <- c(overlap_m6A_num,one_overlap_num)
# 
# }
new_addfeat <- data.frame()
for (i in 1:length(synMALL_pos)) {
  one_feat <- data.frame(WER_sum=WERencodes_sum[WERencodes_sum$POS==synMALL_pos[i],-c(1:2)],
                         MFE=RNA_fold_feat[RNA_fold_feat$POS==synMALL_pos[i],ncol(RNA_fold_feat)],
                         GC_content=GC_contents_rowMean[GC_contents_rowMean$POS==synMALL_pos[i],ncol(GC_contents_rowMean)],

                         synmall_feat_pos[synmall_feat_pos$POS==synMALL_pos[i],-1])
  new_addfeat <- rbind(new_addfeat,one_feat)
}
MFE_scale <- (new_addfeat$MFE-min(new_addfeat$MFE))/(max(new_addfeat$MFE)-min(new_addfeat$MFE))
new_addfeat$MFE <- MFE_scale
GC_scale <- (new_addfeat$GC_content-min(new_addfeat$GC_content))/(max(new_addfeat$GC_content)-min(new_addfeat$GC_content))
new_addfeat$GC_content <- GC_scale
# pos_samples <- synmall_feat[synmall_feat$label==1,]
# neg_samples <- synmall_feat[synmall_feat$label==0,]
pos_samples <- new_addfeat[new_addfeat$label==1,]
neg_samples <- new_addfeat[new_addfeat$label==0,]
all_raw_feat <- rbind(pos_samples[,-ncol(pos_samples)],neg_samples[,-ncol(neg_samples)])
write.csv(all_raw_feat,file = "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\all_raw_feat.csv",row.names = F)
label_df <- data.frame(label=c(rep(1,nrow(pos_samples)),rep(0,nrow(neg_samples))))
write.csv(label_df,file = "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\label_df.csv",row.names = F)
