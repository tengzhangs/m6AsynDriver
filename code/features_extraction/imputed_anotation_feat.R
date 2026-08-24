load("D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\synMALL_feat.Rdata")
feat_datas <- select_synMALL_feat[,5:ncol(select_synMALL_feat)]
library(dplyr)

char_cols <- feat_datas %>% 
  select(where(is.character)) %>% 
  names()

feat_data <- feat_datas[,which(is.na(match(colnames(feat_datas),char_cols)))]

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
###
no_full_feat <- new_feat[,!select_col]
# fill_feat <- data.frame()
# for (i in 1:ncol(no_full_feat)) {
#   one_nofull <- no_full_feat[,i]
#   one_nofull_value <- one_nofull[!is.na(one_nofull)]
#   
# }
library(mice)
new_nofull_feat <- no_full_feat[,-c(grep("erprspval",colnames(no_full_feat)),grep("gnomad",colnames(no_full_feat)))]
imputed <- mice(new_nofull_feat, seed = 123,print = FALSE)
save(imputed,file = "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\imputed_NA_value.Rdata")
