# f1 <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\Input_data\\mDM_synom_mutat.bed"
# mDM_synom<- read.delim2(f1,header = F)
# colnames(mDM_synom) <- c("CHROM","POS","ID","REF","ALT")
#mDM_synom$label <- c(rep(1,1444),rep(0,nrow(mDM_synom)-1444))
fa <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_syn_driver_mutations.bed"
driver_synom_mDM_sites <- read.delim2(fa,header = F)
fb <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_syn_passenger_mutations.bed"
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


library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicFeatures)

pos_mDM_synom <- mDM_synom[mDM_synom$label==1,]
neg_mDM_synom <- mDM_synom[mDM_synom$label==0,]

select_negsamples <- neg_mDM_synom[sample(1:nrow(neg_mDM_synom),nrow(pos_mDM_synom)),]

######################################
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# 定义函数生成突变前后的序列
generate_dnabert_inputs <- function(mutations_df, k = 101,seq_label_name) {
  # k应为奇数，使突变位点在中心
  if (k %% 2 == 0) {
    k <- k + 1
    message(paste("Adjusted k to", k, "for symmetric flanking"))
  }
  
  half_k <- floor(k / 2)
  mut_pos_in_seq <- half_k + 1  # 突变在序列中的位置
  
  # 创建GenomicRanges对象
  gr <- GRanges(
    seqnames = mutations_df$CHROM,
    ranges = IRanges(start = mutations_df$POS,
                     end = mutations_df$POS),
    strand = "*",
    ref = mutations_df$REF,
    alt = mutations_df$ALT
  )
  
  # 扩展获取侧翼序列
  flanking_gr <- resize(gr, width = k, fix = "center")
  
  # 获取参考序列（野生型）
  wt_sequences <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, flanking_gr))
  
  # 生成突变序列（突变型）
  mt_sequences <- sapply(1:length(wt_sequences), function(i) {
    seq_chars <- strsplit(wt_sequences[i], "")[[1]]
    
    # 验证参考碱基匹配
    ref_base <- mutations_df$REF[i]
    if (seq_chars[mut_pos_in_seq] != ref_base) {
      warning(paste("Mismatch at position", i, 
                    "Expected", ref_base, 
                    "but found", seq_chars[mut_pos_in_seq]))
    }
    
    # 应用突变
    seq_chars[mut_pos_in_seq] <- mutations_df$ALT[i]
    return(paste(seq_chars, collapse = ""))
  })
  
  # 创建结果数据框
  result <- data.frame(
    chr = mutations_df$CHROM,
    position = mutations_df$POS,
    ref = mutations_df$REF,
    alt = mutations_df$ALT,
    wt_sequence = wt_sequences,
    mt_sequence = mt_sequences,
    mutation_position_in_sequence = mut_pos_in_seq,
    sequence_length = k,
    stringsAsFactors = FALSE
  )
  
  # 添加序列ID（可选）
  result$sequence_id <- paste0(seq_label_name,1:nrow(result))
  
  return(result)
}


# 生成序列
pos_inputs <- generate_dnabert_inputs(mutations_df=pos_mDM_synom, k = 251,"pos")
#neg_inputs <- generate_dnabert_inputs(mutations_df=select_negsamples, k = 251,"neg")
neg_inputs <- generate_dnabert_inputs(mutations_df=neg_mDM_synom, k = 251,"neg")
dnabert_inputs <- rbind(pos_inputs,neg_inputs)

#output_dir <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\Input_data\\fast_files\\"
output_dir <- "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\RNA_fold\\Mu_VS_WT\\fast_file\\"
write_fasta <- function(sequences, ids, filename) {
  fasta_lines <- character(length(sequences) * 2)
  for (i in seq_along(sequences)) {
    fasta_lines[(i-1)*2 + 1] <- paste0(">", ids[i])
    fasta_lines[(i-1)*2 + 2] <- sequences[i]
  }
  writeLines(fasta_lines, file.path(output_dir, filename))
}
# 野生型序列
write_fasta(dnabert_inputs$wt_sequence, 
            paste0(dnabert_inputs$sequence_id, "_WT"),
            "wildtype_251sequences.fasta")

# 突变型序列
write_fasta(dnabert_inputs$mt_sequence,
            paste0(dnabert_inputs$sequence_id, "_MT"),
            "mutant_251sequences.fasta")





















###################################################
# get_mutat_seq <- function(mDM_synoms,seq_length){
#   mutat_seq <- DNAStringSet()
#   genome <- BSgenome.Hsapiens.UCSC.hg38
#   target_GR <- GRanges(seqnames = as.character(mDM_synoms$CHROM),
#                        IRanges(start = as.numeric(as.character(mDM_synoms$POS)),
#                                width = 1),ref = mDM_synoms$REF,
#                        alt = mDM_synoms$ALT,
#   )
#   half_k <- floor(seq_length / 2)
#   flanking_gr <- resize(target_GR, width = seq_length, fix = "center")
#   # 标记突变位点位置（在k长度序列中的位置）
#   mut_position <- half_k + 1
#   
#   for (i in 1:length(target_GR)) {
#     starts <-( start(target_GR)[i]-(seq_length-1)/2)
#     ends <- ( start(target_GR)[i]+(seq_length-1)/2)
#     one_seq_GR <- GRanges(seqnames = as.character(seqnames(target_GR)[i]),
#                           IRanges(start = starts,end = ends))
#     
#     one_seq <- getSeq(genome, one_seq_GR)
#     mutat_seq <- c(mutat_seq,one_seq)
#     
#   }
#   return(mutat_seq)
# }
# names(neg_mutatselect_seqs) <- paste0("sites_",(length(pos_mutat_seqs)+1):(length(pos_mutat_seqs)+length(neg_mutatselect_seqs)))
# writeXStringSet(neg_mutatselect_seqs,"D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\neg\\selectneg_101_seqs.fasta")
# 
# 
# pos_mutat_seqs <- get_mutat_seq(mDM_synoms = pos_mDM_synom,seq_length = 101)
# neg_mutat_seqs <- get_mutat_seq(mDM_synoms = neg_mDM_synom,seq_length = 101)
# neg_mutatselect_seqs <- get_mutat_seq(mDM_synoms = select_negsamples,seq_length = 101)
# 
# names(pos_mutat_seqs) <- paste0('sites_',1:length(pos_mutat_seqs))
# names(neg_mutat_seqs) <- paste0("sites_",(length(pos_mutat_seqs)+1):(length(pos_mutat_seqs)+length(neg_mutat_seqs)))
# 
# library(Biostrings)
# writeXStringSet(pos_mutat_seqs,"D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\pos\\pos_101_seqs.fasta")
# writeXStringSet(neg_mutat_seqs,"D:\\research\\m6Acancer_prediction\\m6A_mutation\\mutation_data\\new_data\\driver_passenger\\RNA_fold\\neg\\neg_101_seqs.fasta")
