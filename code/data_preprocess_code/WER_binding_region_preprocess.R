###load WER binding sties data
fa <- "/home/zhangteng/Research/m6A_mutation/data/hg38_Human_m6A_ClipSeqResult.txt"
WER_association <- fread(fa)
colname_infor <- as.character(WER_association[7,])
colnames(WER_association) <- colname_infor
WER_association<- WER_association[-7,]
WER_association$Seqname <-  str_remove(as.character( WER_association$chr),"chr")
WER_association$pos <-  with(WER_association, paste0(Seqname, ":", as.character(WER_association$position)))
protein_WERassociation <- WER_association[WER_association$gene_type=="protein-coding",]
protein_WERassociation <- protein_WERassociation %>% distinct()
###
WER_pos_table<-as.data.frame(table(protein_WERassociation$position))
##only one WER target
one_WER <- WER_pos_table[WER_pos_table$Freq==1,]
one_WER_asso <- protein_WERassociation[which(!is.na(match(protein_WERassociation$position,one_WER$Var1))),]
oneWER_bind_sites <- as.character(one_WER_asso$position)
oneWER_onebind_site <- one_WER_asso[!grepl(",", oneWER_bind_sites),]
oneWER_multiple_sites <- one_WER_asso[grepl(",", oneWER_bind_sites),]
#library(splitstackshape)
library(tidyr)
oneWER_multiplesites_split <- oneWER_multiple_sites %>% 
  separate_rows(position, sep = ",")
oneWER_multisites_split <- oneWER_multiplesites_split[,-ncol(oneWER_multiplesites_split)]
oneWER_multisites_split$pos <-  with(oneWER_multisites_split, paste0(Seqname, ":", as.character(oneWER_multisites_split$position)))
oneWER_multisites_split <- as.data.frame(oneWER_multisites_split)
oneWER_bind_sites <- rbind(oneWER_onebind_site,oneWER_multisites_split)
###multiple WER target
multi_WER <- WER_pos_table[WER_pos_table$Freq>1,]
multi_WER_asso <- protein_WERassociation[which(!is.na(match(protein_WERassociation$position,multi_WER$Var1))),]
###select only one bind site from multiple WER
multiWER_bind_sites <- as.character(multi_WER_asso$position)
multiWER_onebind_site <- multi_WER_asso[!grepl(",", multiWER_bind_sites),]
multiWER_multiple_sites <- multi_WER_asso[grepl(",", multiWER_bind_sites),]
multiWER_multiplesites_split <- as.data.frame(multiWER_multiple_sites %>% 
                                                separate_rows(position, sep = ","))

multiWER_multisites_split <- multiWER_multiplesites_split[,-ncol(multiWER_multiplesites_split)]
multiWER_multisites_split$pos <-  with(multiWER_multisites_split, paste0(Seqname, ":", as.character(multiWER_multisites_split$position)))

multi_WER_bind_sites <- rbind(multiWER_onebind_site,multiWER_multisites_split)
###combine WER bind sites
WER_bind_sites <- rbind(oneWER_bind_sites,multi_WER_bind_sites)
save(WER_bind_sites,file = "/home/zhangteng/Research/m6A_mutation/data/hg38_WER_bindsites.Rdata")
