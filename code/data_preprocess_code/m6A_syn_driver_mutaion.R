#####passenger synonymous mutations
#mDM identification and definition
library(readr)
library(dplyr)
library(GenomicRanges)
library(maftools)
library(tidyverse)
library(stringr)
library(tidyr)
library(data.table)
f1 <- "/home/zhangteng/disk1/zhangteng/Reasearch/big_data/Synonymous_mutation_set.tsv"
mutation_set <- fread(f1)

mutation_set <- mutation_set %>%
  distinct(Gene.name, HGVSG, Sample.name, .keep_all = TRUE)
mutation_set <- mutation_set %>%
  filter(grepl("Substitution", Mutation.Description))

## Recurrence screening 
mutation_set <- mutation_set %>%
  add_count(HGVSG, name = "mutation_count")
mDMs <- mutation_set %>%
  filter(mutation_count >=6 )

## m6A screening
mDMs <- mDMs %>%
  filter(grepl("A>", Mutation.CDS))

mDM <- mDMs %>%separate(Mutation.genome.position, into = c("Seqname", "Start", "End"), sep = "[:-]")
rm_row <- c(14245)
rm_row <- c(74378, 78711, 174978, 207296)
mDM <-mDM[-rm_row,]
mDM$Start <- as.numeric(as.character(mDM$Start))
mDM$End <- as.numeric(as.character(mDM$End)) 
# Convert Seqname, Start, and End columns to GenomicRanges1 objects
GenomicRanges1 <-  GRanges(
  seqnames = mDM$Seqname,
  ranges = IRanges(start = mDM$Start, end = mDM$End)
)
####paper provide m6A sites range mapping to m6A sites
fa <- "/home/zhangteng/Research/m6A_mutation/data/m6A_Data_REPIC_m6Aatlas_GEO_hg38.bed.txt"
m6A_data <- read.delim2(fa,header = F)
colnames(m6A_data) <- c("m6A_Seqname","m6A_Start","m6A_End","width","strand")
m6A_data$m6A_Seqname <- str_remove(m6A_data$m6A_Seqname,"chr")
m6A_range_split <- m6A_data
m6A_start_site <- as.numeric(as.character(m6A_range_split$m6A_Start))
m6A_end_site <- as.numeric(as.character(m6A_range_split$m6A_End))
# Convert m6A_Seqname, m6A_Start, and m6A_End columns to GenomicRanges2 objects
GenomicRanges2 <-  GRanges(
  seqnames = m6A_range_split$m6A_Seqname,
  ranges = IRanges(start = m6A_start_site, end = m6A_end_site)
)
# Find intersection of intervals between GenomicRanges1 and GenomicRanges2
GenomicRanges3 <- GenomicRanges1[GenomicRanges1 %over% GenomicRanges2]

# Convert GenomicRanges3 to "Seqname:Start-End" format
overlap_pos <- with(GenomicRanges3, paste0(seqnames(GenomicRanges3), ":", start(GenomicRanges3), "-", end(GenomicRanges3)))

# Filter rows in mDM where Mutation.genome.position values are present in overlap_pos
mDM_filtered <- mDM[with(mDM, paste0(Seqname, ":", Start, "-", End) %in% overlap_pos), ]

## Reliability screening
# Generate the HGVSG2 column
mDM_filtered$HGVSG2 <- gsub(":g\\.", "-", mDM_filtered$HGVSG)
mDM_filtered$HGVSG2 <- gsub(">", "-", mDM_filtered$HGVSG2)
mDM_filtered$HGVSG2 <- gsub("([[:alpha:]])([[:digit:]])", "\\1-\\2", mDM_filtered$HGVSG2)
mDM_filtered$Mutation.genome.position <- with(mDM_filtered, paste0(Seqname, ":", Start, "-", End))

##load A-to-I dataset
fb <- "/home/zhangteng/Research/m6A_mutation/data/atoi_edit_pos.tsv"
atoi_edit_pos <- fread(fb)
conda <- !(mDM_filtered$Mutation.genome.position %in% atoi_edit_pos$pos)
mDM_filtereda <- mDM_filtered[conda,]
### mapping to WER binding sites
load("/home/zhangteng/Research/m6A_mutation/data/hg38_WER_bindsites.Rdata")
WER_bind_sites_pos <- as.character(WER_bind_sites$pos)
WER_bind_sites_pos <- as.data.frame(WER_bind_sites_pos)
WER_sites_split <- separate(WER_bind_sites_pos,WER_bind_sites_pos, into = c("Seqname", "Start", "End"), sep = "[:-]")
###do the overlap
mDM_filtereda$Start <- as.numeric(as.character(mDM_filtereda$Start))
mDM_filtereda$End <- as.numeric(as.character(mDM_filtereda$End)) 
WER_sites_split$Start <- as.numeric(as.character(WER_sites_split$Start))
WER_sites_split$End <- as.numeric(as.character(WER_sites_split$End))

# Convert Seqname, Start, and End columns to GenomicRanges1 objects
mDM_GR <-  GRanges(
  seqnames = mDM_filtereda$Seqname,
  ranges = IRanges(start = mDM_filtereda$Start, end = mDM_filtereda$End)
)

WER_GR <-  GRanges(
  seqnames = WER_sites_split$Seqname,
  ranges = IRanges(start = WER_sites_split$Start, end = WER_sites_split$End)
)
# Find intersection of intervals between GenomicRanges1 and GenomicRanges2
mDM_overlap_GR <- mDM_GR[mDM_GR %over% WER_GR]
# Convert mDM_overlap_GR to "Seqname:Start-End" format
overlap_pos <- with(mDM_overlap_GR, paste0(seqnames(mDM_overlap_GR), ":", start(mDM_overlap_GR), "-", end(mDM_overlap_GR)))
condb <- (mDM_filtereda$Mutation.genome.position %in% overlap_pos)
mDM_filteredb <- mDM_filtereda[condb,]
##gnomAD filtering
load("/home/disk1/zhangteng/Reasearch/big_data/gnomAD_data_process.Rdata")

match_idx <- which(!is.na(match(new_gnomAD$mutation_iden_last,mDM_filteredb$HGVSG2)))
match_gnomAD <- new_gnomAD[match_idx,]
no_match_idx <- which(is.na(match(mDM_filteredb$HGVSG2,new_gnomAD$mutation_iden_last)))
no_match_mDM <- mDM_filteredb[no_match_idx,]
AF_value <- as.numeric( str_remove(match_gnomAD$INFO,"AF="))
match_gnomAD$AF_value <- AF_value
select_match_gnomAD <- match_gnomAD[match_gnomAD$AF_value<0.05,]
select_match_mDM <- mDM_filteredb[which(!is.na(match(mDM_filteredb$HGVSG2,
                                                    select_match_gnomAD$mutation_iden_last))),]

###match label 
gnomAD_match <- !is.na(match(mDM_filteredb$HGVSG2,select_match_gnomAD$mutation_iden_last))
gnomAD_nomatch<- is.na(match(mDM_filteredb$HGVSG2,new_gnomAD$mutation_iden_last))

#info_vec  <- new_gnomAD$INFO[match_idx]  
cond1     <- gnomAD_nomatch | gnomAD_match
somatic_hgvs <- unique(mDM_filteredb$HGVSG[
  grepl("Confirmed somatic variant", mDM_filteredb$Mutation.somatic.status,
        ignore.case = TRUE)
])
extra_keep <- mDM_filteredb$HGVSG %in% somatic_hgvs
cond <- cond1|extra_keep
mDM_filteredc <- mDM_filteredb[cond,]
save(mDM_filteredc,file="/home/zhangteng/Research/m6A_mutation/data/driver_mDM_synom.Rdata")
#####
library(GenomicRanges)
load("D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\driver_mDM_synom.Rdata")
# Read the hg38 genome annotation file (GTF format)
gtf <- read.delim("D:\\GTF_files\\hg38\\genes.gtf", header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

# Extract the CDS positions mapping to CDS region
cds <- gtf[gtf$V3 == "CDS", c("V1", "V4", "V5")]

# Create GenomicRanges objects
cds_gr <- GRanges(seqnames = str_remove(cds$V1,"chr"),
                  ranges = IRanges(start = cds$V4,end = cds$V5))
GenomicRanges1 <-  GRanges(
  seqnames = mDM_filteredc$Seqname,
  ranges = IRanges(start = mDM_filteredc$Start, end = mDM_filteredc$End)
)

cds_synom_muta <- subsetByOverlaps(GenomicRanges1,cds_gr)
overlap_cds_pos <- with(cds_synom_muta, paste0(seqnames(cds_synom_muta), ":", start(cds_synom_muta), "-", end(cds_synom_muta)))
overlap_cds_driver_synom <- mDM_filteredc[with(mDM_filteredc, paste0(Seqname, ":", Start, "-", End) %in% overlap_cds_pos), ]
####select unique driver mutation
driver_synom_mDM_sites <- unique(overlap_cds_driver_synom$Mutation.genome.position)
select_dri_synom_mDM <- data.frame()
for (i in 1:length(driver_synom_mDM_sites)) {
  one_dri_synom_mDM <- overlap_cds_driver_synom[overlap_cds_driver_synom$Mutation.genome.position==driver_synom_mDM_sites[i],]
  one_dri_synom_mDM <- one_dri_synom_mDM[which.max(one_dri_synom_mDM$Gene.CDS.length),]
  select_dri_synom_mDM <- rbind(select_dri_synom_mDM,one_dri_synom_mDM)
}

gene_name <-  str_remove(as.character(select_dri_synom_mDM$Gene.name), "_ENST.*")
select_dri_synom_mDM$Gene.name <- gene_name
# 或者更简单的版本
ref_alleles <- str_extract(as.character(select_dri_synom_mDM$HGVSG2), "[A-Z](?=-[A-Z]$)")
# 提取替代等位基因 (A, C, T)
alt_alleles <- str_extract(as.character(select_dri_synom_mDM$HGVSG2), "(?<=-)[A-Z]$")
Seqnames <- paste0("chr",select_dri_synom_mDM$Seqname)
pos <- as.numeric(as.character(select_dri_synom_mDM$Start))
drive_synom_mDMbed <-  data.frame(chrom=Seqnames,start=pos,end=pos,REF=ref_alleles,ALT=alt_alleles)
write.table(
  drive_synom_mDMbed,
  file = "D:\\research\\m6Acancer_prediction\\m6A_mutation\\paper_wirting\\data\\Input_data\\m6A_synom_driver_mutations.bed",  # 替换为你希望保存的文件路径和名称
  sep = "\t",              # BED文件是制表符分隔
  col.names = FALSE,       # BED文件通常不包含列名行
  row.names = FALSE,       # 不写入行名
  quote = FALSE            # 不写入引号，避免其他软件读取时出错
)
