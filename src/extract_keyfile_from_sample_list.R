rm(list=ls())

library(dplyr)
library(readr)
library(readxl)

setwd("/fs/ess/PAS2444/jignacio/2024/gs-predictions/")

load("data/merged_pheno_dat.Rdata")
merged_dat <- read.csv("data/2024-07-12T214129phenotype_download.csv")
str(merged_dat)
length(unique(merged_dat$germplasmName))
merged_dat$NAME <- merged_dat$germplasmName

sample_column_name <- "GBS_Name"
keyfile_column_names <- c("FullSampleName","original_sample_name","synonym")
sample_file <- "/fs/ess/PAS2444/25NORGRAINS/IL_GBS_LineNames.csv"
out_keyfile <- "/fs/ess/PAS2444/25NORGRAINS/IL_GBS_keyfile.txt"
keyfile_input <- c("/fs/ess/PAS2444/jignacio/2024/gs-predictions/data",
                   "/fs/ess/PAS2444/jignacio/2024/gs-predictions/data/24Historical_Compiled_12K_production_lines_keyfile_05162024.txt",
                   "/fs/ess/PAS2444/25NORGRAINS/25NOR01-09_LIBPREPID_03-26-2025.txt",
                   "/fs/ess/PAS2444/24NORGRAINS_BIG6/HasEmmer_UpdateLines/KeyFile/24NorGrainHasEmmer_23_24Big6_keyfile_041724.txt")

input=keyfile_input[3]
read_input_keyfiles <- function(input){
  if(file.exists(input)){
    if(file_test("-f", input)){
      kf <- read.table(input, header=T, sep="\t", fill=NA, quote="")
      kf$Flowcell <- as.character(kf$Flowcell)
      kf$source_kf_file <- input
    }else if(file_test("-d", input)){
      kf_files_list <- list.files(input,"*keyfile*.txt",full.names=TRUE)
      kf_list <- lapply(kf_files_list,read.table,header=T, sep="\t", fill=NA, quote="")
      kf_list <- lapply(kf_files_list,function(x){
        df <- read.table(x, header=T, sep="\t", fill=NA, quote="")
        df$Flowcell <- as.character(df$Flowcell)
        df$source_kf_file <- x
        return(df)
      })
      kf <- do.call(bind_rows,kf_list) %>%
        distinct(Flowcell, Lane, Barcode, .keep_all = T)
    }
  }else{
    cat("Input file/directory", input, "does not exist. Skipping.")
  }
  return(kf)
}

in_kf <- lapply(keyfile_input, read_input_keyfiles)
all_kf <- do.call(bind_rows, in_kf) %>%
  distinct(Flowcell, Lane, Barcode, .keep_all = T)
unique(all_kf$source_kf_file)

if(grep("*.csv$",sample_file)){
  df1 <- read_csv(sample_file)
}else{
  df1 <- read_table(sample_file)
}
# kf1 <- read.table("data/24Historical_Compiled_12K_production_lines_keyfile_05162024.txt", header=T, sep="\t", fill=NA, quote="")
# str(kf1)
# kf_files_list <- list.files("data","keyfile.*.txt",full.names=TRUE)
# kf_list <- lapply(kf_files_list,read.table,header=T, sep="\t", fill=NA, quote="")
# kf_list <- lapply(kf_files_list,function(x){
#   df <- read.table(x, header=T, sep="\t", fill=NA, quote="")
#   df$source_kf_file <- x
#   return(df)
# })
# lapply(kf_list, dim)
# all_kf <- do.call(bind_rows,kf_list)
# dim(all_kf)

# l1 <- grepl("NY", all_kf$FullSampleName)
# l2 <- grepl("NY", all_kf$FullSampleName)
# ny.idx <- which(grepl("NY",all_kf))
# colnames(all_kf)[ny.idx]
# 
# get_idx <- function(x){
#   grepl("NY", all_kf[,x], ignore.case = T)
# }
# 
# lout <- lapply(ny.idx, get_idx)
# rout <- Reduce("|", lout)
# 
# ny_samples <- all_kf %>%
#   filter(rout)
# 
# write.table(ny_samples, file="output/ny_samples_keyfile_10082024.txt", row.names = F, quote=F, sep="\t")
# 
# ny_samples2 <- ny_samples %>%
#   distinct(FullSampleName, original_sample_name)
# 
# write.table(ny_samples2, file="output/ny_samples_10082024.txt", row.names = F, quote=F, sep="\t")

name="-14933"
# Function to normalize sample names
normalize_name <- function(name) {
  name <- tryCatch({
    name <- gsub("[_\\\\/ -]", "-", name) # Remove spaces, dashes, underscores, backslashes, forward slashes
    name <- gsub("^IL2020-","IL20-", name)
    name <- gsub("^2020-","IL20-", name)
    name <- gsub("^([0-9]{2})-","IL\\1-",name)
    name <- gsub("^PIONEER","", name, ignore.case = T)
    name <- gsub("^Pio","", name, ignore.case = T)
    name <- gsub("^SYN","SY", name, ignore.case = T)
    name <- gsub("^SY-","SY", name, ignore.case = T)
    name <- gsub("^agrimaxx-","agrimaxx", name, ignore.case = T)
    name <- gsub("^-","", name)
    name <- tolower(name)
    return(name)
  }, error = function(e) e, finally = return(name))
  return(name)
}

sample_column_name_idx <- ifelse(!(length(sample_column_name) && nzchar(sample_column_name)),grep(sample_column_name,colnames(df1)),1)
keyfile_column_names_idx <- which(colnames(all_kf) %in% keyfile_column_names)

# Sample vectors
sample_names <- unique(pull(df1, sample_column_name_idx))
normalized_sample_names <- sapply(sample_names, normalize_name)
normalized_kf_names <- lapply(keyfile_column_names_idx, function(idx, df){
  normalized_vector <- sapply(df[,idx], normalize_name)
}, all_kf)

# vector2 <- kf1$FullSampleName
# vector3 <- kf1$original_sample_name

# # Normalize both vectors
# normalized_vector1 <- sapply(vector1, normalize_name)
# normalized_vector2 <- sapply(vector2, normalize_name)
# normalized_vector3 <- sapply(vector3, normalize_name)

# Find matches and their indices in the second vector
name = "kaskaskia"
matches <- sapply(normalized_sample_names, function(name) {
  match_names_list <- lapply(normalized_kf_names, function(kf_name){
    kf_name == name
  })
  match_idx <- which(Reduce("|", match_names_list))
  names(match_idx) <- NULL
  if (length(match_idx) == 0) {
    return(NA)
  } else {
    return(match_idx)
  }
})

head(matches)
sum(unlist(sapply(matches,is.na)))
which(unlist(lapply(sapply(matches,is.na),any)))
ul <- unlist(lapply(matches, length))
names(ul)[which(ul>1)]
sum(table(unlist(lapply(matches, length))))

# grep("4505$",normalized_kf_names[[1]], value = T)
# 
# vector2 <- kf_list[[3]]
# vector1 <- c("SY 100")

new_kf <- all_kf[unlist(matches),]
new_kf$OldFullSampleName <- new_kf$FullSampleName
new_kf$FullSampleName <- new_kf$query_names <- names(unlist(matches))

write_delim(new_kf, out_keyfile, "\t")
# match_names <- function(vector2, vector1) {
#   # Normalize both vectors
#   normalized_vector1 <- sapply(vector1, normalize_name)
#   normalized_vector2 <- sapply(vector2$FullSampleName, normalize_name)
#   normalized_vector3 <- sapply(vector2$original_sample_name, normalize_name)
#   
#   # Find matches and their indices in the second vector
#   matches1 <- match(normalized_vector1, normalized_vector2)
#   names(matches1) <- names(normalized_vector1)
#   
#   matches2 <- match(normalized_vector1, normalized_vector3)
#   names(matches2) <- names(normalized_vector1)
#   
#   # Return results
#   coalesce(matches1, matches2)
#   
# }
# 
# match_names <- function(vector2, vector1) {
#   # Normalize both vectors
#   normalized_vector1 <- sapply(vector1, normalize_name)
#   normalized_vector2 <- sapply(vector2$FullSampleName, normalize_name)
#   normalized_vector3 <- sapply(vector2$original_sample_name, normalize_name)
#   
#   # Find matches and their indices in the second vector
#   matches1 <- match(normalized_vector1, normalized_vector2)
#   names(matches1) <- names(normalized_vector1)
#   
#   matches2 <- match(normalized_vector1, normalized_vector3)
#   names(matches2) <- names(normalized_vector1)
#   
#   # Return results
#   coalesce(matches1, matches2)
#   
# }
# 
# # file_list <- list.files(path = "data", pattern = "_samples.txt", full.names=T)
# # df_list <- lapply(file_list, read.table)
# # lapply(df_list,dim)
# sn <- unique(merged_dat$NAME)
# #sn <- gsub("^OH23-15","15", sn)
# 
# 
# #df_list <- lapply(kf_list,function(x) x$FullSampleName)
# lout <- lapply(kf_list, match_names, sn)
# length(sn)
# dim(merged_dat)
# 
# # lout2 <- lapply(df_list, match_names, merged_dat$NAME)
# # lsubstr <- lapply(lout2, is.na)
# # lsdc <- do.call(bind_cols, lsubstr)
# # merged_dat[rowSums(lsdc)==0,]
# #grep("156",df_list[[1]],value=T)
# #grep("156-19",kf_list[[2]][,"FullSampleName"])
# #grep("156-19",kf_list[[3]][,"FullSampleName"])
# 
# bc <- as.data.frame(do.call(bind_cols, lout))
# dim(bc)
# rownames(bc) <- unique(merged_dat$NAME)
# 
# bc$n.marker.dat <- rowSums(!is.na(bc)) 
# sum(bc$n.marker.dat == 0)
# bc[bc$n.marker.dat == 0,]
# (missing_samples <- rownames(bc[bc$n.marker.dat == 0,]))
# x=trimmed_names[1]
# trimmed_names <- sapply(missing_samples, function(x){
#   len <- nchar(x)
#   #y <- substr(x, 1, 7)
#   y <- substr(x, len-3, len)
#   normalize_name(y)
# })
# all_kf$normFullSampleName <- sapply(all_kf$FullSampleName, normalize_name)
# all_kf$norm_original_sample_name <- sapply(all_kf$original_sample_name, normalize_name)
# sapply(trimmed_names, function(x){
#   grep(x, all_kf$normFullSampleName, value=TRUE, ignore.case = T)
# })
# sapply(trimmed_names, function(x){
#   grep(x, all_kf$norm_original_sample_name, value=TRUE, ignore.case = T)
# })
# sum(bc$n.marker.dat > 0)
# sum(bc$n.marker.dat > 1)
# # View(bc[bc$n.marker.dat == 0,])
# i=1
# sample_list_out <- kf_list[[1]][F,]
# for(i in 1:length(kf_list)){
#   sample_list_tmp <- kf_list[[i]][bc[,i],]
#   sample_list_tmp$queryName <- row.names(bc)
#   infile <- kf_files_list[i]
#   sample_list_tmp$KeyFile2 <- infile
#   sample_list_out <- bind_rows(sample_list_out,sample_list_tmp[!is.na(sample_list_tmp[,"FullSampleName"]),])
#   # outfile <- paste0("output/",gsub("data/(.*).txt","\\1",infile),"_samples.txt")
#   # write.table(sample_list_out, file=outfile, col.names =F, row.names = F, quote=F)
# }
# 
# outfile <- "output/2024-07-12T214129phenotype_download_keyfile_09122024.txt"
# str(sample_list_out)
# sample_list_out2 <- sample_list_out[!with(sample_list_out, duplicated(paste(Flowcell, Lane, Barcode))),]  %>%
#   mutate(KeyFile = gsub("^data\\/","",KeyFile2))
# dim(sample_list_out2)
# length(unique(sample_list_out$FullSampleName))
# length(unique(sample_list_out$original_sample_name))
# write.table(sample_list_out2, file=outfile, row.names = F, quote=F, sep="\t")
# 
# cleaned_keyfile <- read_excel("/fs/ess/PAS2444/jignacio/2024/gs-predictions/data/24Historical_Big6_24_Norgrains_OH_lines_trn_and_pred_keyfile_07092024.xlsx")
# with(cleaned_keyfile, anyDuplicated(paste(Flowcell, Lane, Barcode)))
# 
# cleaned_keyfile[cleaned_keyfile$FullSampleName == "OH18-46-89","Source"]
# # cleaned_outfile <- left_join(cleaned_keyfile, sample_list_out, by=join_by(FullSampleName, plate_id, well), suffix=c("",".remove"), multiple="first") %>%
# #   select(-grep("remove$",colnames(.))) %>%
# #   mutate(KeyFile = gsub("^data\\/","",KeyFile2)) %>%
# #   select(-KeyFile2)
# # write.table(cleaned_outfile,"output/24Historical_Big6_24_Norgrains_OH_lines_trn_and_pred_cleaned_keyfile_07092024.txt",
# # row.names = F, quote=F, sep="\t")