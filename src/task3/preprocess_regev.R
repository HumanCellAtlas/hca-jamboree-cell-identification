rm(list = ls())
library(tidyverse)

#author: AW

regev_df = read_table2("jamboree/task_3/counts/dead_cells_regev.tsv", col_names = FALSE, skip = 1)
regev_colnames = read_table2("jamboree/task_3/counts/dead_cells_regev.tsv", col_names = FALSE, n_max = 1)
colnames(regev_df) = c("symbol", as.vector(t(regev_colnames[1,])))
regev_df$symbol = sapply(regev_df$symbol, toupper)

tmp = strsplit(colnames(regev_df)[-1], "_")
#quick and dirty for the cases of too few columns
tmp2 = lapply(tmp, function(element) if(length(element) ==2) {return(c(element, "1"))} else {return(element)}) 
#ifelse(length(element)==2, c(element, "1"), element)) #ifelse doesn't work for some reason

stopifnot(all(sapply(tmp2, function(x) length(x) == 3)))
#which(sapply(tmp2, function(x) length(x) != 3)) #-->debug
anno = data.frame(matrix(unlist(tmp2), ncol=3, byrow=T))
colnames(anno) = list("barcode", "timepoint", "replicate")
colnames(regev_df) = make.unique(c("symbol", as.character(anno$barcode)))
anno$barcode = colnames(regev_df)[-1]

#according to Kharcehnko's explanation:
anno$batch = "N/A"
anno$batch[anno$timepoint == "T0" & anno$replicate %in% c("1", "2")] = "batch1"
anno$batch[anno$timepoint %in% c("T0.5", "T1", "T2")] = "batch1"
#batch 2 was also sequenced to higher depths
anno$batch[anno$timepoint == "T0" & anno$replicate == "3"] = "batch2"
anno$batch[anno$timepoint == "T4"] = "batch2"
stopifnot(!any(anno$batch == "N/A"))

save(file="shared_scratch/group1/task3/regev_loaded.Rda", regev_df, anno)
