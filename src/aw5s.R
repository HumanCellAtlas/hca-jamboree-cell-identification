rm(list = ls())
library(tidyverse)

library(ggthemes)
library(pryr)
library(tidyverse)
library(scater)
library(idcells)
library(matrixStats)

df_macosko <- readRDS("/home/jovyan/shared_scratch/group1/Macosko/filtered_data.rds")
head(df_macosko)

df2sparse_macosko <- function(d) {
  stopifnot(all(c("barcode","umi","alignment", "chromosome", "count") %in% colnames(d)))
  d = d %>% dplyr::filter(!(alignment %in% c("none-G", "INTERGENIC", "UTR", "INTRONIC"))) %>%
    group_by(barcode, alignment) %>% summarize(num_uniq_umis = length(unique(umi)))
  g <- as.factor(d$alignment)
  b <- as.factor(d$barcode)
  m <- Matrix::sparseMatrix(i=as.numeric(g), j=as.numeric(b), x=d$count)
  colnames(m) <- levels(b)
  rownames(m) <- levels(g)
  m
}

mat_simu_1 <- df2sparse_macosko(df_macosko)
dim(mat_simu_1)
set.seed(101)
mat_simu_1_small <- mat_simu_1[, sample(x = 1:ncol(mat_simu_1), size = 10000)]


pd <- as(data.frame(barcode = colnames(mat_simu_1_small)), 
         "AnnotatedDataFrame")
rownames(pd) <- colnames(mat_simu_1_small)
fd <- as(data.frame(gene_id = rownames(mat_simu_1_small)), 
         "AnnotatedDataFrame")
rownames(fd) <- rownames(mat_simu_1_small)
sce <- newSCESet(countData = as.matrix(mat_simu_1_small), phenoData = pd,
                 featureData = fd)
ambient_bool <- read_tsv("ambient_genes.txt", col_names = FALSE)
ambient_sparse <- read_tsv("ambient_genes.sparse.txt", col_names = FALSE)
ambient_gene <- read_tsv("ambient_genes.gene_names.txt", col_names = FALSE)
ambient <- bind_cols(ambient_gene, ambient_bool, ambient_sparse)
table(ambient[[1]] %in% featureNames(sce))
sce <- sce[(featureNames(sce) %in% ambient[[1]]),]
ambient <- ambient[ambient[[1]] %in% featureNames(sce),]
fData(sce)$ambient_gene_l2 <- ambient[[2]]
fData(sce)$ambient_gene_l1 <- ambient[[3]]
sce <- sce[rowSums(counts(sce)) > 0.5,]
sce







sce = readRDS("shared_scratch/group1/sceset_simu_1.rds")
relevantExpMat = assayData(sce)$counts

m = rowMeans(relevantExpMat)
s = apply(relevantExpMat, 1, sd)
y_for_fit = s#s/m
plot(y_for_fit ~ m)

lo = loess(y_for_fit ~ m)
lo_pred = predict(lo, se = TRUE) #see https://stackoverflow.com/questions/22717930/how-to-get-the-confidence-intervals-for-lowess-fit-using-r
ci_up = lo_pred$fit + 3 * 1.96 * lo_pred$se.fit
ci_lo = lo_pred$fit - 3 * 1.96 * lo_pred$se.fit
#plot(s/m ~ m)
#plot(s ~ m)
plot(y_for_fit ~ m)
ord = order(m)
lines(m[ord], lo_pred$fit[ord], lwd=3, col="red")
lines(m[ord], ci_lo[ord], lty=2, lwd=3, col="purple")
lines(m[ord], ci_up[ord], lty=2, lwd=3, col="purple")
keep = (s/m) > ci_up #| (s/m) < ci_lo stop("need to keep only highly variable...")
relevantExpMat = relevantExpMat[keep,]


RELOAD_DATA = FALSE
if(RELOAD_DATA) {
  original_input = read_tsv(file = "jamboree/macosko_2015/counts/GSM1629192_Pure_HumanMouse.raw", col_names = FALSE)
  #barcode = cell barcode, count = umi_count
  colnames(original_input) = c("barcode", "umi", "alignment", "chromosome", "count")
  
  
  barcode_df = original_input %>% group_by(barcode) %>%
    summarise(total_transcripts = sum(count), num_unique_umis = length(umi))
  
  #remove barcodes that have less than 5 unique UMIs associated with them
  UNIQ_UMI_THRESH = 5
  bad_barcodes = barcode_df %>% filter(num_unique_umis >= UNIQ_UMI_THRESH) %>% 
    select(barcode)
  
  filtered_input = original_input %>% filter(!(barcode %in% bad_barcodes$barcode))
  
  saveRDS(original_input, file="shared_scratch/group1/Macosko/raw_data.rds")
  saveRDS(filtered_input, file="shared_scratch/group1/Macosko/filtered_data.rds")
} else if(FALSE) {
  start.time = Sys.time()
  filtered_input = readRDS(file="shared_scratch/group1/Macosko/filtered_data.rds")
  no_none_input = filtered_input %>% filter(!(alignment %in% c("none-G", "INTERGENIC", "UTR", "INTRONIC")))
  saveRDS(no_none_input, file="shared_scratch/group1/Macosko/no_none_input.rds")
  end.time = Sys.time()
  time.taken = end.time - start.time
} else {
  no_none_input = readRDS(file="shared_scratch/group1/Macosko/no_none_input.rds")
}


NUM_TRUNC = 1e5
input = no_none_input
# input = original_input[1:NUM_TRUNC,]

barcode_df = input %>% group_by(barcode) %>%
  summarise(total_transcripts = sum(count), num_unique_umis = length(umi))

#ggplot(barcode_df, aes(log10(total_transcripts))) + stat_ecdf(geom = "step")
ggplot(barcode_df, aes(log10(num_unique_umis))) + stat_ecdf(geom = "step")

plot(density(log10(barcode_df$num_unique_umis)))
# plot(density(log10(barcode_df$total_genes)))

AMBIENT_LOWER_CUTOFF = 10^0#10^0.25#100
AMBIENT_UPPER_CUTOFF = 10^0.35#10^0.35#2^10

HIGH_QUALITY_LOWER_CUTOFF = 10^4

ambient_df = barcode_df %>% filter(num_unique_umis >= AMBIENT_LOWER_CUTOFF & num_unique_umis <= AMBIENT_UPPER_CUTOFF)
ggplot(ambient_df, aes(log10(num_unique_umis))) + stat_ecdf(geom = "step")
plot(density(log2(1+ambient_df$num_unique_umis)))

# high_quality_df = barcode_df %>% filter(num_unique_umis >= HIGH_QUALITY_LOWER_CUTOFF)
# ggplot(high_quality_df, aes(log10(num_unique_umis))) + stat_ecdf(geom = "step")
# plot(density(log2(1+high_quality_df$num_unique_umis)))


genes_in_ambient_df = input %>% filter(barcode %in% ambient_df$barcode) %>%
  group_by(alignment) %>% summarize(freq = length(alignment) / nrow(ambient_df)) %>% arrange(desc(freq))

hist((genes_in_ambient_df %>% filter(freq >= 0))$freq, breaks = 100)
AMBIENT_GENES_NAIVE_CUTOFF = 0#0.1#0.4
genes_in_ambient_df = genes_in_ambient_df %>% filter(freq >= AMBIENT_GENES_NAIVE_CUTOFF)

#genes_in_high_quality_df = input %>% filter(barcode %in% high_quality_df$barcode) %>%
#  group_by(gene) %>% summarize(freq = length(gene) / nrow(high_quality_df)) %>% arrange(desc(freq))


library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
att = listAttributes(ensembl) %>% filter(page == "feature_page") %>% dplyr::select(-page)
grep("symbol", att$description, ignore.case = TRUE, value = TRUE)
biomart_res1 = getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), 
                     filters = 'ensembl_gene_id', 
                     values = genes_in_ambient_df$alignment, 
                     mart = ensembl)

genes_in_ambient_df$hgnc_symbol = sapply(genes_in_ambient_df$alignment, toupper)
#biomart_res1[match(sapply(genes_in_ambient_df$alignment, toupper),
#     biomart_res1[,1]), 
#      2]

write.csv(file="shared_scratch/group1/Macosko/ambient_genes.csv", genes_in_ambient_df)

biomart_res2 = getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), 
                     filters = 'ensembl_gene_id', 
                     values = genes_in_high_quality_df$gene, 
                     mart = ensembl)

genes_in_high_quality_df$hgnc_symbol =  biomart_res2[match(sapply(genes_in_high_quality_df$gene, toupper),
                                                           biomart_res2[,1]), 
                                                     2]

genesInAmbientSimulationAndrew = read.delim("shared_scratch/group1/ambient_genes.stringent.gene_names.txt", header = FALSE)
temp = read.delim("shared_scratch/group1/ambient_genes.stringent.sparse.txt", header = FALSE)
genesInAmbientSimulationAndrew = cbind(genesInAmbientSimulationAndrew, temp)
colnames(genesInAmbientSimulationAndrew) = c("ensemble_id", "are_ambient")

biomart_res1 = getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), 
                     filters = 'ensembl_gene_id', 
                     values = genesInAmbientSimulationAndrew %>% filter(are_ambient) %>% dplyr::select(ensemble_id), 
                     mart = ensembl)
write.csv(file="shared_scratch/group1/Macosko/aaa.csv", biomart_res1)


stop("until here)
     
     # total_umis = gene_df_raw %>% summarise_at(vars(count), funs(total_umis = sum))
     # total_barcodes = gene_df_raw %>% summarise_at(vars(barcode), funs(total_barcodes = n()))
     
     plot(density(log10(1+gene_df$total_umis)))
     plot(density(log10(1+gene_df$total_barcodes)))
     
     
     
     hist(log10(1+gene_df$total_umis), breaks = 100)
     hist(log10(1+gene_df$total_barcodes), breaks = 100)
     
     
     transcripts_per_cell = input %>% group_by(barcode) %>% 
     summarise_at(vars(count), funs(sum = sum))
     
     #how do I do that in dplyr? with n() doesn't work
     count_cells = table(input$barcode)
     