rm(list = ls())

library(tidyverse)
original_input = readRDS("shared_scratch/group1/task1_simu.filtered/simu_1.rds")

NUM_TRUNC = 1e5
input = original_input
#input = original_input[1:NUM_TRUNC,]

barcode_df = input %>% group_by(barcode) %>%
  summarise(total_umis = sum(count), total_genes = length(gene))

ggplot(barcode_df, aes(log10(total_umis))) + stat_ecdf(geom = "step")

plot(density(log2(1+barcode_df$total_umis)))
# plot(density(log10(barcode_df$total_genes)))

AMBIENT_LOWER_CUTOFF = 100
AMBIENT_UPPER_CUTOFF = 2^10

HIGH_QUALITY_LOWER_CUTOFF = 10^4

ambient_df = barcode_df %>% dplyr::filter(total_umis >= AMBIENT_LOWER_CUTOFF & total_umis <= AMBIENT_UPPER_CUTOFF)
ggplot(ambient_df, aes(log10(total_umis))) + stat_ecdf(geom = "step")
plot(density(log2(1+ambient_df$total_umis)))

high_quality_df = barcode_df %>% dplyr::filter(total_umis >= HIGH_QUALITY_LOWER_CUTOFF)
ggplot(high_quality_df, aes(log10(total_umis))) + stat_ecdf(geom = "step")
plot(density(log2(1+high_quality_df$total_umis)))


genes_in_ambient_df = input %>% dplyr::filter(barcode %in% ambient_df$barcode) %>%
  group_by(gene) %>% summarize(freq = length(gene) / nrow(ambient_df)) %>% dplyr::arrange(desc(freq))

hist((genes_in_ambient_df %>% dplyr::filter(freq >= 0.1))$freq, breaks = 100)
AMBIENT_GENES_NAIVE_CUTOFF = 0.1#0.4
genes_in_ambient_df = genes_in_ambient_df %>% dplyr::filter(freq >= AMBIENT_GENES_NAIVE_CUTOFF)

genes_in_high_quality_df = input %>% dplyr::filter(barcode %in% high_quality_df$barcode) %>%
  group_by(gene) %>% summarize(freq = length(gene) / nrow(high_quality_df)) %>% dplyr::arrange(desc(freq))


library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
att = listAttributes(ensembl) %>% dplyr::filter(page == "feature_page") %>% dplyr::select(-page)
grep("symbol", att$description, ignore.case = TRUE, value = TRUE)
biomart_res1 = getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), 
      filters = 'ensembl_gene_id', 
      values = genes_in_ambient_df$gene, 
      mart = ensembl)

genes_in_ambient_df$hgnc_symbol =  biomart_res1[match(sapply(genes_in_ambient_df$gene, toupper),
                                        biomart_res1[,1]), 
                                        2]

write.csv(file = "shared_scratch/group1/aw_simu1_ambient_genes.csv",  genes_in_ambient_df)



biomart_res2 = getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), 
                    filters = 'ensembl_gene_id', 
                    values = genes_in_high_quality_df$gene, 
                    mart = ensembl)

genes_in_high_quality_df$hgnc_symbol =  biomart_res2[match(sapply(genes_in_high_quality_df$gene, toupper),
                                                     biomart_res2[,1]), 
                                               2]
write.csv(file = "shared_scratch/group1/aw_high_quality_genes.csv",  genes_in_high_quality_df)

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
