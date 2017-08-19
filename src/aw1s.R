rm(list = ls())

library(tidyverse)
original_input = read_csv(file = "jamboree/task1_simu/simu_1.gz")
#barcode = cell barcode, count = umi_count
colnames(original_input) = c("gene", "barcode", "count")

NUM_TRUNC = 1e5
input = original_input[1:NUM_TRUNC,]

barcode_df = input %>% group_by(barcode) %>%
  summarise(total_umis = sum(count), total_genes = length(gene))

ggplot(barcode_df, aes(log10(total_genes))) + stat_ecdf(geom = "step")

plot(density(log10(barcode_df$total_umis)))
plot(density(log10(barcode_df$total_genes)))

AMBIENT_LOWER_CUTOFF = 100
AMBIENT_UPPER_CUTOFF = 

really_bad_df = barcode_df %>% filter(total_genes <= 1)

genes_in_really_bad = input %>% filter(barcode %in% really_bad_df$barcode) %>%
  group_by(gene) %>% summarize(freq_perc = length(gene) / nrow(really_bad_df)) %>% arrange(desc(freq_perc))

gene_df = input %>% group_by(gene) %>%
  summarise(total_umis = sum(count), total_barcodes = length(barcode))

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
