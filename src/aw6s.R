rm(list = ls())
library(tidyverse)

load(file="shared_scratch/group1/task3/regev_loaded.Rda")

#let's focus on batch2 --> simpler
anno = anno %>% filter(batch == "batch2")
regev_df = regev_df[, colnames(regev_df) %in% union(anno$barcode, "symbol")]

regev_mat = data.matrix(regev_df[,-1])
rownames(regev_mat) = regev_df$symbol
transcripts_per_cell = colSums(regev_mat)
stopifnot(all(transcripts_per_cell > 0))
hist(log10(transcripts_per_cell), breaks = 100)

umi_per_cell = colSums(regev_mat)
hist(umi_per_cell, breaks = 1000)
THRESH_TOO_LARGE_CELLS = 40000
too_large = sum(umi_per_cell[umi_per_cell >= THRESH_TOO_LARGE_CELLS]) / sum(umi_per_cell)
cat(sprintf("%.2f%% of the umis were taken by the too large cells", 100*too_large))

regev_mat = regev_mat[, umi_per_cell > THRESH_TOO_LARGE_CELLS]
anno = anno[umi_per_cell > THRESH_TOO_LARGE_CELLS, ]

regev_mat_cpm = 1e6 * sweep(regev_mat, 2, colSums(regev_mat), "/")
regev_mat_log_cpm = log1p(regev_mat_cpm)


# write.csv(file = "shared_scratch/group1/task3/subset_regev_log_cpm.csv", regev_mat_log_cpm)
# write.csv(file = "shared_scratch/group1/task3/subset_regev_anno.csv", anno)

#remove genes that were detected in at most 5 cells
badGenes = rowSums(regev_mat_cpm > 0) <= 5

testMats = list("T0" = regev_mat_log_cpm[!badGenes, anno$timepoint == "T0"],
                "T4" = regev_mat_log_cpm[!badGenes, anno$timepoint == "T4"])

res = list()
for(time_label in names(testMats)) {
  relevantExpMat = testMats[[time_label]] #regev_mat_log_cpm[!badGenes,]      #[, anno$timepoint == "T4"]
  
  #plot(log(mean),log(var/mean))
  
  m = rowMeans(relevantExpMat)
  s = apply(relevantExpMat, 1, sd)
  plot(s ~ m)
  
  aaa = m < 6
  m = m[aaa]
  s = s[aaa]
  
  lo = loess(s ~ m)
  PREDICT_SE = TRUE #FALSE --> to save time
  lo_pred = predict(lo, se = PREDICT_SE) #see https://stackoverflow.com/questions/22717930/how-to-get-the-confidence-intervals-for-lowess-fit-using-r
  if(PREDICT_SE) {
    ci_up = lo_pred$fit + 3 * 1.96 * lo_pred$se.fit
    ci_lo = lo_pred$fit - 3 * 1.96 * lo_pred$se.fit
  }
  #plot(s/m ~ m)
  #plot(s ~ m)
  plot(s ~ m)
  ord = order(m)
  
  if(PREDICT_SE) {
    lines(m[ord], lo_pred$fit[ord], lwd=3, col="red")
    lines(m[ord], ci_lo[ord], lty=2, lwd=3, col="purple")
    lines(m[ord], ci_up[ord], lty=2, lwd=3, col="purple")
    
    res[[time_label]]$lo_var = s <= ci_lo
    res[[time_label]]$hi_var = s >= ci_up
  } else {
    lines(m[ord], lo_pred[ord], lwd=3, col="red")
    res[[time_label]]$lo_var = s <= 0.75 * lo_pred
    res[[time_label]]$hi_var = s >= 1.25 * lo_pred
  }

}
intersect(names(which(res[["T0"]]$lo_var)), names(which(res[["T4"]]$hi_var)))

keep = (s/m) > ci_up #| (s/m) < ci_lo stop("need to keep only highly variable...")
relevantExpMat = relevantExpMat[keep,]

res_pca <- prcomp(t(relevantExpMat), center = TRUE, scale = TRUE)
res_pca$variances = 100 * (res_pca$sdev^2) / sum( (res_pca$sdev)^2 )
names(res_pca$variances) = colnames(res_pca$rotation)
res_pca$scaled_rotation = (scale((res_pca$rotation)))
#sanity
#(scale(t(relevantExpMat), center = TRUE, scale = TRUE)) %*% res_pca$rotation - res_pca$x
animalsToInclude = TRUE
x_axis = "PC1"
y_axis = "PC2"
g1 = ggplot(, aes(x = res_pca$x[, x_axis], y = res_pca$x[, y_axis])) +
  geom_point(size=4, aes(color = anno$MD_SPLIT_GENOTYPE[animalsToInclude])) +
  #ggtitle(matName) + 
  #ggtitle("PCA (all gene space)") + scale_color_brewer(type = "qual", palette = "Paired") +
  #scale_color_manual(values = STIMULUS_COLORS) +
  scale_shape_discrete(solid = TRUE) + 
  xlab(sprintf("%s [%.d%% variance]", x_axis, round(res_pca$variances[x_axis]))) +
  ylab(sprintf("%s [%.d%% variance]", y_axis, round(res_pca$variances[y_axis]))) +
  guides(color = guide_legend(title = "genotype"))
#guides(color = guide_legend(title = "stimulus"), shape = guide_legend(title = "genotype"))
if(TRUE || REWRITE_OUTPUT) print(g1)





SUSPECTED_EMPTY_CUTOFF = 10^3.5
ambiance_mat = regev_mat[, transcripts_per_cell < SUSPECTED_EMPTY_CUTOFF]
regev_mat = regev_mat[, transcripts_per_cell >= SUSPECTED_EMPTY_CUTOFF]


s = apply(mat_10x_1, 1, sd)
badRows = s < 1e-1


set.seed(101)
relevantExpMat <- as.matrix(mat_10x_1[, sample(x = 1:ncol(mat_10x_1), size = 10000)])


# relevantExpMat = mat_10x_1

m = rowMeans(relevantExpMat)

relevantExpMat = relevantExpMat[!badRows,]
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
     
