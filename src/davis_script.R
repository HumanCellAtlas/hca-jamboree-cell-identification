library(tidyverse)
simu_1 <- read_csv("jamboree/task1_simu/simu_1.gz")
head(simu_1)
length(unique(simu_1[[2]]))
length(unique(simu_1[[1]]))
colnames(simu_1) <- c("gene", "barcode", "count")
simu_1_stats <- group_by(simu_1, barcode) %>%
    summarise(total_counts = sum(count),
              mean_counts = mean(count),
              var_counts = var(count),
              genes_detected = sum(count > 0),
              prop_dropout = mean(count > 0)
              )

ggplot(simu_1_stats, aes(x = total_counts, genes_detected)) +
    geom_point(alpha = 0.3)

ggplot(simu_1_stats, aes(x = log10(total_counts))) +
    geom_density()
