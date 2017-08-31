# example run:
# Rscript task1.R simu_1.gz 100

library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# define the file name
file_name <- strsplit(args[1], "\\.")[[1]][1]

# read data
d <- read.csv(args[1])

# sum umis by barcodes
d1 <- d %>% 
  group_by(column) %>%
  summarise(
    sum = sum(value)
  )
d1 <- as.data.frame(d1)

# order barcodes by the sum of umis
d1 <- d1[order(d1$sum, decreasing = T), ]
write.csv(d1, paste0(file_name, "_barcodes.txt"), quote = F, row.names = F)

# fit linear models to the first N points and save the R-squared
# use the step parameter to differentiate between linear models

# define a step size for linear models
step <- as.numeric(args[2])

d2 <- log10(d1$sum)
tmp <- data.frame(x = log10(1:length(d2)), umi = d2)
rsq <- NULL
for(i in seq(1, nrow(tmp), by = step)) {
  m <- lm(tmp$umi[1:i] ~ tmp$x[1:i])
  rsq <- c(rsq, summary(m)$r.squared)
}

# all barcodes before the R-squared global minimum are real
# others are empty
d1$non_empty_cell <- "n"
d1[1:(which.min(rsq[2:length(rsq)]) * step), ]$non_empty_cell <- "y"
write.csv(d1, paste0(file_name, "_result.txt"), quote = F, row.names = F)
