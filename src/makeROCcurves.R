
# Take as input (1) truth files from simulation and (2) files with algorithm scores
# and create ROC plots for each simulated dataset

setwd("/home/jovyan/scratch/")

library(data.table)
library(dplyr)
library(ROCR)

# (( repeat this for each SIM number ))
for (SIM in 1:16){
  
# read in output file
output <-  read.csv(paste0("~/shared_scratch/group1/numpy_arrays/GOOD_CELL_SCORES/simu_",
                           SIM, "_scores_revised.txt"), stringsAsFactors=FALSE,
                           header=FALSE)

# read in truth file
truth <- fread(paste0("~/shared_scratch/group1/truth/ssimu/real_", SIM), 
               header=TRUE, stringsAsFactors = FALSE)

# First check error rate in cells we threw away (don't even have predictions on)
# want to make sure AUC is 1 for these
perc.in <- round((sum(truth$barc %in% output$V1) / length(truth$barc))*100,3)
message(paste0(perc.in, 
               " % of the true cells were scored"))

# add the truth for the cells we have predictions on
x <- match(truth$barc, output$V1)
x <- x[!is.na(x)]
output$truth <- 0
output$truth[x] <- 1

# rank on total UMI counts
raw <- fread(paste0("simu_", SIM))
umis <- dplyr::group_by(raw, column) %>% dplyr::summarise(UMI=sum(value))
rm(raw)
gc()

# subset UMIs to include only the ones we scored
x <- match(output$V1, umis$column)
umis <- umis[x,]

# Roc curve
pred.obj <- prediction(output$V2, 1-output$truth)
perf.obj <- performance(pred.obj, 'tpr','fpr')
AUC.obj <- performance(pred.obj,"auc")@y.values[[1]]

pred.naive <- prediction(-umis$UMI, 1-output$truth)
perf.naive <- performance(pred.naive, 'tpr', 'fpr')
AUC.naive <- performance(pred.naive,"auc")@y.values[[1]]

#plot in pdf format with a grey line with theoretical random results
pdf(paste0("/home/jovyan/shared_scratch/group2/ROC/ROC_Simulation_", SIM, ".pdf"))
  plot(perf.obj, main=paste0("Simulation ", SIM, " ROC Curve (", perc.in, "% scored)"),
       col="purple", lwd=1)
  plot(perf.naive, col="blue", add=TRUE)
  legend(0.5, 0.4, legend=c(paste0("Wolball (AUC=", round(AUC.obj,3), ")"), 
                            paste0("Total UMIs (AUC=", round(AUC.naive,3), ")")), 
         col=c("purple", "blue"), lty=c(1,1))
  abline(0,1,col="grey", lty=2)
  
  plot(perf.obj, main=paste0("Simulation ", SIM, " ROC Curve (", perc.in, "% scored)"),
       col="purple",
       xlim=c(0,0.2))
  plot(perf.naive, col="blue", add=TRUE)
  legend(0.1, 0.4, legend=c(paste0("Wolball (AUC=", round(AUC.obj,3), ")"), 
                             paste0("Total UMIs (AUC=", round(AUC.naive,3), ")")),
         col=c("purple", "blue"), lty=c(1,1))
  abline(0,1,col="grey", lty=2)
  
  # same thing but just for top 20,000 ranked by UMI count
  top20K <- which(rank(-umis$UMI, ties="random") <= 20000)
  output <- output[top20K,]
  umis <- umis[top20K,]
  
  perc.in <- round((sum(truth$barc %in% umis$column)/ length(truth$barc))*100, 3)
  message(paste0(perc.in, 
                 " % of the true cells were in the top 20K by total UMI"))
  
  pred.obj <- prediction(output$V2, 1-output$truth)
  perf.obj <- performance(pred.obj, 'tpr','fpr')
  AUC.obj <- performance(pred.obj,"auc")@y.values[[1]]
  
  pred.naive <- prediction(-umis$UMI, 1-output$truth)
  perf.naive <- performance(pred.naive, 'tpr', 'fpr')
  AUC.naive <- performance(pred.naive,"auc")@y.values[[1]]
  
  plot(perf.obj, main=paste0("Simulation ", SIM, " - top 20K cells (",
                             perc.in, "% in top 20K by UMI total)"),
       col="purple", lwd=1)
  plot(perf.naive, col="blue", add=TRUE)
  legend(0.5, 0.4, legend=c(paste0("Wolball (AUC=", round(AUC.obj,3), ")"), 
                            paste0("Total UMIs (AUC=", round(AUC.naive,3), ")")),
         col=c("purple", "blue"), lty=c(1,1))
  abline(0,1,col="grey", lty=2)
  
dev.off()
}