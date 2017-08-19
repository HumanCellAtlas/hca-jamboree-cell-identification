README
by Will Townes

This folder contains the results of model based clustering only on log2 of total UMI counts of the simulated datasets. Each file represents the results of the fit for a single simulated dataset. The column names are:
"barcode" ID for the "cell", same as original data
"tc" total count (number of UMIs)
"h" entropy of the count distribution
"l2tc" log2 of the total count
"simulation" ID of the simulation, same as the filename
"cl" cluster prediction indicator, the numbers are meaningless except as category labels
"pred_cell" is this barcode predicted to be a real cell (True/False)
