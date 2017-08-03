---
title: Task 1 README, Group 4
author: Aaron Lun, Samantha Riesenfeld, Tallulah Andrews, The Phuong Dao, Tomas Gomes
date: 2 August 2017
---

# Algorithm overview

We assumed that empty droplets would contain transcripts that were randomly sampled from the ambient RNA pool.
We obtained an estimate of transcript proportions in the ambient pool, by summing UMI counts for each gene across all low-count barcodes (i.e., all barcodes with the lowest 85% of total UMI counts).
We applied the Good-Turing algorithm to the counts to obtain a posterior estimate of the proportion for each gene.
This ensures that all proportions were non-zero, avoiding situations with log-likelihoods of negative infinity.

For each barcode, we tested whether its expression profile was significantly different to the ambient pool.
We computed the expected count for each gene by scaling the ambient proportions by the total UMI count for that barcode.
The count for each gene was considered to be sampled from a Poisson distribution with mean equal to the expected count.
We computed the sum of log-likelihoods across genes for each barcode, and subtracted it from the log-likelihood for a fully saturated model to obtain a likelihood ratio (LR).

We derived a p-value corresponding to the LR of each cell using a simulation approach.
For a simulated barcode, we randomly sampled Poisson-distributed counts using the ambient proportions and a given total UMI count.
We performed this for a range of total counts, and fitted a trend to the simulated LRs against the total UMI count.
The variability of simulated LRs around the trend was modelled with a log-normal distribution with total count-dependent parameters.
The upper tail of this distribution was used to obtain a p-value for each observed LR of each barcode.

We applied the Benjamini-Hochberg correction to adjust for multiple testing across barcodes.
All barcodes with fewer than 250 UMIs were not considered during correction (and were assigned adjusted p-values of unity).
This reduces the burden of the multiple testing correction for barcodes corresponding to real cells.
All barcodes detected at a FDR of 5% were called as real cells. 

# File manifest

- `model_fitter.R`, a script that computes adjusted p-values for all cells.
- `sparse_maker.R`, a script that generates serialized sparse matrices for use in model fitting.
