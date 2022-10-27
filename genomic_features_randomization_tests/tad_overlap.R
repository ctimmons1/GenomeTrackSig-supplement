# load packages
library(tidyverse)
library(GenomeTrackSig)
devtools::load_all()

# import helpful functions not in package
source("../helper_functions.R")

# get chromosome lengths
txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- base::suppressMessages(GenomicFeatures::genes(txdb))
chrom_lengths <- stats::setNames(object = all.genes@seqinfo@seqlengths[1:24], all.genes@seqinfo@seqnames[1:24])
chrom_lengths_mb_bins <- floor(chrom_lengths/1e6)+1

# Find the proportion of changepoints in a given sample that
# overlap with a TAD boundary
generateSampleOverlaps <- function(traj, tads) {
  # get changepoints
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)

  starts <- changepoints$mb_bin_start
  ends <- changepoints$mb_bin_end
  lower <- sapply(starts, function(x){which.min(abs(tads$mb_bin_start - x))})
  upper <- sapply(ends, function(x){which.min(abs(tads$mb_bin_end - x))})
  tf <- abs(starts - tads$mb_bin_end[lower]) > abs(ends - tads$mb_bin_start[upper])
  closest <- replace(lower, which(tf==TRUE), upper[which(tf==TRUE)])

  cp_range <- mapply(function(x, y){IRanges::IRanges(x:y)}, starts, ends)
  tad_range <- sapply(closest, function(x){rep(tads$mb_bin_end[x]:tads$mb_bin_after[x], each=1)})
  tad_range <- sapply(tad_range, function(x){IRanges::IRanges(x)})
  overlaps <- mapply(function(x, y){IRanges::overlapsAny(x, y)}, cp_range, tad_range)
  overlaps <- sapply(overlaps, function(x){if (TRUE %in% x) {x = 1} else {x = 0}})
  overlap_prop <- length(overlaps[overlaps>0])/length(overlaps)

  return (overlap_prop)
}

# Generate a null distribution of 10,000 samples with randomly-placed changepoints
# Find proportion of random changepoints that overlap with TAD boundary
# Report proportion
generateRandomOverlaps <- function(traj, blank_genome, tads) {
  # get changepoints -- to know number of cps in the sample
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)

  n = 100

  binData <- traj[['binData']]
  if (!is.null(binData$actual_bin)) {
    binData <- binData[order(nrow(binData):1),]
  }
  binData$mb_bin_start <- 1
  binData$mb_bin_end <- 1

  for (i in 1:nrow(binData)) {
    if (binData$end[i] %% 1e6 == 0) {
      if (binData$end[i] != (binData$start[i+1]-1)) {
        binData$start[i+1] = binData$end[i]+1
      }
    }
    else {
      if (binData$end_chrom[i] == 24 & binData$start_chrom[i] == 23) {
        if (binData$end[i] != 59373566) {
          binData$end[i] <- binData$start[i+1] - 1
        }
      }
    }
    if (i == nrow(binData)) {
      if (binData$end_chrom[i] == 24) {
        binData$end[i] <- 59373566
        }
    }
    mb_bin_start <- blank_genome %>%
      dplyr::filter(start_chrom == binData$start_chrom[i]) %>%
      dplyr::filter(start == binData$start[i])
    mb_bin_end <- blank_genome %>%
      dplyr::filter(end_chrom == binData$end_chrom[i]) %>%
      dplyr::filter(end == binData$end[i])

    binData$mb_bin_start[i] = mb_bin_start$mb_bin
    binData$mb_bin_end[i] = mb_bin_end$mb_bin

  }

  indices <- t(replicate(n, sample(1:nrow(binData), nrow(changepoints), replace=TRUE)))
  ends <- binData$mb_bin_end[indices]
  starts <- binData$mb_bin_start[indices]
  lower <- sapply(starts, function(x){which.min(abs(tads$mb_bin_end - x))})
  upper <- sapply(ends, function(x){which.min(abs(tads$mb_bin_start - x))})
  tf <- abs(starts - tads$mb_bin_end[lower]) > abs(ends - tads$mb_bin_start[upper])
  closest <- replace(lower, which(tf==TRUE), upper[which(tf==TRUE)])

  cp_range <- mapply(function(x, y){IRanges::IRanges(x:y)}, starts, ends)
  tad_range <- sapply(closest, function(x){tads$mb_bin_end[x]:tads$mb_bin_after[x]})
  tad_range <- sapply(tad_range, function(x){IRanges::IRanges(x)})
  overlaps <- mapply(function(x, y){IRanges::overlapsAny(x, y)}, cp_range, tad_range)
  overlaps <- sapply(overlaps, function(x){if (TRUE %in% x) {x = 1} else {x = 0}})
  overlaps <- matrix(unlist(overlaps), ncol = nrow(changepoints), nrow =n)
  overlap_props <- apply(overlaps, 1, function(x){length(x[x>0])/length(x)})

  return (overlap_props)
}

# calculate the fraction of null proportions that are greater than sample proportion
testOverlapSignificance <- function(random_overlaps, sample_overlap) {
  p <- length(random_overlaps[random_overlaps>=sample_overlap])/length(random_overlaps)
  return (p)
}

# perform overlap test on a single sample
overlapTest <- function(traj, blank_genome, tads) {
  random_overlaps <- generateRandomOverlaps(traj, blank_genome, tads)
  sample_overlaps <- generateSampleOverlaps(traj, tads)
  pval <- testOverlapSignificance(random_overlaps, sample_overlaps)
  return (pval)
}

# Automates testing
# Input folder containing signature profiles for a cancer type and template genome binned by Mb
# Performs overlap test on changepoints in original sample
# Generates null distribution of 10,000 random samples and
# the proportion of their changepoints that overlap with a TAD boundary
# Compares sample proportion to null proportion
# Reports the fraction of null proportions that are lower than sample proportion
testAllSamples <- function(type, blank_genome, tads) {
  pvals_all <- c()
  ids <- c()


  for (i in c(1:330)) {
    if (file.exists(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))) {
        load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))
        temp <- load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))
        # check for changepoints
        if (!is.null(get(temp)[[2]])||!is.null(get(temp)[[6]])||!is.null(get(temp)[[10]])||!is.null(get(temp)[[14]])||!is.null(get(temp)[[18]])) {
            ids <- c(ids, i)
            pval <- overlapTest(get(temp), blank_genome, tads)
            pvals_all <- c(pvals_all, pval)
        }
    }
  }

  results <- data.frame('id' = ids, 'pval' = pvals_all)
  return (results)
}





# make dataframe of information about all changepoints in a trajectory
summarizeChangepoints <- function(trajectory, sampleID) {
  cpPos <- assignChangepoints(trajectory,0)
  if (nrow(cpPos) > 0) {
    cpPos$start_chrom <- 1
    cpPos$end_chrom <- 1
    cpPos$before_chrom <- 1
    cpPos$after_chrom <- 1
    cpPos$start <- 1
    cpPos$end <- 1
    cpPos$mb_bin_start <- 1
    cpPos$mb_bin_end <- 1
    cpPos$start_before <- 1
    cpPos$end_after <- 1
    cpPos$mb_bin_before <- 1
    cpPos$mb_bin_after <- 1
    cpPos$sampleID <- sampleID
    if (nrow(cpPos) == 1) {
      if (!is.null(trajectory[['binData']]$actual_bin)) {
        cpPos$start_before[1] <- 1
        cpPos$before_chrom[1] <- 1
        cpPos$end_after[1] <- 59373566
        cpPos$after_chrom[1] <- 24
        cpPos$start[1] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==cpPos$cpPos1[1]]
        cpPos$end[1] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==cpPos$cpPos2[1]]
        cpPos$start_chrom[1] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[1]]
        cpPos$end_chrom[1] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos2[1]]
        cpPos <- cpPos %>%
          dplyr::mutate(start_before = if_else(before_chrom != start_chrom, 1, start_before),
                        end_after = if_else(after_chrom != end_chrom, max(blank_genome$end[blank_genome$start_chrom==end_chrom]), end_after),
                        before_chrom = if_else(before_chrom != start_chrom, start_chrom, before_chrom),
                        after_chrom = if_else(after_chrom != end_chrom, end_chrom, after_chrom))
      }
      else {
        cpPos$start_before[1] <- 1
        cpPos$before_chrom[1] <- 1
        cpPos$end_after[1] <- 59373566
        cpPos$after_chrom[1] <- 24
        cpPos$start[1] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==cpPos$cpPos1[1]]
        cpPos$end[1] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==cpPos$cpPos2[1]]
        cpPos$start_chrom[1] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[1]]
        cpPos$end_chrom[1] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==cpPos$cpPos2[1]]
        cpPos <- cpPos %>%
          dplyr::mutate(start_before = if_else(before_chrom != start_chrom, 1, start_before),
                        end_after = if_else(after_chrom != end_chrom, max(blank_genome$end[blank_genome$start_chrom==end_chrom]), end_after),
                        before_chrom = if_else(before_chrom != start_chrom, start_chrom, before_chrom),
                        after_chrom = if_else(after_chrom != end_chrom, end_chrom, after_chrom))
      }
      mb_bin_start <- blank_genome %>%
        dplyr::filter(start_chrom == cpPos$start_chrom[1] & start == cpPos$start[1])
      cpPos$mb_bin_start[1] <- mb_bin_start$mb_bin
      mb_bin_end <- blank_genome %>%
        dplyr::filter(end_chrom == cpPos$end_chrom[1] & end == cpPos$end[1])
      cpPos$mb_bin_end[1] <- mb_bin_end$mb_bin
    }
    else {
      cpPos <- cpPos %>%
        dplyr::arrange(start_chrom, start)
      for (i in 1:nrow(cpPos)) {
        if (!is.null(trajectory[['binData']]$actual_bin)) {
          if (i == 1) {
            cpPos$start_before[i] <- 1
            cpPos$end_after[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])-1]
            cpPos$before_chrom[i] <- 1
            cpPos$after_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])-1]

          }
          else if (i == nrow(cpPos)) {
            cpPos$start_before[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
            cpPos$end_after[i] <- 59373566
            cpPos$before_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
            cpPos$after_chrom[i] <- 24
          }
          else {

            if (cpPos$mb_bin_start[i] > cpPos$mb_bin_start[i-1]) {
              cpPos$start_before[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
              cpPos$end_after[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])-1]
              cpPos$before_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
              cpPos$after_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])-1]

            }
            else {
              cpPos$start_before[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
              cpPos$end_after[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])]
              cpPos$before_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos2[i-1])]
              cpPos$after_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==(cpPos$cpPos1[i+1])-1]

            }

          }
          cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
          cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==cpPos$cpPos2[i]]
          cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
          cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos2[i]]


        }
        else {
          cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
          cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==cpPos$cpPos2[i]]
          cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
          cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==cpPos$cpPos2[i]]

          if (i == 1) {
            cpPos$start_before[i] <- 1
            cpPos$end_after[i] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==(cpPos$cpPos1[i+1]-1)]
            cpPos$before_chrom[i] <- 1
            cpPos$after_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==(cpPos$cpPos1[i+1]-1)]

          }
          else if (i == nrow(cpPos)) {
            cpPos$start_before[i] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==(cpPos$cpPos2[i-1])]
            cpPos$end_after[i] <- 59373566
            cpPos$before_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==(cpPos$cpPos2[i-1])]
            cpPos$after_chrom[i] <- 24

          }
          else {
            cpPos$start_before[i] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==(cpPos$cpPos2[i-1])]
            cpPos$end_after[i] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==(cpPos$cpPos1[i+1]-1)]
            cpPos$before_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==(cpPos$cpPos2[i-1])]
            cpPos$after_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==(cpPos$cpPos1[i+1]-1)]

          }
        }
        mb_bin_start <- blank_genome %>%
          dplyr::filter(start_chrom == cpPos$start_chrom[i] & start == cpPos$start[i])
        cpPos$mb_bin_start[i] <- mb_bin_start$mb_bin
        mb_bin_end <- blank_genome %>%
          dplyr::filter(end_chrom == cpPos$end_chrom[i] & end == cpPos$end[i])
        cpPos$mb_bin_end[i] <- mb_bin_end$mb_bin
      }
    }
    return (cpPos)
  }
  else {
    return (NULL)
  }
}

main <- function() {
    devtools::load_all()
    args <- commandArgs(trailingOnly=TRUE)
    type <- as.character(args[1])
    blank_genome <- readr::read_csv("/home/timmonsc/testing/Thy-AdenoCA_pooled.csv")[,1:4]
    blank_genome$mb_bin <- rep(1:nrow(blank_genome))
    txdb <-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    all.genes <- base::suppressMessages(GenomicFeatures::genes(txdb))
    chrom_lengths <- stats::setNames(object = all.genes@seqinfo@seqlengths[1:24], all.genes@seqinfo@seqnames[1:24])
    chrom_lengths_mb_bins <- floor(chrom_lengths/1e6)+1

    tads <- read_csv("/home/timmonsc/testing/TADMap_geneset_hs.csv") %>%
    dplyr::select(-genelist) %>%
    tidyr::separate(tad, into=c('tad', 'start_chrom', 'start', 'end'), sep="\\|") %>%
    tidyr::separate(start_chrom, into=c('x', 'start_chrom'), sep="r") %>%
    dplyr::select(-x) %>%
    dplyr::mutate(start_chrom = dplyr::case_when(start_chrom == "X" ~ 23,
                                               start_chrom == "Y" ~ 24,
                                               TRUE ~ as.numeric(start_chrom)),
                tad = as.numeric(tad),
                start = as.numeric(start),
                end = as.numeric(end)) %>%
    dplyr::mutate(mb_bin_start = dplyr::case_when(start_chrom == 1 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6),
                                                start_chrom == 2 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                                start_chrom == 3 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                                start_chrom == 4 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                                start_chrom == 5 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                                start_chrom == 6 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                                start_chrom == 7 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                                start_chrom == 8 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                                start_chrom == 9 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                                start_chrom == 10 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                                start_chrom == 11 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                                start_chrom == 12 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                                start_chrom == 13 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                                start_chrom == 14 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                                start_chrom == 15 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                                start_chrom == 16 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                                start_chrom == 17 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                                start_chrom == 18 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                                start_chrom == 19 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                                start_chrom == 20 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                                start_chrom == 21 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                                start_chrom == 22 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                                start_chrom == 23 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                                start_chrom == 24 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23])),
                mb_bin_end = dplyr::case_when(start_chrom == 1 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6),
                                              start_chrom == 2 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                              start_chrom == 3 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                              start_chrom == 4 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                              start_chrom == 5 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                              start_chrom == 6 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                              start_chrom == 7 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                              start_chrom == 8 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                              start_chrom == 9 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                              start_chrom == 10 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                              start_chrom == 11 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                              start_chrom == 12 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                              start_chrom == 13 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                              start_chrom == 14 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                              start_chrom == 15 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                              start_chrom == 16 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                              start_chrom == 17 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                              start_chrom == 18 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                              start_chrom == 19 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                              start_chrom == 20 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                              start_chrom == 21 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                              start_chrom == 22 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                              start_chrom == 23 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                              start_chrom == 24 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23]))) %>%
    dplyr::arrange(start_chrom, start)
    tads <- tads %>%
    dplyr::mutate(mb_bin_after = dplyr::lead(tads$mb_bin_start)) %>%
    tidyr::replace_na(list(mb_bin_after = 3054)) %>%
    dplyr::mutate(width = mb_bin_after - mb_bin_end)
    pvals <- testAllSamples(type, blank_genome, tads)
    pvals <- pvals %>%
    dplyr::mutate("p_adj_bf" = p.adjust(pvals$pval, method="bonferroni"))
    write.csv(pvals, file=paste0("/home/timmonsc/testing/",type,"_tadpvals.txt"), row.names=FALSE)
}

main()
