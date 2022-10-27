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

# make dataset of copy number across the genome
makeCnaDataset <- function(fileID, blank_genome) {
  cna_dataset <- data.frame()

    file <- read.delim(paste0("/home/timmonsc/testing/consensus.20170119.somatic.cna.annotated/",fileID,".consensus.20170119.somatic.cna.annotated.txt")) %>%
      dplyr::select(chromosome, start, end, total_cn) %>%
      dplyr::mutate(chromosome = dplyr::case_when(chromosome == "X" ~ 23,
                                                  chromosome == "Y" ~ 24,
                                                  TRUE ~ as.numeric(chromosome)),
                    total_cn = dplyr::case_when(total_cn == NA ~ 2,
                                                 TRUE ~ as.numeric(total_cn)),
           mb_bin_start = dplyr::case_when(chromosome == 1 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6),
                                           chromosome == 2 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                           chromosome == 3 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                           chromosome == 4 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                           chromosome == 5 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                           chromosome == 6 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                           chromosome == 7 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                           chromosome == 8 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                           chromosome == 9 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                           chromosome == 10 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                           chromosome == 11 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                           chromosome == 12 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                           chromosome == 13 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                           chromosome == 14 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                           chromosome == 15 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                           chromosome == 16 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                           chromosome == 17 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                           chromosome == 18 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                           chromosome == 19 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                           chromosome == 20 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                           chromosome == 21 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                           chromosome == 22 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                           chromosome == 23 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                           chromosome == 24 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23])),
           mb_bin_end = dplyr::case_when(chromosome == 1 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6),
                                         chromosome == 2 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                         chromosome == 3 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                         chromosome == 4 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                         chromosome == 5 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                         chromosome == 6 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                         chromosome == 7 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                         chromosome == 8 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                         chromosome == 9 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                         chromosome == 10 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                         chromosome == 11 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                         chromosome == 12 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                         chromosome == 13 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                         chromosome == 14 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                         chromosome == 15 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                         chromosome == 16 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                         chromosome == 17 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                         chromosome == 18 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                         chromosome == 19 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                         chromosome == 20 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                         chromosome == 21 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                         chromosome == 22 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                         chromosome == 23 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                         chromosome == 24 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23]))) %>%
      tidyr::drop_na()
    return (file)

}

# Find the proportion of changepoints in a given sample that
# overlap with a CNA
generateSampleOverlaps <- function(fileID, traj, blank_genome) {
  # get changepoints
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)
  # get CNAs
  cnas <- makeCnaDataset(fileID, blank_genome) %>%
    dplyr::mutate_at(c("total_cn"), ~replace(., is.na(.), 2)) %>%
    dplyr::filter(total_cn != 2) %>%
    tidyr::drop_na()

  if (nrow(cnas) > 0) {
    # vector to store overlaps of all changepoints in the sample
    sample_overlaps <- c()
    # find the proportion of overlap between each changepoint and its closest CNA
    for (i in 1:nrow(changepoints)) {
      # find closest CNA on either side of changepoint
      sample_cna_lower <- which.min(abs(changepoints$mb_bin_start[i]-cnas$mb_bin_end))
      sample_cna_higher <- which.min(abs((changepoints$mb_bin_end[i])-cnas$mb_bin_start))
      # find which CNA on either side of the changepoint is closest to the changepoint
      if (abs((cnas$mb_bin_start[sample_cna_higher]-changepoints$mb_bin_end[i])) < abs((cnas$mb_bin_end[sample_cna_lower]-changepoints$mb_bin_start[i]))) {
        sample_cna_index <- sample_cna_higher
      } else {
        sample_cna_index <- sample_cna_lower
      }
      # range of coordinates containing the changepoint and the CNA
      sample_cna_range <- IRanges::IRanges(rep(cnas$mb_bin_start[sample_cna_index]:cnas$mb_bin_end[sample_cna_index]))
      sample_cp_range <- IRanges::IRanges(rep(changepoints$mb_bin_start[i]:changepoints$mb_bin_end[i]))
      cp_overlaps <- IRanges::overlapsAny(sample_cna_range, sample_cp_range)
      if (TRUE %in% cp_overlaps) {
        cp_overlaps <- 1
      }
      # proportions of overlaps for each changepoint in the sample
      sample_overlaps <- c(sample_overlaps, cp_overlaps)
    }
    # we have the overlaps for all the changepoints in the sample
    # find the proportion of changepoints in the sample that overlap with a CNA
    sample_prop <- length(sample_overlaps[sample_overlaps>0])/length(sample_overlaps)
    return (sample_prop)
  }
  else {
    return (rep(0, each=nrow(changepoints)))
  }

}

# Generate a null distribution of 10,000 samples with randomly-placed changepoints
# Find proportion of random changepoints that overlap with CNA
# Report proportion
generateRandomOverlaps <- function(fileID, traj, blank_genome) {
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)

  # get cnas
  cnas <- makeCnaDataset(fileID, blank_genome) %>%
    dplyr::mutate_at(c("total_cn"), ~replace(., is.na(.), 2)) %>%
    dplyr::filter(total_cn != 2) %>%
    tidyr::drop_na()

  n = 10000

  if (nrow(cnas) > 0) {
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
    lower <- sapply(starts, function(x){which.min(abs(cnas$mb_bin_end - x))})
    upper <- sapply(ends, function(x){which.min(abs(cnas$mb_bin_start - x))})
    tf <- abs(starts - cnas$mb_bin_end[lower]) > abs(ends - cnas$mb_bin_start[upper])
    closest <- replace(lower, which(tf==TRUE), upper[which(tf==TRUE)])

    cp_range <- mapply(function(x, y){IRanges::IRanges(x:y)}, starts, ends)
    cna_range <- sapply(closest, function(x){cnas$mb_bin_start[x]:cnas$mb_bin_end[x]})
    cna_range <- sapply(cna_range, function(x){IRanges::IRanges(x)})
    overlaps <- mapply(function(x, y){IRanges::overlapsAny(x, y)}, cp_range, cna_range)
    overlaps <- sapply(overlaps, function(x){if (TRUE %in% x) {x = 1} else {x = 0}})
    overlaps <- matrix(unlist(overlaps), ncol = nrow(changepoints), nrow =n)
    overlap_props <- apply(overlaps, 1, function(x){length(x[x>0])/length(x)})

    return (overlap_props)
  }
  else {
    return (rep(0, each=n))
  }
}

# calculate the fraction of null proportions that are greater than sample proportion
testOverlapSignificance <- function(random_overlaps, sample_overlap) {
  p <- length(random_overlaps[random_overlaps>=sample_overlap])/length(random_overlaps)
  return (p)
}

# perform overlap test on a single sample
overlapTest <- function(fileID, traj, blank_genome) {
  random_overlaps <- generateRandomOverlaps(fileID, traj, blank_genome)
  sample_overlaps <- generateSampleOverlaps(fileID, traj, blank_genome)
  pval <- testOverlapSignificance(random_overlaps, sample_overlaps)
  return (pval)
}

# Automates testing
# Input folder containing signature profiles for a cancer type and template genome binned by Mb
# Performs overlap test on changepoints in original sample
# Generates null distribution of 10,000 random samples and
# the proportion of their changepoints that overlap with a CNA
# Compares sample proportion to null proportion
# Reports the fraction of null proportions that are lower than sample proportion
testAllSamples <- function(type, blank_genome) {
  pvals_all <- c()
  ids <- c()


  for (i in c(1:330)) {
    if (file.exists(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))) {
      if (file.exists(paste0("/home/timmonsc/testing/bsub_logs/",type,"_ids.txt"))) {
        print(i)
        file_ids <- read.delim(paste0("/home/timmonsc/testing/bsub_logs/",type,"_ids.txt"), header=FALSE)
        colnames(file_ids) <- c("id")

        load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))
        temp <- load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))

        # check for changepoints
        if (!is.null(get(temp)[[2]])||!is.null(get(temp)[[6]])||!is.null(get(temp)[[10]])||!is.null(get(temp)[[14]])||!is.null(get(temp)[[18]])) {
          # ||!is.null(get(temp)[[22]])||!is.null(get(temp)[[26]])||!is.null(get(temp)[[30]])||!is.null(get(temp)[[34]])||!is.null(get(temp)[[38]])) {
          fileID <- trimws(as.character(file_ids$id[i]), which='both')
          if (file.exists(paste0("/home/timmonsc/testing/consensus.20170119.somatic.cna.annotated/",fileID,".consensus.20170119.somatic.cna.annotated.txt"))) {
            # test for overlap
            ids <- c(ids, i)
            pval <- overlapTest(fileID, get(temp), blank_genome)
            pvals_all <- c(pvals_all, pval)
          }
        }
      }
    }
  }

  results <- data.frame('id' = ids, 'pval' = pvals_all)
  return (results)
}

main <- function() {
    args <- commandArgs(trailingOnly=TRUE)
    type <- as.character(args[1])
    blank_genome <- readr::read_csv("/home/timmonsc/testing/Thy-AdenoCA_pooled.csv")[,1:4]
    blank_genome$mb_bin <- rep(1:nrow(blank_genome))
    pvals <- testAllSamples(type, blank_genome)
    pvals <- pvals %>%
    dplyr::mutate("p_adj_bf" = p.adjust(pvals$pval, method="bonferroni"))
    write.csv(pvals, file=paste0("/home/timmonsc/testing/",type,"_pvals.txt"), row.names=FALSE)
}

main()
