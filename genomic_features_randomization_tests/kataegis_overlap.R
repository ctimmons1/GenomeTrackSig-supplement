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

# Find the proportion of changepoints in a given sample that
# overlap with a kataegis event
generateSampleOverlaps <- function(fileID, traj, blank_genome) {
  # get changepoints
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)
  # get kat
  kat <- makeKatDataset(fileID)

  if (nrow(kat) > 0) {

    starts <- changepoints$mb_bin_start
    ends <- changepoints$mb_bin_end
    lower <- sapply(starts, function(x){which.min(abs(kat$mb_bin_end - x))})
    upper <- sapply(ends, function(x){which.min(abs(kat$mb_bin_start - x))})
    tf <- abs(starts - kat$mb_bin_end[lower]) > abs(ends - kat$mb_bin_start[upper])
    closest <- replace(lower, which(tf==TRUE), upper[which(tf==TRUE)])

    cp_range <- mapply(function(x, y){IRanges::IRanges(x:y)}, starts, ends)
    kat_range <- sapply(closest, function(x){kat$mb_bin_start[x]:kat$mb_bin_end[x]})
    kat_range <- sapply(kat_range, function(x){IRanges::IRanges(x)})
    overlaps <- mapply(function(x, y){IRanges::overlapsAny(x, y)}, cp_range, kat_range)
    overlaps <- sapply(overlaps, function(x){if (TRUE %in% x) {x = 1} else {x = 0}})
    overlap_prop <- length(overlaps[overlaps>0])/length(overlaps)

    return (overlap_prop)
  }
  else {
    return (rep(0, each=nrow(changepoints)))
  }
}

# Generate a null distribution of 10,000 samples with randomly-placed changepoints
# Find proportion of random changepoints that overlap with kataegis event
# Report proportion
generateRandomOverlaps <- function(fileID, traj, blank_genome) {
  # get changepoints -- to know number of cps in the sample
  changepoints <- summarizeChangepoints(traj, 1) %>%
    dplyr::mutate(width = mb_bin_end - mb_bin_start)
  # get kat
  kat <- makeKatDataset(fileID)

  n = 10000

  if (nrow(kat) > 0) {
    binData <- traj[['binData']]
    if (!is.null(binData$actual_bin)) {
      binData <- binData[order(nrow(binData):1),]
    }
    binData$mb_bin_start <- 1
    binData$mb_bin_end <- 1
    binData$mb_bin_before <- 1
    binData$mb_bin_after <- 1
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
      if (is.null(binData$actual_bin)) {
        if (i == nrow(binData)) {
          if (binData$end_chrom[i] == 24) {
            binData$end[i] <- 59373566
          }
        }
      }
      else {
        if (i == 1) {
          if (binData$end_chrom[i] == 24) {
            binData$end[i] <- 59373566
          }
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
    lower <- sapply(starts, function(x){which.min(abs(kat$mb_bin_end - x))})
    upper <- sapply(ends, function(x){which.min(abs(kat$mb_bin_start - x))})
    tf <- abs(starts - kat$mb_bin_end[lower]) > abs(ends - kat$mb_bin_start[upper])
    closest <- replace(lower, which(tf==TRUE), upper[which(tf==TRUE)])

    cp_range <- mapply(function(x, y){IRanges::IRanges(x:y)}, starts, ends)
    kat_range <- sapply(closest, function(x){kat$mb_bin_start[x]:kat$mb_bin_end[x]})
    kat_range <- sapply(kat_range, function(x){IRanges::IRanges(x)})
    overlaps <- mapply(function(x, y){IRanges::overlapsAny(x, y)}, cp_range, kat_range)
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
# the proportion of their changepoints that overlap with a kataegis event
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
          fileID <- trimws(as.character(file_ids$id[i]), which='both')
          if (file.exists(paste0("/home/timmonsc/testing/pcawg-kataegis/",fileID,"_Kataegis.tsv.gz"))) {
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

# make dataset of kataegis events in a sample
makeKatDataset <- function(fileID) {
  sample <- read.delim(paste0("/home/timmonsc/testing/pcawg-kataegis/",fileID,"_Kataegis.tsv.gz"))
  sample <- sample[,1:7] %>%
    dplyr::select(Chromosome, Start_Position, End_Position) %>%
    dplyr::mutate(Chromosome = dplyr::case_when(Chromosome == "X" ~ 23,
                                                Chromosome == "Y" ~ 24,
                                                TRUE ~ as.numeric(Chromosome)),
                  mb_bin_start = dplyr::case_when(Chromosome == 1 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6),
                                                  Chromosome == 2 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                                  Chromosome == 3 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                                  Chromosome == 4 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                                  Chromosome == 5 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                                  Chromosome == 6 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                                  Chromosome == 7 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                                  Chromosome == 8 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                                  Chromosome == 9 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                                  Chromosome == 10 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                                  Chromosome == 11 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                                  Chromosome == 12 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                                  Chromosome == 13 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                                  Chromosome == 14 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                                  Chromosome == 15 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                                  Chromosome == 16 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                                  Chromosome == 17 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                                  Chromosome == 18 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                                  Chromosome == 19 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                                  Chromosome == 20 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                                  Chromosome == 21 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                                  Chromosome == 22 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                                  Chromosome == 23 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                                  Chromosome == 24 ~ floor(Start_Position/1e6)+1 + ((Start_Position - (trunc((Start_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23])),
                  mb_bin_end = dplyr::case_when(Chromosome == 1 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6),
                                                Chromosome == 2 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                                Chromosome == 3 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                                Chromosome == 4 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                                Chromosome == 5 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                                Chromosome == 6 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                                Chromosome == 7 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                                Chromosome == 8 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                                Chromosome == 9 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                                Chromosome == 10 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                                Chromosome == 11 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                                Chromosome == 12 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                                Chromosome == 13 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                                Chromosome == 14 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                                Chromosome == 15 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                                Chromosome == 16 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                                Chromosome == 17 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                                Chromosome == 18 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                                Chromosome == 19 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                                Chromosome == 20 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                                Chromosome == 21 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                                Chromosome == 22 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                                Chromosome == 23 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                                Chromosome == 24 ~ (ceiling(End_Position/1e6)) + ((End_Position - (trunc((End_Position/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23]))) %>%
    tidyr::drop_na()

  return (sample)
}

main <- function() {
    args <- commandArgs(trailingOnly=TRUE)
    type <- as.character(args[1])
    blank_genome <- readr::read_csv("/home/timmonsc/testing/Thy-AdenoCA_pooled.csv")[,1:4]
    blank_genome$mb_bin <- rep(1:nrow(blank_genome))
    pvals <- testAllSamples(type, blank_genome)
    pvals <- pvals %>%
    dplyr::mutate("p_adj_bf" = p.adjust(pvals$pval, method="bonferroni"))
    write.csv(pvals, file=paste0("/home/timmonsc/testing/",type,"_kataegis_overlap_pvals.txt"), row.names=FALSE)
}

main()
