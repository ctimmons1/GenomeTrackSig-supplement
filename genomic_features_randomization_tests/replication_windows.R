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

# make replication timing dataset across the genome
makeRep <- function() {
  rep <- read.delim("/home/timmonsc/testing/iPSC_individual_level_data.txt") %>%
    dplyr::rename(chr = X.chr) %>%
    dplyr::rowwise(chr, start, end, ID)
  rep <- rep %>%
    dplyr::mutate(median = stats::median(dplyr::c_across(colnames(rep)[5]:colnames(rep)[304]), na.rm=TRUE)) %>%
    dplyr::select(chr, start, end, ID, median) %>%
    mutate(chr = as.numeric(chr),
           mb_bin_start = dplyr::case_when(chr == 1 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6),
                                           chr == 2 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                           chr == 3 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                           chr == 4 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                           chr == 5 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                           chr == 6 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                           chr == 7 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                           chr == 8 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                           chr == 9 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                           chr == 10 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                           chr == 11 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                           chr == 12 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                           chr == 13 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                           chr == 14 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                           chr == 15 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                           chr == 16 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                           chr == 17 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                           chr == 18 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                           chr == 19 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                           chr == 20 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                           chr == 21 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                           chr == 22 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                           chr == 23 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                           chr == 24 ~ floor(start/1e6)+1 + ((start - (trunc((start/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23])),
           mb_bin_end = dplyr::case_when(chr == 1 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6),
                                         chr == 2 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1]),
                                         chr == 3 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:2]),
                                         chr == 4 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:3]),
                                         chr == 5 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:4]),
                                         chr == 6 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:5]),
                                         chr == 7 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:6]),
                                         chr == 8 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:7]),
                                         chr == 9 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:8]),
                                         chr == 10 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:9]),
                                         chr == 11 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:10]),
                                         chr == 12 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:11]),
                                         chr == 13 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:12]),
                                         chr == 14 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:13]),
                                         chr == 15 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:14]),
                                         chr == 16 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:15]),
                                         chr == 17 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:16]),
                                         chr == 18 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:17]),
                                         chr == 19 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:18]),
                                         chr == 20 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:19]),
                                         chr == 21 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:20]),
                                         chr == 22 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:21]),
                                         chr == 23 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:22]),
                                         chr == 24 ~ (ceiling(end/1e6)) + ((end - (trunc((end/1e6),digits=0)*1e6))/1e6) + sum(chrom_lengths_mb_bins[1:23]))) %>%
    dplyr::arrange(chr, start)
  return (rep)
}

# Automates testing
# Input folder containing signature profiles for a cancer type and replication timing dataset
# Performs KS test on changepoints in original sample
# Generates null distribution of 10,000 random changepoints and their KS p-values
# Compares sample p-values to null p-values
# Reports the proportion of null p-values that are lower than sample p-values
repTest <- function(type, rep) {
  pvals <- data.frame()

  for (i in c(1:330)) {
    if (file.exists(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))) {
      print(i)

      load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))
      temp <- load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))

      # check for changepoints
      if (!is.null(get(temp)[[2]])||!is.null(get(temp)[[6]])||!is.null(get(temp)[[10]])||!is.null(get(temp)[[14]])||!is.null(get(temp)[[18]])) {
        changepoints_all <- summarizeChangepoints(get(temp), i) %>%
          dplyr::filter(start_chrom < 23 & end_chrom < 23) %>%
          dplyr::arrange(cpPos1)
        if (nrow(changepoints_all) > 0) {
              if (nrow(changepoints_all) == 1) {
                changepoints_all$mb_bin_before[1] <- min(blank_genome$mb_bin[blank_genome$start_chrom==changepoints_all$start_chrom[1]])
                changepoints_all$mb_bin_after[1] <- max(blank_genome$mb_bin[blank_genome$start_chrom==changepoints_all$end_chrom[1]])
              }
              else {
              for (n in 2:nrow(changepoints_all)) {
                if (changepoints_all$mb_bin_start[n] > changepoints_all$mb_bin_end[n-1]) {
                  changepoints_all$mb_bin_before[n] = changepoints_all$mb_bin_end[n-1]
                }
                else {
                  changepoints_all$mb_bin_before[n] = min(blank_genome$mb_bin[blank_genome$start_chrom==changepoints_all$start_chrom[n]])
                }
              }
              for (n in 1:((nrow(changepoints_all))-1)) {
                if (changepoints_all$mb_bin_end[n] < changepoints_all$mb_bin_start[n+1]) {
                  changepoints_all$mb_bin_after[n] = changepoints_all$mb_bin_start[n+1]
                }
                else {
                  changepoints_all$mb_bin_after[n] = max(blank_genome$mb_bin[blank_genome$start_chrom==changepoints_all$end_chrom[n]])
                }
              }
              changepoints_all$mb_bin_after[nrow(changepoints_all)] <- max(blank_genome$mb_bin[blank_genome$start_chrom==changepoints_all$end_chrom[nrow(changepoints_all)]])
              }

          changepoints_full <- changepoints_all %>%
            dplyr::filter(before_chrom == start_chrom & after_chrom == end_chrom)
          changepoints_truc <- changepoints_all %>%
            dplyr::filter(before_chrom != start_chrom | after_chrom != end_chrom) %>%
            dplyr::mutate(before_chrom = start_chrom,
                          after_chrom = end_chrom)
          if (nrow(changepoints_truc) > 0) {
            for (j in 1:nrow(changepoints_truc)) {
              changepoints_truc$start_before[j] <- min(blank_genome$start[blank_genome$start_chrom==changepoints_truc$before_chrom[j]])
              changepoints_truc$end_after[j] <- max(blank_genome$end[blank_genome$end_chrom==changepoints_truc$after_chrom[j]])
              mb_bin_before = blank_genome %>%
                dplyr::filter(start_chrom == changepoints_truc$before_chrom[j]) %>%
                dplyr::filter(start == changepoints_truc$start_before[j])
              changepoints_truc$mb_bin_before[j] = mb_bin_before$mb_bin
              mb_bin_after = blank_genome %>%
                dplyr::filter(end_chrom == changepoints_truc$after_chrom[j]) %>%
                dplyr::filter(end == changepoints_truc$end_after[j])
              changepoints_truc$mb_bin_after[j] = mb_bin_after$mb_bin
            }
          }

          changepoints_all <- rbind(changepoints_full, changepoints_truc) %>%
            dplyr::arrange(cpPos1)

          null_dist <- getRandomChangepoints(get(temp), rep)

          all <- testSampleChangepoints(changepoints_all, rep) %>%
            dplyr::mutate("start_bin" = changepoints_all$cpPos1,
                          "sampleID" = i,
                          "null_before" = 1,
                          "null_after" = 1)

          for (j in 1:nrow(all)) {
            all$null_before[j] = length(which(null_dist$before_within<=all$before_within[j]))/nrow(null_dist)
            all$null_after[j] = length(which(null_dist$after_within<=all$after_within[j]))/nrow(null_dist)
          }
          pvals <- rbind(pvals, all)
        }
      }
    }
  }
  return (pvals)
}

# Generate a null distribution of 10,000 randomly-placed changepoints
# Find the replication timing distributions within random changepoint-containing segment
# and regions surrounding changepoint
# Perform KS test between changepoint region and surrounding regions for
# each random changepoint
# Return p-values
getRandomChangepoints <- function(traj, rep) {

  n = 10000

  if (nrow(rep) > 0) {
    binData <- traj[['binData']] %>%
      dplyr::filter(start_chrom < 23)
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
    if (!is.null(binData$actual_bin)) {
      bin_boundaries = assignChromosomeBounds(traj, T)
    }
    else {
      bin_boundaries = assignChromosomeBounds(traj, F)
    }

    indices <- sample(2:(nrow(binData)-1), n, replace=TRUE)
    chroms <- sapply(indices, function(x){which.max(bin_boundaries[bin_boundaries<=x])})

    truncate_before <- function(indices, chroms) {
      min = bin_boundaries[chroms]
      if ((indices-2) > min) {
        return (sample(rep(min:(indices-2)), size=1))
      }
      else {
        return(min)
      }
    }

    truncate_after <- function(indices, chroms) {
      max = (bin_boundaries[chroms+1])-1
      max = rapply(list(max), function(x)ifelse(is.na(x),nrow(binData),x), how = "replace")
      if ((indices+2) < unlist(max)) {
        return (sample(rep((indices+2):unlist(max)), size=1))
      }
      else {
        return(max)
      }
    }

    indices_before <- unlist(mapply(truncate_before, indices, chroms))
    indices_after <- unlist(mapply(truncate_after, indices, chroms))
    ends <- binData$mb_bin_end[indices]
    starts <- binData$mb_bin_start[indices]
    middles <- starts + floor((ends-starts)/2)
    befores <- binData$mb_bin_start[indices_before]
    afters <- binData$mb_bin_end[indices_after]
    midpoints_before <- middles - floor((middles-befores)/2)
    midpoints_after <- afters - floor((afters-middles)/2)
    seg_before_widths <- floor((middles-befores)/2)
    seg_after_widths <- floor((afters-middles)/2)

    seg_function <- function(x, y) {
      sequence <- c((x-floor((y/2))),(x+floor((y/2))))
      for (i in 1:length(sequence)) {
        if (sequence[i] < 1) {
          sequence[i] <- 1
        }
        if (sequence[i] > 3113) {
          sequence[i] <- 3113
        }
      }
      return (sequence)
    }

    before_segment <- t(mapply(seg_function, midpoints_before, seg_before_widths))
    cp_segment_before <- t(mapply(seg_function, middles, seg_before_widths))
    cp_segment_after <- t(mapply(seg_function, middles, seg_after_widths))
    after_segment <- t(mapply(seg_function, midpoints_after, seg_after_widths))

    before_segment_vals <- apply(before_segment, 1, function(x){rep$median[(which.min(abs(x[1]-rep$mb_bin_start))):((which.min(abs((x[2]+1)-rep$mb_bin_end)))-1)]})
    cp_segment_before_vals <- apply(cp_segment_before, 1, function(x){rep$median[(which.min(abs(x[1]-rep$mb_bin_start))):((which.min(abs((x[2]+1)-rep$mb_bin_end)))-1)]})
    cp_segment_after_vals <- apply(cp_segment_after, 1, function(x){rep$median[(which.min(abs(x[1]-rep$mb_bin_start))):((which.min(abs((x[2]+1)-rep$mb_bin_end)))-1)]})
    after_segment_vals <- apply(after_segment, 1, function(x){rep$median[(which.min(abs(x[1]-rep$mb_bin_start))):((which.min(abs((x[2]+1)-rep$mb_bin_end)))-1)]})

    before_within <- mapply(function(x,y){ks.test(x, y)$p.value}, before_segment_vals, cp_segment_before_vals)
    after_within <- mapply(function(x,y){ks.test(x, y)$p.value}, after_segment_vals, cp_segment_after_vals)

    results <- data.frame("before_within" = before_within,
                          "after_within" = after_within)

    return (results)
  }

  else {
    return (data.frame("before_within" = rep(1, each=n),
                       "after_within" = rep(1, each=n)))
  }
}

# Applied to a single sample with changepoints
# Find the replication timing distributions within each changepoint-containing segment
# and regions surrounding changepoints
# Perform KS test between changepoint region and surrounding regions
# Return p values
testSampleChangepoints <- function(changepoints, rep) {
  changepoints <- changepoints %>%
    filter(end_chrom < 23)
  before_within <- rep(NA, nrow(changepoints))
  after_within <- rep(NA, nrow(changepoints))
  midpoints_before <- changepoints$mb_bin_start - ((changepoints$mb_bin_start-changepoints$mb_bin_before)/2)
  midpoints_within <- changepoints$mb_bin_end - ((changepoints$mb_bin_end-changepoints$mb_bin_start)/2)
  midpoints_after <- changepoints$mb_bin_after - ((changepoints$mb_bin_after-changepoints$mb_bin_end)/2)
  starts <- changepoints$mb_bin_start
  ends <- changepoints$mb_bin_end
  midpoints_after <- changepoints$mb_bin_after - ((changepoints$mb_bin_after-changepoints$mb_bin_end)/2)
  lengths_before <- (starts - midpoints_before) / 2
  lengths_after <- (midpoints_after - ends) / 2
  for (n in 1:nrow(changepoints)) {
    before_segment1 <- which.min(abs((midpoints_before[n]-(lengths_before[n]/2))-rep$mb_bin_start))
    before_segment2 <- which.min(abs((midpoints_before[n]+(lengths_before[n]/2))-rep$mb_bin_start))
    before_segment <- rep$median[before_segment1:before_segment2]
    before_cp1 <- which.min(abs((midpoints_within[n]-(lengths_before[n]/2))-rep$mb_bin_start))
    before_cp2 <- which.min(abs((midpoints_within[n]+(lengths_before[n]/2))-rep$mb_bin_start))
    before_cp <- rep$median[before_cp1:before_cp2]
    # find dist of atac medians after cp vs within cp
    after_segment1 <- which.min(abs((midpoints_after[n]-(lengths_after[n]/2))-rep$mb_bin_start))
    after_segment2 <- which.min(abs((midpoints_after[n]+(lengths_after[n]/2))-rep$mb_bin_start))
    after_segment <- rep$median[after_segment1:after_segment2]
    after_cp1 <- which.min(abs((midpoints_within[n]-(lengths_after[n]/2))-rep$mb_bin_start))
    after_cp2 <- which.min(abs((midpoints_within[n]+(lengths_after[n]/2))-rep$mb_bin_start))
    after_cp <- rep$median[after_cp1:after_cp2]

    # ks test
    before_within[n] <- ks.test(before_segment, before_cp)$p.value
    after_within[n] <- ks.test(after_segment, after_cp)$p.value

  }

  return (data.frame("before_within" = before_within,
                     "after_within" = after_within))
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
  rep <- makeRep()
  rep_test <- repTest(type, rep)
  write.csv(rep_test, file=paste0("/home/timmonsc/testing/",type,"_replicationpvals.txt"), row.names=FALSE)
}

main()

