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

# Automates testing
# Input folder containing signature profiles for a cancer type
# Performs KS test on changepoints in original sample
# Generates null distribution of 10,000 random changepoints and their KS p-values
# Compares sample p-values to null p-values
# Reports the proportion of null p-values that are lower than sample p-values
katDist <- function(type) {
  file_ids <- read.delim(paste0("/home/timmonsc/testing/bsub_logs/",type,"_ids.txt"), header=FALSE)
  colnames(file_ids) <- c("id")
  pvals <- data.frame()

  for (i in c(1:330)) {
    if (file.exists(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))) {
      if (file.exists(paste0("/home/timmonsc/testing/bsub_logs/",type,"_ids.txt"))) {
        print(i)

        load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))
        temp <- load(paste0("/home/timmonsc/testing/",type,"/",type,i,".Rdata"))

        if (!is.null(get(temp)[[2]])||!is.null(get(temp)[[6]])||!is.null(get(temp)[[10]])||!is.null(get(temp)[[14]])||!is.null(get(temp)[[18]])) {
          fileID <- trimws(as.character(file_ids$id[i]), which='both')
          if (file.exists(paste0("/home/timmonsc/testing/pcawg-kataegis/",fileID,"_Kataegis.tsv.gz"))) {
            # test for overlap

            changepoints_all <- summarizeChangepoints(get(temp), i) %>%
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

            kat <- makeKatDataset(fileID)
            kat_bins <- sort(unique(unlist(mapply(function(x,y){rep(floor(x):floor(y))}, kat$mb_bin_start, kat$mb_bin_end))))
            kat_presence <- blank_genome %>%
              dplyr::mutate(kat = dplyr::if_else(mb_bin %in% kat_bins, 1, 0))

            null_dist <- getRandomChangepoints(get(temp), kat_presence)

            all <- testSampleChangepoints(changepoints_all, kat_presence) %>%
              dplyr::mutate(
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
  }

  }
  return (pvals)
}

# Generate a null distribution of 10,000 randomly-placed changepoints
# Find the kataegis distributions within random changepoint-containing segment
# and regions surrounding changepoint
# Perform KS test between changepoint region and surrounding regions for
# each random changepoint
# Return p-values
getRandomChangepoints <- function(traj, kat_presence) {

  n = 10000

  if (nrow(kat_presence) > 0) {
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
      sequence <- (x-floor((y/2))):(x+floor((y/2)))
      sequence <- sequence[which(sequence >= 1)]
      sequence <- sequence[which(sequence <= 3113)]
      return (sequence)
    }

    before_segment <- mapply(seg_function, midpoints_before, seg_before_widths)
    cp_segment_before <- mapply(seg_function, middles, seg_before_widths)
    cp_segment_after <- mapply(seg_function, middles, seg_after_widths)
    after_segment <- mapply(seg_function, midpoints_after, seg_after_widths)


    before_segment_vals <- sapply(before_segment, function(x){kat_presence$kat[x]})
    cp_segment_before_vals <- sapply(cp_segment_before, function(x){kat_presence$kat[x]})
    cp_segment_after_vals <- sapply(cp_segment_after, function(x){kat_presence$kat[x]})
    after_segment_vals <- sapply(after_segment, function(x){kat_presence$kat[x]})

    before_segment_vals <- sapply(before_segment_vals, function(x){c(length(x[x==1]), length(x[x==0]))})
    cp_segment_before_vals <- sapply(cp_segment_before_vals, function(x){c(length(x[x==1]), length(x[x==0]))})
    cp_segment_after_vals <- sapply(cp_segment_after_vals, function(x){c(length(x[x==1]), length(x[x==0]))})
    after_segment_vals <- sapply(after_segment_vals, function(x){c(length(x[x==1]), length(x[x==0]))})

    before_within <- apply(before_segment_vals, 2, function(x,y){fisher.test(rbind(c(x[1],x[2]), c(y[1],y[2])))$p.value}, y = cp_segment_before_vals)
    after_within <- apply(after_segment_vals, 2, function(x,y){fisher.test(rbind(c(x[1],x[2]), c(y[1],y[2])))$p.value}, y = cp_segment_after_vals)

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
# Find the kataegis distributions within each changepoint-containing segment
# and regions surrounding changepoints
# Perform KS test between changepoint region and surrounding regions
# Return p values
testSampleChangepoints <- function(changepoints, kat_presence) {
  if (nrow(changepoints) > 0) {
    before_within <- rep(NA, nrow(changepoints))
    after_within <- rep(NA, nrow(changepoints))
    midpoints_before <- changepoints$mb_bin_start - ((changepoints$mb_bin_start-changepoints$mb_bin_before)/2)
    midpoints_within <- changepoints$mb_bin_end - ((changepoints$mb_bin_end-changepoints$mb_bin_start)/2)
    midpoints_after <- changepoints$mb_bin_after - ((changepoints$mb_bin_after-changepoints$mb_bin_end)/2)
    starts <- changepoints$mb_bin_start
    ends <- changepoints$mb_bin_end
    lengths_before <- (starts - midpoints_before) / 2
    lengths_after <- (midpoints_after - ends) / 2
    for (n in 1:nrow(changepoints)) {
      # find dist of atac scores before cp vs within cp
      before_segment1 <- which.min(abs((midpoints_before[n]-(lengths_before[n]/2))-kat_presence$mb_bin))
      before_segment2 <- which.min(abs((midpoints_before[n]+(lengths_before[n]/2))-kat_presence$mb_bin))
      before_segment <- kat_presence$kat[before_segment1:before_segment2]
      before_segment <- c(length(before_segment[before_segment==1]), length(before_segment[before_segment==0]))
      before_cp1 <- which.min(abs((midpoints_within[n]-(lengths_before[n]/2))-kat_presence$mb_bin))
      before_cp2 <- which.min(abs((midpoints_within[n]+(lengths_before[n]/2))-kat_presence$mb_bin))
      before_cp <- kat_presence$kat[before_cp1:before_cp2]
      before_cp <- c(length(before_cp[before_cp==1]), length(before_cp[before_cp==0]))
      # find dist of atac scores after cp vs within cp
      after_segment1 <- which.min(abs((midpoints_after[n]-(lengths_after[n]/2))-kat_presence$mb_bin))
      after_segment2 <- which.min(abs((midpoints_after[n]+(lengths_after[n]/2))-kat_presence$mb_bin))
      after_segment <- kat_presence$kat[after_segment1:after_segment2]
      after_segment <- c(length(after_segment[after_segment==1]), length(after_segment[after_segment==0]))
      after_cp1 <- which.min(abs((midpoints_within[n]-(lengths_after[n]/2))-kat_presence$mb_bin))
      after_cp2 <- which.min(abs((midpoints_within[n]+(lengths_after[n]/2))-kat_presence$mb_bin))
      after_cp <- kat_presence$kat[after_cp1:after_cp2]
      after_cp <- c(length(after_cp[after_cp==1]), length(after_cp[after_cp==0]))

      # fisher test
      before_within[n] <- fisher.test(rbind(before_cp, before_segment))$p.value
      after_within[n] <- fisher.test(rbind(after_cp, after_segment))$p.value
    }

    return (data.frame("before_within" = before_within,
                       "after_within" = after_within))
  }
  else {
    return (data.frame())
  }
}

# make clean dataset of kataegis events in a sample
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
    kat_test <- katDist(type)
    print(head(kat_test))
    write.csv(kat_test, file=paste0("/home/timmonsc/testing/",type,"_kataegis_windows_pvals.txt"), row.names=FALSE)
}

main()
