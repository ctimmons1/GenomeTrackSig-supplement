########
# AUTHOR: Caitlin Timmons (caitlintimmons811@gmail.com)
# Contains helpful datasets and functions to analyze signature profiles and changepoints
########


# 1 Mb binned genomic coordinates
blank_genome <- readr::read_csv("~/Desktop/CBSP2021/Thy-AdenoCA_pooled.csv")[,1:4]
blank_genome$mb_bin <- rep(1:nrow(blank_genome))
change_bins <- c(1)
for (i in 2:nrow(blank_genome)-1){
  if (blank_genome$start_chrom[i] < blank_genome$start_chrom[i+1]) {
    change_bins <- c(change_bins, blank_genome$mb_bin[i+1])
  }
}
chr_labels <- as.character(c(1:22, "X", "Y"))

# statistical mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# return the cosine distances for each changepoint in a profile
cosineDist <- function(traj) {
  bootstraps <- length(traj) / 4
  distances <- data.frame(cpPos1 = c(), dist = c())
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        cosine_dist <- 1 - (as.numeric(lsa::cosine(traj[[i]][,j], traj[[i]][,j+1])))
        distances <- rbind(distances, data.frame(cpPos1=j, dist = cosine_dist))
      }
    }
    i <- i + 4
  }
  return (distances)
}

# find the amount by which each signature changes at all the
# changepoints in a profile
findSigChanges <- function(traj) {
  bootstraps <- length(traj) / 4
  sig_changes <- data.frame()
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        changes <- as.data.frame(traj[[i]][,j+1]-traj[[i]][,j])
        changes_flip <- data.table::transpose(changes)
        colnames(changes_flip) <- rownames(changes)
        changes_flip$cpPos1 <- j
        sig_changes <- rbind(sig_changes, changes_flip)
      }
    }
    i <- i + 4
  }
  return (sig_changes)
}

# find the signatures that have the highest activity on either side of
# a changepoint in a profile
findDominantSigs <- function(traj) {
  bootstraps <- length(traj) / 4
  dominant_sigs <- data.frame()
  i <- 1
  for (k in 1:bootstraps) {
    if (!is.null(traj[[i+1]])) {
      for (j in traj[[i+1]]) {
        dominant_pre <- rownames(traj[[i]])[as.numeric(which(traj[[i]][,j] == max(traj[[i]][,j])))]
        dominant_post <- rownames(traj[[i]])[as.numeric(which(traj[[i]][,j+1] == max(traj[[i]][,j+1])))]
        sample <- data.frame(dominant_pre = dominant_pre, dominant_post = dominant_post, cpPos1 = j)
        dominant_sigs <- rbind(dominant_sigs, sample)
      }
    }
    i <- i + 4
  }
  return (dominant_sigs)
}

# make dataframe of information about all changepoints in a trajectory
summarizeChangepoints <- function(trajectory, sampleID) {
  cpPos <- assignChangepoints(trajectory,0)
  if (nrow(cpPos) > 0) {
    cpPos$start_chrom <- 1
    cpPos$end_chrom <- 1
    cpPos$start <- 1
    cpPos$end <- 1
    cpPos$mb_bin_start <- 1
    cpPos$mb_bin_end <- 1
    cpPos$sampleID <- sampleID
    for (i in 1:nrow(cpPos)) {
      if (!is.null(trajectory[['binData']]$actual_bin)) {
        cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
        cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$actual_bin==cpPos$cpPos1[i]]
      }

      else {
        cpPos$start[i] <- trajectory[['binData']]$start[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$end[i] <- trajectory[['binData']]$end[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$start_chrom[i] <- trajectory[['binData']]$start_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
        cpPos$end_chrom[i] <- trajectory[['binData']]$end_chrom[trajectory[['binData']]$bin==cpPos$cpPos1[i]]
      }
      mb_bin_start <- blank_genome %>%
        dplyr::filter(start_chrom == cpPos$start_chrom[i] & start == cpPos$start[i])
      cpPos$mb_bin_start[i] <- mb_bin_start$mb_bin
      mb_bin_end <- blank_genome %>%
        dplyr::filter(end_chrom == cpPos$end_chrom[i] & end == cpPos$end[i])
      cpPos$mb_bin_end[i] <- mb_bin_end$mb_bin
    }
    return (cpPos)
  }
  else {
    return (NULL)
  }
}

# build dataset of information for all changepoints for a cancer type
buildChangepointData <- function(type, blank_genome) {
  # initialize dataframes to combine
  cpPos <- data.frame()
  cosine_distances <- data.frame(cpPos1 = c(), dist = c(), sampleID = c())
  sig_changes_full <- data.frame()
  dominant_sigs_full <- data.frame()
  mutation_counts <- c()
  j <- 1
  cp_id <- 1
  for (i in c(1:330)) {
    cps_sample <- data.frame()
    full_traj <- list()
    for (l in c('a', 'b', 'c', 'd', 'e')) {
      if (file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,l,".Rdata"))) {
        load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,l,".Rdata"))
        temp <- load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,l,".Rdata"))
        full_traj <- c(full_traj, get(temp))
      }
    }
    if (!file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,l,".Rdata"))) {
      if (file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))) {
        load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
        temp <- load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
        full_traj <- get(temp)
      }
    }
    cps_sample <- summarizeChangepoints(full_traj, i)
    if (!is.null(cps_sample)) {
      dominant_sigs <- findDominantSigs(full_traj)
      if (nrow(dominant_sigs)>0) {
        dominant_sigs$sampleID <- i
        dominant_sigs_full <- rbind(dominant_sigs_full, dominant_sigs)
      }
      sig_changes <- findSigChanges(full_traj)
      if (nrow(sig_changes)>0) {
        sig_changes$sampleID <- i
        sig_changes_full <- rbind(sig_changes_full, sig_changes)
      }
      dist <- cosineDist(full_traj)
      if (nrow(dist)>0) {
        dist$sampleID <- i
        cosine_distances <- rbind(cosine_distances, dist)
      }
    }
    cpPos <- rbind(cpPos, cps_sample)
    j <- j + 1
  }
  if (nrow(cpPos)>0) {
    # find average cosine distances at each cp across bootstraps
    avg_distances <- cosine_distances %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarize(avg_dist = mean(dist))
    # find average signature exposure changes at each cp across bootstraps
    avg_sig_changes <- sig_changes_full %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarise_at(dplyr::vars(colnames(sig_changes_full)[1]:colnames(sig_changes_full)[ncol(sig_changes_full)-2]), mean) %>%
      dplyr::select(-c(sampleID,cpPos1))
    # find average of dominant signatures at each cp across bootstraps
    avg_dominant_sigs <- dominant_sigs_full %>%
      dplyr::group_by(sampleID, cpPos1) %>%
      dplyr::summarize_at(dplyr::vars(dominant_pre:dominant_post), Mode)

    changepoint_summaries <- cpPos %>%
      dplyr::rename(start_bin = cpPos1,
                    end_bin = cpPos2,
                    bootstrap_support = prob) %>%
      dplyr::mutate("avg_cosine_dist" = avg_distances$avg_dist,
                    "dominant_pre" = avg_dominant_sigs$dominant_pre,
                    "dominant_post" = avg_dominant_sigs$dominant_post)
    changepoint_summaries <- cbind(changepoint_summaries, avg_sig_changes[,2:length(colnames(avg_sig_changes))])
    changepoint_summaries$cp_id <- rep(1:nrow(changepoint_summaries))
    return (changepoint_summaries)
  }
  else {
    return (NULL)
  }
}

# identify possible range of locations for each changepoint
# based on kernel density function with user-specified bandwidth
# bandwidth of 4 was most appropriate for all PCAWG cancer types
buildOverlapsData <- function(cpPos, bandwidth) {
  changepoints <- c()
  overlaps_full <- data.frame('sampleID' = c(),'sd_min' = c(), 'sd_max' = c())
  for (i in 1:nrow(cpPos)) {
    cp_range <- rep(cpPos$mb_bin_start[i]:cpPos$mb_bin_end[i])
    if (length(cp_range)>1) {
      kde <- ks::kde(cp_range, bandwidth)
      overlaps <- data.frame(sampleID = cpPos$sampleID[i],
                             sd_min = kde$eval.points[which.min(abs(kde$estimate-max(kde$estimate)))] - sd(kde$eval.points),
                             sd_max =  kde$eval.points[which.max(abs(kde$estimate+max(kde$estimate)))] + +sd(kde$eval.points))
      overlaps_full <- rbind(overlaps_full, overlaps)
    }
    else {
      cp_range <- c(cpPos$mb_bin_start[i], (cpPos$mb_bin_start[i]+.5), (cpPos$mb_bin_start[i]+.99))
      kde <- ks::kde(cp_range, bandwidth)
      overlaps <- data.frame(sampleID = cpPos$sampleID[i],
                             sd_min = kde$eval.points[which.min(abs(kde$estimate-max(kde$estimate)))] - sd(kde$eval.points),
                             sd_max =  kde$eval.points[which.max(abs(kde$estimate+max(kde$estimate)))] + +sd(kde$eval.points))
      overlaps_full <- rbind(overlaps_full, overlaps)
    }

  }

  return (overlaps_full)
}


# Plotting code to visualize changepoint distributions across the genome
# for all samples from a particular tissue
plotChangepointSummary <- function(cpPos, overlaps_full=NULL, blank_genome, title, subtitle) {

  change_bins <- c(1)
  blank_genome <- readr::read_csv("~/Desktop/CBSP2021/Thy-AdenoCA_pooled.csv")[,1:4]
  blank_genome$mb_bin <- rep(1:nrow(blank_genome))
  for (i in 2:nrow(blank_genome)-1){
    if (blank_genome$start_chrom[i] < blank_genome$start_chrom[i+1]) {
      change_bins <- c(change_bins, blank_genome$mb_bin[i+1])
    }
  }
  chr_labels <- as.character(c(1:22, "X", "Y"))

  g <- ggplot2::ggplot(data = blank_genome, ggplot2::aes(x = mb_bin)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
    ggplot2::labs(x = "Chromosome", title = title,
                  subtitle = subtitle, y = "Sample") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())

  for (i in 1:length(change_bins)-1) {
    if (i %% 2 != 0) {
      g <- g + ggplot2::annotate("rect", xmin=change_bins[i], xmax = change_bins[i+1],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill='lightgrey')
    }
  }
  if (!is.null(overlaps_full)) {
    for (i in 1:nrow(overlaps_full)) {
      g <- g + ggplot2::annotate("rect", xmin = overlaps_full$sd_min[i], xmax = overlaps_full$sd_max[i],
                                 ymin=overlaps_full$tempID[i], ymax=overlaps_full$tempID[i]+1, fill="orange", alpha=0.3)
    }
  }
  for (i in 1:nrow(cpPos)) {
    g <- g + ggplot2::annotate("rect", xmax=cpPos$mb_bin_end[i], xmin=cpPos$mb_bin_start[i],
                               ymin=cpPos$tempID[i], ymax=cpPos$tempID[i]+1, alpha=cpPos$bootstrap_support[i], fill = "red")
  }
  g <- g + ggplot2::annotate("rect", xmax=3100, xmin=5,
                             ymin=max(cpPos$tempID), ymax = max(cpPos$tempID)+1, alpha=0, fill='white')
  g
}

# IDENTIFY CHANGEPOINT PEAKS (basis for recurrent changepoints)
calculatePeaks <- function(overlaps) {
  peaks <- data.frame("start" = c(), "stop" = c(), "n_cps" = c())
  for (k in c(1:3109)) {
    peak <- nrow(overlaps %>%
                   dplyr::filter((sd_min>=k & sd_min<=(k+5)) | (sd_max>=k & sd_max<=(k+5)) | (sd_min<=k & sd_max >= (k+5))))
    peaks <- rbind(peaks, data.frame("start" = k, "stop" = k+5, "n_cps" = peak))
  }
  peaks <- peaks %>%
    dplyr::filter(n_cps > 2)

  return (peaks)
}

# match each changepoint to its corresponding recurrent changepoint
# (if applicable)
wranglePeaks <- function(peaks, overlaps, changepoints, n, y) {
  peaks <- peaks %>%
    dplyr::mutate('end' = start+5,
                  'sample_count' = 0)
  overlaps$peak <- 0
  for (k in 1:nrow(peaks)) {
    peak <- overlaps %>%
      dplyr::filter((sd_min>=peaks$start[k] & sd_min<=peaks$end[k]) | (sd_max>=peaks$start[k] & sd_max<=peaks$end[k]) | (sd_min<=peaks$start[k] & sd_max >=peaks$end[k]))
    sample_count <- peak %>%
      dplyr::group_by(sampleID) %>%
      dplyr::summarize(n = dplyr::n_distinct())
    peaks$sample_count[k] <- nrow(sample_count)
    for (j in 1:nrow(overlaps)) {
      if (overlaps$sd_min[j] %in% peak$sd_min) {
        overlaps$peak[j] <- peaks$start[k]
      }
    }
  }
  sig_changes <- changepoints[,15:ncol(changepoints)]
  max_increase <- names(sig_changes)[apply(sig_changes, 1, which.max)]
  max_decrease <- names(sig_changes)[apply(sig_changes, 1, which.min)]
  sig_changes <- sig_changes %>%
    dplyr::mutate(max_increase = max_increase,
                  max_decrease = max_decrease)


  plot_table <- changepoints %>%
    cbind(overlaps[,c(2,3,4)])
    cbind(sig_changes[,c(ncol(sig_changes),ncol(sig_changes)-1)])
  return (plot_table)
}


# make dataframes with locations of repeated changepoint (found in >1)
# start locations for all applicable cancer types
melanoma_peaks <- data.frame('start' = c(47,191,481,565,703,898,908,1008,1050,1113,1195,1243,
                                         1389,1521,1531,1167,1760,1797,1881,1900,1928,2055,2141,2198,2294,
                                         2496,2560,2696,2826))
lung_scc_peaks <- data.frame('start' = c(47,539,671,689,913,1109,1561,2114,2149,2180,2188,2214,2241,2873,2885,2897))
bladder_peaks <- data.frame('start' = c(47,539,671,689,913,1109,1561,2149,2162,2180,2188,2214,2241,2672,2788,2873,2885,2897,2928))
eso_peaks <- data.frame('start' = c(42,256,674,689,707,1532,1562,1576,1620,2051,2084,2092,2138,2181,2189,2256,2290,2348,2352,2358,
                                    2488,2495,2563,2581,2662,2790,2824,2836,2876,2895,2922,2938,2946))
lymph_bnhl_peaks <- data.frame('start' = c(707,911,2181,2194,2201,2208,2248,2275,2290,2308,2840,2864,2906))
prostate_peaks <- data.frame('start' = c(1392,1703))
stomach_peaks <- data.frame('start' = c(47,2292,2907,2917))
lymph_cll_peaks <- data.frame('start' = c(426,457,485,524,580,601,612,624,656,668,675,832,846,861,867,875,908,1097,1367,
                                    1435,1514,1536,1557,1581,1607,1776,2182,2238,2264,2331,2494,2623,2657,2734,2835))
cns_peaks <- data.frame('start' = c(66,251,704,905,1170,1290,1351,1875,1926,2142,2185,2249,2269,2284,2335,
                                    2855,2913))
colorect_peaks <- data.frame('start'= c(42,58,104,191,388,478,567,596,604,702,754,826,897,908,
                                        1054,1110,1193,1214,1237,1249,1293,1391,1452,1532,1669,1725,
                                        1807,1899,1931,1967,2062,2142,2179,2198,2254,2294,2480,2615,
                                        2659,2696,2830,2967))
uterus_peaks <- data.frame('start' = c(47,191,481,705,891,1052,1110,1195,1243,1405,1536,1669,1736,1818,1895,2015,
                                       2062,2145,2188,2297,2419,2480,2695,2782,2824))
cervix_peaks <- data.frame('start' = c(2933))
breast_adenoca_peaks <- data.frame('start' = c(1205,2894))
panc_adenoca_peaks <- data.frame('start' = c(2022,2494,2536,2664))
bone_osteosarc_peaks <- data.frame('start' =c(1088,1102,2068))
