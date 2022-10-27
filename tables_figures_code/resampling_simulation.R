# load packages
library(tidyverse)
devtools::load_all()
# import changepoint analysis functions not defined in package
source("../helper_functions.R")

# original cll recurrent changepoints
cll_cps <- buildChangepointData('lymph-cll', blank_genome)
cll_cps <- cll_cps %>%
  dplyr::mutate("type" = 'lymph-cll',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
cll_cps <- cll_cps[,c(1:16, (ncol(cll_cps)-2), (ncol(cll_cps)-1),ncol(cll_cps))]
cll_all <- wranglePeaks(cll_peaks, buildOverlapsData(cll_cps, 4), cll_cps, 95) %>%
  filter(sampleID != 83)

cll_sample_counts = cll_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

cll_all_recurrents <- cll_all %>%
  filter(peak != 0) %>%
  left_join(cll_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped cll changepoints
cll_bootstrap_cps <- buildChangepointData('cll_20bootstrap', blank_genome)
cll_all_bootstrap <- wranglePeaks(cll_peaks, buildOverlapsData(cll_bootstrap_cps, 4),
                              cll_bootstrap_cps, 95) %>%
  mutate(type='lymph-cll')

# look at changepoint confidence levels
hist(cll_bootstrap_cps$bootstrap_support)

ggplot(data = cll_bootstrap_cps, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = cll_cps, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "N Changepoints across all samples", title="Lymph-CLL changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
cll_sample_counts_bootstrap = cll_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

cll_all_recurrents_bootstrap <- cll_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(cll_sample_counts_bootstrap, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# 34/35 peaks  were recovered at least once
# Mb 668 (3q) peak is missing from bootstrap samples, but was
# only hit once

# all recurrent changepoints with >= 7 samples hit were recovered with >= 7
# samples hit in the bootstrap analysis
cll_recurrent_matches = cll_sample_counts %>%
  left_join(cll_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
cll_to_join = cll_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = cll_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

cll_all_recurrents <- cll_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(cll_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# all TRUE values: all bootstrap samples recover at least the same recurrent
# changepoints as the original samples

# compare changepoints with high confidence
cll_bootstrap_highconf = cll_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

cll_orig_highconf = cll_all %>%
  filter(sampleID %in% cll_all_bootstrap$sampleID & bootstrap_support > 0.1) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(cll_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(cll_orig_highconf$match), "of", nrow(cll_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(cll_recurrent_matches$sample_count.y)), "of", nrow(cll_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(cll_all_recurrents$match), "of", nrow(cll_all_recurrents), "sample-recurrent changepoint matches"))

# original melanoma recurrent changepoints
melanoma_cps <- buildChangepointData('melanoma', blank_genome)
melanoma_cps <- melanoma_cps %>%
  dplyr::mutate("type" = 'melanoma',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
melanoma_cps <- melanoma_cps[,c(1:16, (ncol(melanoma_cps)-2), (ncol(melanoma_cps)-1),ncol(melanoma_cps))]
melanoma_all <- wranglePeaks(melanoma_peaks, buildOverlapsData(melanoma_cps, 4), melanoma_cps, 107) %>%
  filter(sampleID %in% melanoma_all_bootstrap$sampleID)

melanoma_sample_counts = melanoma_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

melanoma_all_recurrents <- melanoma_all %>%
  filter(peak != 0) %>%
  left_join(melanoma_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped cll changepoints
melanoma_bootstrap_cps <- buildChangepointData('melanoma_20bootstrap', blank_genome)
melanoma_all_bootstrap <- wranglePeaks(melanoma_peaks, buildOverlapsData(melanoma_bootstrap_cps, 4),
                              melanoma_bootstrap_cps, 107) %>%
  mutate(type='melanoma')

# look at changepoint confidence levels
hist(melanoma_bootstrap_cps$bootstrap_support)

ggplot(data = melanoma_bootstrap_cps, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = melanoma_cps, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "N Changepoints across all samples", title="Melanoma changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
melanoma_sample_counts_bootstrap = melanoma_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

melanoma_all_recurrents_bootstrap <- melanoma_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(melanoma_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# 23/23 peaks  were recovered at least once
# Mb 668 (3q) peak is missing from bootstrap samples, but was
# only hit once

# all recurrent changepoints with >= 7 samples hit were recovered with >= 7
# samples hit in the bootstrap analysis
melanoma_recurrent_matches = melanoma_sample_counts %>%
  left_join(melanoma_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
melanoma_to_join = melanoma_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = melanoma_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

melanoma_all_recurrents <- melanoma_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(melanoma_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

melanoma_all_recurrents_final <- melanoma_all_recurrents_final %>%
  filter(!match)

# compare changepoints with high confidence
melanoma_bootstrap_highconf = melanoma_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

melanoma_orig_highconf = melanoma_all %>%
  filter(sampleID %in% melanoma_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(melanoma_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(melanoma_orig_highconf$match, na.rm=TRUE), "of", nrow(melanoma_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(melanoma_recurrent_matches$sample_count.y)), "of", nrow(melanoma_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(melanoma_all_recurrents$match), "of", nrow(melanoma_all_recurrents), "sample-recurrent changepoint matches"))

# all TRUE values: all bootstrap samples recover at least the same recurrent
# changepoints as the original samples

# original melanoma recurrent changepoints
eso_cps <- buildChangepointData('eso-adenoca', blank_genome)
eso_cps <- eso_cps %>%
  dplyr::mutate("type" = 'eso-adenoca',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
eso_cps <- eso_cps[,c(1:16, (ncol(eso_cps)-2), (ncol(eso_cps)-1),ncol(eso_cps))]
eso_all <- wranglePeaks(eso_peaks, buildOverlapsData(eso_cps, 4), eso_cps, 98) %>%
  filter(sampleID %in% eso_bootstrap_cps$sampleID)

eso_sample_counts = eso_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

eso_all_recurrents <- eso_all %>%
  filter(peak != 0) %>%
  left_join(eso_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped cll changepoints
eso_bootstrap_cps <- buildChangepointData('eso_20bootstrap', blank_genome)
eso_all_bootstrap <- wranglePeaks(eso_peaks, buildOverlapsData(eso_bootstrap_cps, 4),
                              eso_bootstrap_cps, 98) %>%
  mutate(type='eso-adenoca')

# look at changepoint confidence levels
hist(eso_bootstrap_cps$bootstrap_support)

ggplot(data = eso_bootstrap_cps, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = eso_cps, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "N Changepoints across all samples", title="Eso-AdenoCA changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
eso_sample_counts_bootstrap = eso_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

eso_all_recurrents_bootstrap <- eso_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(eso_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# 23/23 peaks  were recovered at least once
# Mb 668 (3q) peak is missing from bootstrap samples, but was
# only hit once

# all recurrent changepoints with >= 7 samples hit were recovered with >= 7
# samples hit in the bootstrap analysis
eso_recurrent_matches = eso_sample_counts %>%
  left_join(eso_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
eso_to_join = eso_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = eso_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

eso_all_recurrents_final <- eso_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(eso_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# compare changepoints with high confidence
eso_bootstrap_highconf = eso_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

eso_orig_highconf = eso_all %>%
  filter(sampleID %in% eso_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(eso_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(eso_orig_highconf$match, na.rm=TRUE), "of", nrow(eso_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(eso_recurrent_matches$sample_count.y)), "of", nrow(eso_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(eso_all_recurrents_final$match), "of", nrow(eso_all_recurrents_final), "sample-recurrent changepoint matches"))

######### CNS_GBM #############


# original cns-gbm recurrent changepoints
cns_cps <- buildChangepointData('cns-gbm', blank_genome)
cns_cps <- cns_cps %>%
  dplyr::mutate("type" = 'cns-gbm',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
cns_cps <- cns_cps[,c(1:16, (ncol(cns_cps)-2), (ncol(cns_cps)-1),ncol(cns_cps))]
cns_all <- wranglePeaks(cns_peaks, buildOverlapsData(cns_cps, 4), cns_cps, 41) %>%
  filter(sampleID %in% cns_all_bootstrap$sampleID)

cns_sample_counts = cns_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

cns_all_recurrents <- cns_all %>%
  filter(peak != 0) %>%
  left_join(cns_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped cns-gbm changepoints
cns_bootstrap_cps <- buildChangepointData('cns_20bootstrap', blank_genome)
cns_all_bootstrap <- wranglePeaks(cns_peaks, buildOverlapsData(cns_bootstrap_cps, 4),
                                  cns_bootstrap_cps, 41) %>%
  mutate(type='cns-gbm')

# look at changepoint confidence levels
hist(cns_bootstrap_cps$bootstrap_support)

ggplot(data = cns_all_bootstrap, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = cns_all, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "Changepoint density across all samples", title="CNS-GBM changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
cns_sample_counts_bootstrap = cns_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

cns_all_recurrents_bootstrap <- cns_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(cns_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

cns_recurrent_matches = cns_sample_counts %>%
  left_join(cns_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
cns_to_join = cns_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = cns_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

cns_all_recurrents_final <- cns_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(cns_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# compare changepoints with high confidence
cns_bootstrap_highconf = cns_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

cns_orig_highconf = cns_all %>%
  filter(sampleID %in% cns_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(cns_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(cns_orig_highconf$match, na.rm=TRUE), "of", nrow(cns_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(cns_recurrent_matches$sample_count.y)), "of", nrow(cns_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(cns_all_recurrents_final$match), "of", nrow(cns_all_recurrents_final), "sample-recurrent changepoint matches"))

######### LUNG-SCC #############


# original lung-scc recurrent changepoints
scc_cps <- buildChangepointData('lung-scc', blank_genome)
scc_cps <- scc_cps %>%
  dplyr::mutate("type" = 'lung-scc',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
scc_cps <- scc_cps[,c(1:16, (ncol(scc_cps)-2), (ncol(scc_cps)-1),ncol(scc_cps))]
scc_all <- wranglePeaks(lung_scc_peaks, buildOverlapsData(scc_cps, 4), scc_cps, 48) %>%
  filter(sampleID %in% scc_all_bootstrap$sampleID)

scc_sample_counts = scc_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

scc_all_recurrents <- scc_all %>%
  filter(peak != 0) %>%
  left_join(scc_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped scc changepoints
scc_bootstrap_cps <- buildChangepointData('lung-scc_20bootstrap', blank_genome)
scc_all_bootstrap <- wranglePeaks(lung_scc_peaks, buildOverlapsData(scc_bootstrap_cps, 4),
                                       scc_bootstrap_cps, 60) %>%
  mutate(type='lung-scc')

# look at changepoint confidence levels
hist(scc_bootstrap_cps$bootstrap_support)

ggplot(data = scc_bootstrap_cps, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = scc_cps, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "Changepoint density across all samples", title="Lung-SCC changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
scc_sample_counts_bootstrap = scc_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

scc_all_recurrents_bootstrap <- scc_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(scc_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

scc_recurrent_matches = scc_sample_counts %>%
  left_join(scc_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
scc_to_join = scc_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = scc_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

scc_all_recurrents_final <- scc_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(scc_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# compare changepoints with high confidence
scc_bootstrap_highconf = scc_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

scc_orig_highconf = scc_all %>%
  filter(sampleID %in% scc_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(scc_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(scc_orig_highconf$match, na.rm=TRUE), "of", nrow(scc_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(scc_recurrent_matches$sample_count.y)), "of", nrow(scc_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(scc_all_recurrents_final$match), "of", nrow(scc_all_recurrents_final), "sample-recurrent changepoint matches"))

######### ULYMPH-BNHL #############

# original bnhl recurrent changepoints
bnhl_cps <- buildChangepointData('lymph-bnhl', blank_genome)
bnhl_cps <- bnhl_cps %>%
  dplyr::mutate("type" = 'lymph-bnhl',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
bnhl_cps <- bnhl_cps[,c(1:16, (ncol(bnhl_cps)-2), (ncol(bnhl_cps)-1),ncol(bnhl_cps))]
bnhl_all <- wranglePeaks(lymph_bnhl_peaks, buildOverlapsData(bnhl_cps, 4), bnhl_cps, 106) %>%
  mutate(type='lymph-bnhl') %>%
  filter(sampleID %in% bnhl_all_bootstrap$sampleID)

bnhl_sample_counts = bnhl_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

bnhl_all_recurrents <- bnhl_all %>%
  filter(peak != 0) %>%
  left_join(bnhl_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped bnhl changepoints
bnhl_bootstrap_cps <- buildChangepointData('bnhl_20bootstrap', blank_genome)
bnhl_all_bootstrap <- wranglePeaks(lymph_bnhl_peaks, buildOverlapsData(bnhl_bootstrap_cps, 4),
                                  bnhl_bootstrap_cps, 106) %>%
  mutate(type='lymph-bnhl')

# look at changepoint confidence levels
hist(bnhl_bootstrap_cps$bootstrap_support)

ggplot(data = bnhl_all_bootstrap, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = bnhl_all, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "Changepoint density across all samples", title="Lymph-BNHL changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
bnhl_sample_counts_bootstrap = bnhl_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

bnhl_all_recurrents_bootstrap <- bnhl_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(bnhl_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

bnhl_recurrent_matches = bnhl_sample_counts %>%
  left_join(bnhl_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
bnhl_to_join = bnhl_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = bnhl_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

bnhl_all_recurrents_final <- bnhl_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(bnhl_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# compare changepoints with high confidence
bnhl_bootstrap_highconf = bnhl_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

bnhl_orig_highconf = bnhl_all %>%
  filter(sampleID %in% bnhl_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(bnhl_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(bnhl_orig_highconf$match, na.rm=TRUE), "of", nrow(bnhl_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(bnhl_recurrent_matches$sample_count.y)), "of", nrow(bnhl_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(bnhl_all_recurrents_final$match), "of", nrow(bnhl_all_recurrents_final), "sample-recurrent changepoint matches"))

######### UTERUS-ADENOCA #############

# original uterus recurrent changepoints
uterus_cps <- buildChangepointData('uterus-adenoca', blank_genome)
uterus_cps <- uterus_cps %>%
  dplyr::mutate("type" = 'uterus-adenoca',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
uterus_cps <- uterus_cps[,c(1:16, (ncol(uterus_cps)-2), (ncol(uterus_cps)-1),ncol(uterus_cps))]
uterus_all <- wranglePeaks(uterus_peaks, buildOverlapsData(uterus_cps, 4), uterus_cps, 51) %>%
  mutate(type='uterus-adenoca')

uterus_sample_counts = uterus_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

uterus_all_recurrents <- uterus_all %>%
  filter(peak != 0) %>%
  left_join(uterus_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped uterus changepoints
uterus_bootstrap_cps <- buildChangepointData('uterus-adenoca20bootstrap', blank_genome)
uterus_all_bootstrap <- wranglePeaks(uterus_peaks, buildOverlapsData(uterus_bootstrap_cps, 4),
                                   uterus_bootstrap_cps, 51) %>%
  mutate(type='uterus-adenoca')

# look at changepoint confidence levels
hist(uterus_bootstrap_cps$bootstrap_support)

ggplot(data = uterus_all_bootstrap, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = uterus_all, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "Changepoint density across all samples", title="Uterus-AdenoCA changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
uterus_sample_counts_bootstrap = uterus_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID))
  filter(sample_count >= 7)

uterus_all_recurrents_bootstrap <- uterus_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(uterus_sample_counts, by = c("type", "peak"))
  filter(sample_count >= 7)

uterus_recurrent_matches = uterus_sample_counts %>%
  left_join(uterus_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
uterus_to_join = uterus_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = uterus_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

uterus_all_recurrents_final <- uterus_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(uterus_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

###### COLORECTAL CANCER ######

# original colorect recurrent changepoints
colorect_cps <- buildChangepointData('colorect', blank_genome)
colorect_cps <- colorect_cps %>%
  dplyr::mutate("type" = 'colorect',
                width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
colorect_cps <- colorect_cps[,c(1:16, (ncol(colorect_cps)-2), (ncol(colorect_cps)-1),ncol(colorect_cps))]
colorect_all <- wranglePeaks(colorect_peaks, buildOverlapsData(colorect_cps, 4), colorect_cps, 60) %>%
  mutate(type='colorect') %>%
  filter(sampleID %in% colorect_all_bootstrap$sampleID)

colorect_sample_counts = colorect_all %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
  filter(sample_count >= 7)

colorect_all_recurrents <- colorect_all %>%
  filter(peak != 0) %>%
  left_join(colorect_sample_counts, by = c("type", "peak")) %>%
  filter(sample_count >= 7)

# bootstrapped colorect changepoints
colorect_bootstrap_cps <- buildChangepointData('colorect_20bootstrap', blank_genome)
colorect_all_bootstrap <- wranglePeaks(colorect_peaks, buildOverlapsData(colorect_bootstrap_cps, 4),
                                     colorect_bootstrap_cps, 60) %>%
  mutate(type='colorect-adenoca')

# look at changepoint confidence levels
hist(colorect_bootstrap_cps$bootstrap_support)

ggplot(data = colorect_all_bootstrap, aes(x = mb_bin_start)) +
  geom_density(bw=4, fill="darkgray", alpha=0.7, color=NA) +
  geom_density(data = colorect_all, aes(x = mb_bin_start), fill='red', alpha=0.5, bw=4, color=NA) +
  theme_bw() +
  labs(x = "Genomic Location", y = "Changepoint density across all samples", title="Colorect changepoint density comparison") +
  scale_x_continuous(breaks=change_bins, labels=chr_labels, limits=c(1,3113))

# see how many recurrent changepoints were recovered in the simulations
colorect_sample_counts_bootstrap = colorect_all_bootstrap %>%
  filter(peak != 0) %>%
  group_by(type, peak) %>%
  summarize(sample_count = n_distinct(sampleID)) %>%
filter(sample_count >= 7)

colorect_all_recurrents_bootstrap <- colorect_all_bootstrap %>%
  filter(peak != 0) %>%
  left_join(colorect_sample_counts, by = c("type", "peak"))
  filter(sample_count >= 7)

colorect_recurrent_matches = colorect_sample_counts %>%
  left_join(colorect_sample_counts_bootstrap, by="peak")

# check if the same samples recovered the recurrent changepoint
colorect_to_join = colorect_all_recurrents_bootstrap %>%
  select(peak) %>%
  mutate(bootstrap_sampleID = colorect_all_recurrents_bootstrap$sampleID) %>%
  group_by(bootstrap_sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n)

colorect_all_recurrents_final <- colorect_all_recurrents %>%
  group_by(sampleID, peak) %>%
  summarize(n=n()) %>%
  select(-n) %>%
  left_join(colorect_to_join, by="peak") %>%
  group_by(peak, sampleID) %>%
  summarize(match = any(sampleID == bootstrap_sampleID))

# all TRUE values: all bootstrap samples recover at least the same recurrent
# changepoints as the original samples

# compare changepoints with high confidence
colorect_bootstrap_highconf = colorect_all_bootstrap %>%
  filter(bootstrap_support > 0.05) %>%
  select(start_bin, sampleID)

colorect_orig_highconf = colorect_all %>%
  filter(sampleID %in% colorect_all_bootstrap$sampleID & bootstrap_support > 0.2) %>%
  select(sampleID, start_bin, bootstrap_support) %>%
  left_join(colorect_bootstrap_highconf, by=c("start_bin")) %>%
  group_by(sampleID.x, start_bin) %>%
  summarize(match = any(sampleID.x == sampleID.y))

# print out report
print(paste(sum(colorect_orig_highconf$match, na.rm=TRUE), "of", nrow(colorect_orig_highconf), "high-confidence changepoints recovered"))
print(paste(sum(!is.na(colorect_recurrent_matches$sample_count.y)), "of", nrow(colorect_recurrent_matches), "recurrent changepoints recovered"))
print(paste(sum(colorect_all_recurrents_final$match), "of", nrow(colorect_all_recurrents_final), "sample-recurrent changepoint matches"))
