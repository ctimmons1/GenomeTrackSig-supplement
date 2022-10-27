library(tidyverse)
source("../helper_functions.R")
devtools::load_all()

# make changepoint datasets at different bin sizes
bin200 <- buildChangepointData('lung-scc', blank_genome)
bin100 <- buildChangepointData('lung-scc_100bin', blank_genome)
bin150 <- buildChangepointData('lung-scc_150bin', blank_genome)
bin300 <- buildChangepointData('lung-scc_300bin', blank_genome)

x100 <- c(6,7,10,22,27,30,31,34,39,40,48)
x200 <- c()

# set bin size variable
bin100$binsize=100
bin150$binsize=150
bin200$binsize=200
bin300$binsize=300

# bind all cps together
cps_all <- rbind(bin100,bin150,bin200,bin300) %>%
  dplyr::filter(!sampleID %in% x100)

# make sure the same samples are in all datasets
samples_all <- cps_all %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarize(bins = dplyr::n_distinct(binsize)) %>%
  dplyr::filter(bins == 4)

which_samples <- cps_all %>%
  dplyr::group_by(binsize) %>%
  dplyr::summarize(n_cps = dplyr::n_distinct(sampleID))

bin300 <- bin300 %>%
  dplyr::filter(sampleID %in% samples_all$sampleID)

bin100 <- bin100 %>%
  filter(!sampleID %in% x100)
bin150 <- bin150 %>%
  filter(!sampleID %in% x100)
bin200 <- bin200 %>%
  filter(!sampleID %in% x100)
bin300 <- bin300 %>%
  filter(!sampleID %in% x100)

# find changepoint locations in samples with different bin sizes
changepoint_positions100 <- c()
for (i in 1:nrow(bin100)) {
  changepoint_positions100 <- c(changepoint_positions100, seq(from=bin100$mb_bin_start[i], to=bin100$mb_bin_end[i]))
}
changepoint_positions150 <- c()
for (i in 1:nrow(bin150)) {
  changepoint_positions150 <- c(changepoint_positions150, seq(from=bin150$mb_bin_start[i], to=bin150$mb_bin_end[i]))
}
changepoint_positions200 <- c()
for (i in 1:nrow(bin200)) {
  changepoint_positions200 <- c(changepoint_positions200, seq(from=bin200$mb_bin_start[i], to=bin200$mb_bin_end[i]))
}
changepoint_positions300 <- c()
for (i in 1:nrow(bin300)) {
  changepoint_positions300 <- c(changepoint_positions300, seq(from=bin300$mb_bin_start[i], to=bin300$mb_bin_end[i]))
}

# make changepoint positions into dataframes
changepoint_positions100 <- as.data.frame(changepoint_positions100) %>% dplyr::mutate("Binsize" = "100")
changepoint_positions150 <- as.data.frame(changepoint_positions150) %>% dplyr::mutate("Binsize" = "150")
changepoint_positions200 <- as.data.frame(changepoint_positions200) %>% dplyr::mutate("Binsize" = "200")
changepoint_positions300 <- as.data.frame(changepoint_positions300) %>% dplyr::mutate("Binsize" = "300")

# bind together
changepoint_positions = rbind(changepoint_positions100, changepoint_positions150,
                              changepoint_positions200, changepoint_positions300)

# perform KDE to find changepoint density across the genome at different bin sizes
kde <- ks::kde(changepoint_positions, h=4)
kde_data <- data.frame("eval" = kde$eval.points, "estimate" = kde$estimate)

# plot changepoint density
ggplot2::ggplot() +
  ggplot2::geom_density(data = changepoint_positions100, mapping=ggplot2::aes(x=changepoint_positions100, fill=Binsize), alpha=0.3, adjust=0.1) +
  #ggplot2::geom_density(data = eso_cps, mapping=ggplot2::aes(x=eso_cps, col=binsize), alpha=0.3, adjust=0.1) +
  ggplot2::geom_density(data = changepoint_positions150, mapping=ggplot2::aes(x=changepoint_positions150, fill=Binsize), alpha=0.3, adjust=0.1) +
  ggplot2::geom_density(data = changepoint_positions200, mapping=ggplot2::aes(x=changepoint_positions200, fill=Binsize), alpha=0.3, adjust=0.1) +
  ggplot2::geom_density(data = changepoint_positions300, mapping=ggplot2::aes(x=changepoint_positions300, fill=Binsize), alpha=0.3, adjust=0.1) +
  ggplot2::theme_bw() +
  ggplot2::scale_color_brewer(type='qual', palette="Set1") +
  ggplot2::scale_x_continuous(breaks = change_bins, labels = chr_labels) +
  ggplot2::labs(x = "Chromosome", y = "Changepoint Density", col = 'Bin size') +
  ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank())
ggplot2::ggsave(filename="binsize-comparison.png", plot=ggplot2::last_plot(), device="png", path="~/Desktop/CBSP2021/plots",
                scale=1, width = 35, height = 10, units = "cm", dpi=600)

