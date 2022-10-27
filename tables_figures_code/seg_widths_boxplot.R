# build dataset of all changepoints across all cancer types
types <- c('melanoma', 'lung-scc', 'eso-adenoca', 'colorect', 'uterus-adenoca',
           'lung-adenoca', 'kidney-rcc', 'kidney-chrcc', 'cns-gbm', 'breast-adenoca',
           'prost-adenoca', 'bladder', 'head-scc', 'panc-adenoca', 'lymph-bnhl',
           'lymph-cll', 'cervix', 'bone-osteosarc', 'thy-adenoca', 'stomach')

all_cps <- data.frame()
for (i in 1:length(types)) {
  cps <- buildChangepointData(types[i], blank_genome)
  cps <- cps %>%
    dplyr::mutate("type" = types[i],
                  width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
  cps <- cps[,c(1:17, (ncol(cps)-2), (ncol(cps)-1),ncol(cps))]
  all_cps <- rbind(all_cps,  cps)
}

# build dataset of all segments across all samples
segment_widths <- data.frame()
for (type in types) {
  for (i in c(1:330)) {
    if (file.exists(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))) {
      print(i)
      load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
      temp <- load(paste0("~/Desktop/CBSP2021/",type,"/",type,i,".Rdata"))
      binData <- get(temp)[['binData']]
      binData <- binData[,1:102] %>%
        mutate("type" = type,
               "sampleID" = i)
      segment_widths <- rbind(segment_widths, binData)
    }
  }
}

# find geometric mean tmb per cancer type
n_segs <- segment_widths %>%
  group_by(type, sampleID) %>%
  summarize(n_segs = n())
n_segs$tmb = tmb$tmb
geom_mean_tmb <- n_segs %>%
  group_by(type) %>%
  summarize(mean_tmb = prod(tmb)^(1/n()))
tmb <- segment_widths %>%
  group_by(type, sampleID) %>%
  summarize_at(vars(C_A_ACA:T_G_TTT), sum)
tmb$tmb <- rowSums(tmb[,3:98])/3113
tmb <- tmb %>%
  group_by(type) %>%
  summarize(median_tmb = median(tmb)/3113)

# exploratory data analysis
all_cps <- all_cps %>%
  mutate(width = if_else(mb_bin_end-mb_bin_start == 0, 1, mb_bin_end-mb_bin_start))
summary <- all_cps %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(mean_width = mean(width),
                   median_width = median(width),
                   q25_width = quantile(width, 0.25),
                   q75_width = quantile(width, 0.75),
                   min_width = min(width),
                   max_width = max(width)) %>%
  dplyr::mutate(iqr = q75_width - q25_width)

segment_summary <- segment_widths %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(mean_width = mean(width)/1e6,
                   median_width = median(width)/1e6,
                   q25_width = quantile(width, 0.25)/1e6,
                   q75_width = quantile(width, 0.75)/1e6,
                   min_width = min(width)/1e6,
                   max_width = max(width)/1e6) %>%
  dplyr::mutate(iqr = q75_width - q25_width)

library(tidyverse)
all_cps$type <- factor(all_cps$type, levels=c("colorect", "uterus-adenoca","melanoma", "lung-scc", "lung-adenoca", "stomach",
                                              "cns-gbm", "lymph-bnhl", "eso-adenoca", "bladder", "breast-adenoca", 'kidney-rcc',
                                              "head-scc", "cervix", "bone-osteosarc", "prost-adenoca", "panc-adenoca",
                                              "lymph-cll", "kidney-chrcc", "thy-adenoca"))
ggplot(all_cps, aes(x = width, y = type, fill=type)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_continuous(limits=c(0,250)) +
  labs(x = "Changepoint Segment Width (Mb)", y = "Cancer Type") +
  theme(legend.position="none")

# make dataset for plotting
all_cps2 <- all_cps %>%
  dplyr::select(type, width) %>%
  mutate(segment = "Changepoint Bins")
all_segments2 <- segment_widths %>%
  dplyr::select(type, width, segment) %>%
  mutate(width=width/1e6,
         segment="All Bins")

plotting_data <- rbind(all_cps2, all_segments2) %>%
  arrange(type) %>%
  mutate(
    type = case_when(type=="melanoma" ~ "Melanoma (N = 107, TMB = 16.68)",
                     type=="lung-scc" ~ "Lung-SCC (N = 48, TMB = 12.81)",
                     type=="eso-adenoca" ~ "Eso-AdenoCA (N = 97, TMB = 7.12)",
                     type=="colorect" ~ "Colorect-AdenoCA (N = 60, TMB = 9.65)",
                     type=="bladder" ~ "Bladder-TCC (N = 23, TMB = 5.69)",
                     type=="stomach" ~ "Stomach-AdenoCA (N = 67, TMB = 4.19)",
                     type=="head-scc" ~ "Head-SCC (N = 54, TMB = 3.31)",
                     type=="lymph-bnhl" ~ "Lymph-BNHL (N = 106, TMB = 2.54)",
                     type=="uterus-adenoca" ~ "Uterus-AdenoCA (N = 51, TMB = 5.91)",
                     type=="cns-gbm" ~ "CNS-GBM (N = 41, TMB = 2.30)",
                     type=="lung-adenoca" ~ "Lung-AdenoCA (N = 33, TMB = 5.11)",
                     type=="breast-adenoca" ~ "Breast-AdenoCA (N =193, TMB = 1.64)",
                     type=="kidney-rcc" ~ "Kidney-RCC (N = 144, TMB = 1.79)",
                     type=="cervix" ~ "Cervix (N = 20, TMB = 1.79)",
                     type=="bone-osteosarc" ~ "Bone-Osteosarc (N = 39, TMB = 1.16)",
                     type=="prost-adenoca" ~ "Prost-AdenoCA (N = 145, TMB = 1.06)",
                     type=="panc-adenoca" ~ "Panc-AdenoCA (N = 238, TMB = 1.78)",
                     type=="lymph-cll" ~ "Lymph-CLL (N = 95, TMB = 0.73)",
                     type=="kidney-chrcc" ~ "Kidney-ChRCC (N = 38, TMB = 0.64)",
                     type=="thy-adenoca" ~ "Thy-AdenoCA (N = 29, TMB = 0.57)"))
plotting_data$type <- factor(plotting_data$type, levels=c("Melanoma (N = 107, TMB = 16.68)", "Lung-SCC (N = 48, TMB = 12.81)",
                                                          "Colorect-AdenoCA (N = 60, TMB = 9.65)", "Eso-AdenoCA (N = 97, TMB = 7.12)",
                                                          "Uterus-AdenoCA (N = 51, TMB = 5.91)", "Bladder-TCC (N = 23, TMB = 5.69)",
                                                          "Lung-AdenoCA (N = 33, TMB = 5.11)", "Stomach-AdenoCA (N = 67, TMB = 4.19)",
                                                          "Head-SCC (N = 54, TMB = 3.31)", "Lymph-BNHL (N = 106, TMB = 2.54)",
                                                          "CNS-GBM (N = 41, TMB = 2.30)", "Kidney-RCC (N = 144, TMB = 1.79)",
                                                          "Cervix (N = 20, TMB = 1.79)", "Panc-AdenoCA (N = 238, TMB = 1.78)",
                                                          "Breast-AdenoCA (N =193, TMB = 1.64)", "Bone-Osteosarc (N = 39, TMB = 1.16)",
                                                          "Prost-AdenoCA (N = 145, TMB = 1.06)","Lymph-CLL (N = 95, TMB = 0.73)",
                                                          "Kidney-ChRCC (N = 38, TMB = 0.64)", "Thy-AdenoCA (N = 29, TMB = 0.57)"
))

plotting_data <- plotting_data %>%
  arrange(type)
lengths <- plotting_data %>%
  group_by(type) %>%
  summarize(n = n()) %>%
  .$n
tmbs <- rep(geom_mean_tmb$mean_tmb, times=lengths)
plotting_data$tmb <- tmbs
plotting_data$type <- factor(plotting_data$type, levels=c("Colorect-AdenoCA (N = 60, TMB = 9.65)", "Uterus-AdenoCA (N = 51, TMB = 5.91)",
                                                          "Melanoma (N = 107, TMB = 16.68)", "Lung-SCC (N = 48, TMB = 12.81)",
                                                          "Lung-AdenoCA (N = 33, TMB = 5.11)", "Stomach-AdenoCA (N = 67, TMB = 4.19)",
                                                          "CNS-GBM (N = 41, TMB = 2.30)", "Lymph-BNHL (N = 106, TMB = 2.54)",
                                                          "Eso-AdenoCA (N = 97, TMB = 7.12)", "Bladder-TCC (N = 23, TMB = 5.69)",
                                                          "Breast-AdenoCA (N =193, TMB = 1.64)", "Kidney-RCC (N = 144, TMB = 1.79)",
                                                          "Head-SCC (N = 54, TMB = 3.31)", "Cervix (N = 20, TMB = 1.79)",
                                                          "Bone-Osteosarc (N = 39, TMB = 1.16)",
                                                          "Prost-AdenoCA (N = 145, TMB = 1.06)", "Panc-AdenoCA (N = 238, TMB = 1.78)", "Lymph-CLL (N = 95, TMB = 0.73)",
                                                          "Kidney-ChRCC (N = 38, TMB = 0.64)", "Thy-AdenoCA (N = 29, TMB = 0.57)"
))
plotting_data$segment <- factor(plotting_data$segment, levels=c("All Bins", "Changepoint Bins"))
vals <- scales::rescale(plotting_data$tmb, to=c(0,1))

# make plot
ggplot(plotting_data, aes(x = width, y = type, col=segment)) +
  geom_boxplot(position=position_dodge(), aes(fill=tmb)) +
  scale_color_manual(values=c("black", "saddlebrown")) +
  scale_fill_distiller(palette="YlOrRd", direction=-1, values=vals) +
  theme_bw() +
  labs(x = "Bin Width (Mb)", y = "Cancer Type", fill = "TMB (Geometric Avg.)",
       color="Bin Type") +
  theme(axis.text.x = ggplot2::element_text(size=12),
        axis.text.y = ggplot2::element_text(size=12),
        axis.title.x = ggplot2::element_text(size=16),
        axis.title.y = ggplot2::element_text(size=16),
        legend.title = ggplot2::element_text(size=16),
        legend.text = ggplot2::element_text(size=12))

ggplot2::ggsave(filename=paste0("seg_widths_revised.pdf"), path="~/Desktop/CBSP2021/plots", plot=ggplot2::last_plot(), device="pdf",
                scale=1, width=40, height=30, unit="cm", dpi=600)
