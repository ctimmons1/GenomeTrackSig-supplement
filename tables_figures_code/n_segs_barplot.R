# build dataset of all segments across all cancer samples
types <- c('melanoma', 'lung-scc', 'eso-adenoca', 'colorect', 'uterus-adenoca',
           'lung-adenoca', 'kidney-rcc', 'kidney-chrcc', 'cns-gbm', 'breast-adenoca',
           'prost-adenoca', 'bladder', 'head-scc', 'panc-adenoca', 'lymph-bnhl',
           'lymph-cll', 'cervix', 'bone-osteosarc', 'thy-adenoca', 'stomach')
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
# find number of segments per sample
n_segs <- segment_widths %>%
  group_by(type, sampleID) %>%
  summarize(n_segs = n()) %>%
  mutate(n_bin = case_when(n_segs <=50 ~ 1,
                           n_segs > 50 & n_segs <=100 ~ 2,
                           n_segs >100 & n_segs <= 200 ~ 3,
                           n_segs > 200 & n_segs <=400 ~ 4,
                           n_segs >400 & n_segs <= 800 ~ 5,
                           n_segs > 800 ~ 6),
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

# find geometric mean tmb per cancer type
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

# plot number of segments per sample
breaks <- c(1,2,3,4,5,6)
labels <- c("23-50","51-100","101-200","201-400","401-800",">800")
n_segs$type <- factor(n_segs$type, levels=c("Melanoma (N = 107, TMB = 16.68)", "Lung-SCC (N = 48, TMB = 12.81)",
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
ggplot(data = n_segs, aes(x = n_bin, fill=type)) +
  geom_bar(position="stack") +
  theme_bw() +
  scale_x_continuous(breaks = breaks, labels=labels) +
  scale_y_continuous(breaks=c(0,100,200,300,400,500,600,700)) +
  ggsci::scale_fill_ucscgb() +
  labs(x = "N Bins", y = "N Samples", fill = "Cancer Type")
ggsave(filename="n_bins_vs_samples.pdf", path="~/Desktop/CBSP2021/plots", plot=ggplot2::last_plot(), device="pdf",
       scale=1, width=30, height=20, unit="cm", dpi=600)

segment_widths$segment <- "all"
segment_widths$type <- factor(segment_widths$type, levels=c("colorect", "uterus-adenoca","melanoma", "lung-scc", "lung-adenoca", "stomach",
                                                            "cns-gbm", "lymph-bnhl", "eso-adenoca", "bladder", "breast-adenoca", 'kidney-rcc',
                                                            "head-scc", "cervix", "bone-osteosarc", "prost-adenoca", "panc-adenoca",
                                                            "lymph-cll", "kidney-chrcc", "thy-adenoca"))
ggplot(segment_widths, aes(x = width/1e6, y = type, fill=type)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Segment Width (Mb)", y = "Cancer Type") +
  theme(legend.position="none")
