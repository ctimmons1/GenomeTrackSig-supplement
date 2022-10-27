# GENE DENSITY ASSOCIATION TEST
gene_pvals_all <- data.frame()
types <- c("colorect", "eso-adenoca", "lung-scc", "melanoma", "prost-adenoca", "stomach",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "cns-gbm", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/genedens_pvals/",t,"_genedenspvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  gene_pvals_all <- rbind(gene_pvals_all, pvals)
}

summary_gene <- gene_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_gene <- summary_gene %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),
                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "genedens_windows")

# REPLICATION TIMING ASSOCIATION TEST
rep_pvals_all <- data.frame()
types <- c("colorect", "eso-adenoca", "lung-scc", "melanoma", "prost-adenoca", "stomach",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "cns-gbm", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/replication_windows/",t,"_replicationpvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  rep_pvals_all <- rbind(rep_pvals_all, pvals)
}

summary_rep <- rep_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_rep <- summary_rep %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),
                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "replication_windows")


# KATAEGIS ASSOCIATION TEST
kat_pvals_all <- data.frame()
types <- c("colorect", "eso-adenoca", "lung-scc", "melanoma", "prost-adenoca", "stomach",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "cns-gbm", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/kat_windows_new/",t,"_kataegis_windows_pvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  kat_pvals_all <- rbind(kat_pvals_all, pvals)
}

summary_kat <- kat_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_kat <- summary_kat %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),
                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "kat_windows")

# CNA ASSOCIATION TEST
cna_pvals_all <- data.frame()
types <- c("colorect", "eso-adenoca", "lung-scc", "melanoma", "prost-adenoca", "stomach",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "cns-gbm", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/cna_windows_new/",t,"_cnapvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  cna_pvals_all <- rbind(cna_pvals_all, pvals)
}

summary_cna <- cna_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_cna <- summary_cna %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),

                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "cna_windows")

# A/B COMPARTMENT ASSOCIATION TEST
compartment_pvals_all <- data.frame()
types <- c("colorect", "lung-scc", "prost-adenoca",
           "uterus-adenoca", "bladder", "breast-adenoca",
           "kidney-rcc", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/compartment_windows_new/",t,"_compartmentpvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  compartment_pvals_all <- rbind(compartment_pvals_all, pvals)
}

summary_compartment <- compartment_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_compartment <- summary_compartment %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),

                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "compartment_windows")

# ATAC ASSOCIATION TEST
atac_pvals_all <- data.frame()
types <- c("colorect", "cervix", "cns-gbm", "eso-adenoca", "kidney-rcc-kirc",
          "melanoma", "stomach",  "lung-scc", "prost-adenoca",
           "uterus-adenoca", "bladder", "breast-adenoca",
           "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/atac_pvals/",t,"_atacpvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  atac_pvals_all <- rbind(atac_pvals_all, pvals)
}

summary_atac <- atac_pvals_all %>%
  dplyr::group_by(type, sampleID) %>%
  dplyr::summarize(n_cps = dplyr::n(),
                   n_sig_before = length(null_before[null_before<0.05]),
                   n_sig_after = length(null_after[null_after<0.05]),
                   null_before = mean(null_before),
                   null_after = mean(null_after))

for_table_atac <- summary_atac %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_samples = dplyr::n(),
                   n_samples_signif = dplyr::n_distinct(sampleID[n_sig_before > 0 | n_sig_after > 0]),
                   n_cps = sum(n_cps),
                   n_cps_signif = sum(n_sig_before[n_sig_before>0], n_sig_after[n_sig_after>0])) %>%
  dplyr::mutate("test" = "atac_windows")


results <- rbind(for_table_atac, for_table_cna, for_table_compartment, for_table_kat)


# CNA OVERLAP TEST
cna_overlap_pvals_all <- data.frame()
types <- c("eso-adenoca", "lung-scc", "melanoma", "prost-adenoca",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/cna_overlap_pvals/",t,"_pvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  cna_overlap_pvals_all <- rbind(cna_overlap_pvals_all, pvals)
}

summary_cna_overlap <- cna_overlap_pvals_all %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_sig = n_distinct(id[p_adj_bf < 0.05]))

# TAD OVERLAP TEST
tad_overlap_pvals_all <- data.frame()
types <- c("colorect", "eso-adenoca", "lung-scc", "melanoma", "prost-adenoca", "stomach",
           "uterus-adenoca", "bladder", "breast-adenoca", "cervix", "cns-gbm", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/tad_pvals/",t,"_tadpvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  tad_overlap_pvals_all <- rbind(tad_overlap_pvals_all, pvals)
}

summary_tad <- tad_overlap_pvals_all %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_sig = n_distinct(id[p_adj_bf < 0.05]))

# KATAEGIS OVERLAP TEST
kat_overlap_pvals_all <- data.frame()
types <- c("eso-adenoca", "lung-scc", "prost-adenoca", "stomach",
           "bladder", "breast-adenoca", "cervix", "bone-osteosarc", "kidney-chrcc",
           "kidney-rcc", "lymph-bnhl", "panc-adenoca", "lymph-cll", "lung-adenoca", "head-scc", "thy-adenoca")
for (t in types) {
  pvals <- read.delim(paste0("~/Desktop/CBSP2021/kataegis_overlap/",t,"_kataegis_overlap_pvals.txt"), sep=",") %>%
    dplyr::mutate("type" = t)
  kat_overlap_pvals_all <- rbind(kat_overlap_pvals_all, pvals)
}

summary_kat <- kat_overlap_pvals_all %>%
  dplyr::group_by(type) %>%
  dplyr::summarize(n_sig = n_distinct(id[p_adj_bf < 0.05]))


# DNASE REPORT
dnase_pvals <- read.delim("~/Desktop/CBSP2021/melanoma_dnasepvals.txt", sep=",") %>%
  filter(null_before < 0.05 | null_after < 0.05)
