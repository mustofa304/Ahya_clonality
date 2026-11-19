setwd("/Users/mustofa/Desktop/ibsdist211bams")

library(vcfR)
library(adegenet)
library(poppr)
library(dplyr)
library(ggplot2)
library(hierfstat)

# -------------------- INPUT FILES --------------------

vcf_file <- "ahya_bams211_noout.bcf"
ibs_file <- "ahya_bams211_noout.ibsMat"
metadata_file <- "bams2nooutliers.txt"
coverage_file <- "/Users/mustofa/Desktop/Ahya_rerun/qranksAhya.txt"

# IBS distance cutoff for defining MLGs/MLLs
# Change this if your final IBS cutoff is different.
ibs_threshold <- 0.12

# -------------------- READ VCF AND CREATE GENCLONE --------------------

vcf <- read.vcfR(vcf_file, convertNA = TRUE)
genind <- vcfR2genind(vcf)
genclone <- as.genclone(genind)
genclone
sample_names <- indNames(genclone)

# -------------------- READ SAMPLE AND SITE METADATA --------------------
#
# bams2nooutliers.txt has no header:
# column 1 = sample name
# column 2 = site/population name

sample_info <- read.table(
  metadata_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

if (ncol(sample_info) < 2) {
  stop("bams2nooutliers.txt must contain sample names in column 1 and sites in column 2.")
}

sample_info <- sample_info[, 1:2]
colnames(sample_info) <- c("sample", "POP")

sample_info$sample <- basename(as.character(sample_info$sample))
sample_info$POP <- as.character(sample_info$POP)

if (anyDuplicated(sample_info$sample)) {
  stop("Duplicate sample names found in column 1 of bams2nooutliers.txt.")
}

if (!all(sample_names %in% sample_info$sample)) {
  missing_metadata <- setdiff(sample_names, sample_info$sample)
  stop(
    "These VCF samples are missing from bams2nooutliers.txt: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# Arrange metadata in the same order as the VCF samples.
sample_info <- sample_info[
  match(sample_names, sample_info$sample),
  ,
  drop = FALSE
]

rownames(sample_info) <- sample_info$sample

# Add site and replicate information to genclone.
strata_df <- data.frame(
  IND = sample_names,
  POP = sample_info$POP,
  SAMPLE_ID = basename(sample_names),
  REPLICATE = ifelse(
    grepl("redo\\.bam$", basename(sample_names), ignore.case = TRUE),
    "Technical replicate",
    "Original"
  ),
  stringsAsFactors = FALSE,
  row.names = sample_names
)

strata(genclone) <- strata_df
setPop(genclone) <- ~POP

head(strata(genclone))
table(strata(genclone)$POP)
table(strata(genclone)$REPLICATE)

# -------------------- READ ANGSD IBS MATRIX --------------------
#
# ANGSD .ibsMat is a numeric distance matrix without sample names.
# Its row and column order must match the sample order used for ANGSD.
# Here that order is taken from column 1 of bams2nooutliers.txt.

bam_order <- sample_info$sample

ibs_raw <- read.table(
  ibs_file,
  header = FALSE,
  check.names = FALSE,
  fill = TRUE,
  comment.char = "",
  stringsAsFactors = FALSE
)

# Remove only columns that are completely empty, if present.
ibs_raw <- ibs_raw[, colSums(!is.na(ibs_raw)) > 0, drop = FALSE]

ibs_mat <- as.matrix(ibs_raw)
storage.mode(ibs_mat) <- "numeric"

if (nrow(ibs_mat) != length(bam_order) ||
    ncol(ibs_mat) != length(bam_order)) {
  stop(
    "The IBS matrix dimensions do not match the number of samples in ",
    metadata_file,
    ". Check that the metadata file contains the same samples, in the same ",
    "order, as the BAM list used for ANGSD."
  )
}

# Assign ANGSD sample order to the matrix.
rownames(ibs_mat) <- bam_order
colnames(ibs_mat) <- bam_order

# Reorder matrix to match VCF/genclone sample order.
ibs_mat <- ibs_mat[sample_names, sample_names, drop = FALSE]

if (!isTRUE(all.equal(ibs_mat, t(ibs_mat), tolerance = 1e-8))) {
  stop("The IBS matrix is not symmetric. Check the input matrix.")
}

diag(ibs_mat) <- 0
ibs_dist <- as.dist(ibs_mat)

# -------------------- CHECK IBS DISTANCES AND ESTIMATE CUTOFFS --------------------

hist(
  ibs_dist,
  main = "Pairwise ANGSD IBS distances",
  xlab = "IBS distance"
)

gac <- genotype_curve(
  genclone,
  sample = 100,
  maxloci = 100,
  quiet = TRUE
)

genclone_filtered <- filter_stats(
  genclone,
  distance = ibs_dist,
  plot = TRUE
)

farthest_thresh <- cutoff_predictor(
  genclone_filtered$farthest$THRESHOLDS
)
average_thresh <- cutoff_predictor(
  genclone_filtered$average$THRESHOLDS
)
nearest_thresh <- cutoff_predictor(
  genclone_filtered$nearest$THRESHOLDS
)

cat("Farthest-neighbor cutoff estimate:", farthest_thresh, "\n")
cat("Average-neighbor cutoff estimate:", average_thresh, "\n")
cat("Nearest-neighbor cutoff estimate:", nearest_thresh, "\n")
cat("IBS threshold used below:", ibs_threshold, "\n")

# -------------------- TECHNICAL REPLICATE DISTANCES --------------------

TR2s <- sample_names[
  grepl("redo\\.bam$", basename(sample_names), ignore.case = TRUE)
]

TR1s <- sub(
  "redo(?=\\.bam$)",
  "",
  TR2s,
  perl = TRUE,
  ignore.case = TRUE
)

complete_pairs <- TR1s %in% sample_names
TR1s <- TR1s[complete_pairs]
TR2s <- TR2s[complete_pairs]
TRs <- c(TR1s, TR2s)

TR_pairs <- data.frame(
  original = TR1s,
  technical_replicate = TR2s,
  stringsAsFactors = FALSE
)

print(TR_pairs)

ibs_matrix <- as.matrix(ibs_dist)

if (length(TR1s) > 0) {
  TRpw_dist <- mapply(
    FUN = function(original, redo) ibs_matrix[original, redo],
    original = TR1s,
    redo = TR2s
  )
  
  names(TRpw_dist) <- TR2s
  print(TRpw_dist)
  
  cat("Mean technical replicate distance:", mean(TRpw_dist, na.rm = TRUE), "\n")
  cat("Median technical replicate distance:", median(TRpw_dist, na.rm = TRUE), "\n")
  cat("Maximum technical replicate distance:", max(TRpw_dist, na.rm = TRUE), "\n")
  cat("Minimum technical replicate distance:", min(TRpw_dist, na.rm = TRUE), "\n")
} else {
  TRpw_dist <- numeric(0)
  message("No complete original/redo technical replicate pairs were found.")
}

# -------------------- ASSIGN MLGs/MLLs USING IBS DISTANCE --------------------
#
# Algorithm "f" is the farthest-neighbor algorithm.

mlg.filter(
  genclone,
  distance = ibs_dist,
  algorithm = "f"
) <- ibs_threshold

cat("Number of MLLs at IBS threshold", ibs_threshold, ": ",
    nmll(genclone), "\n", sep = "")

MLLtab <- data.frame(
  sample = indNames(genclone),
  MLL = as.character(mll(genclone)),
  stringsAsFactors = FALSE
)

# Inspect technical replicate assignments.
MLLtab[MLLtab$sample %in% TRs, ]

write.csv(
  MLLtab,
  "MLL_assignments_IBS_pt12.csv",
  row.names = FALSE
)

# -------------------- IBS DISTANCE HISTOGRAM --------------------

ibs_values <- as.vector(ibs_dist)

ibs_hist_data <- data.frame(IBS_distance = ibs_values)

ibs_hist <- ggplot(ibs_hist_data, aes(x = IBS_distance)) +
  geom_histogram(
    color = "black",
    fill = "lightgray",
    binwidth = 0.005
  ) +
  geom_vline(
    xintercept = ibs_threshold,
    color = "blue",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = farthest_thresh,
    color = "red",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = TRpw_dist,
    linetype = "dashed"
  ) +
  coord_cartesian(xlim = c(0, 0.40)) +
  theme_classic(base_size = 14) +
  labs(
    x = "ANGSD IBS distance",
    y = "Number of pairwise comparisons"
  )

print(ibs_hist)

ggsave(
  filename = "IBS_distance_histogram_pt12.png",
  plot = ibs_hist,
  width = 8,
  height = 5,
  dpi = 300
)

# -------------------- IBS DENDROGRAM --------------------

hc <- hclust(ibs_dist, method = "complete")

plot(
  hc,
  labels = FALSE,
  main = "Hierarchical clustering using ANGSD IBS distance",
  ylab = "IBS distance"
)

abline(h = farthest_thresh, col = "red", lty = 2)
abline(h = ibs_threshold, col = "blue", lty = 2)

if (length(TRpw_dist) > 0) {
  abline(h = TRpw_dist, col = "black", lty = 3)
}

# -------------------- POPPR DIVERSITY: ALL ORIGINAL SAMPLES --------------------
#
# Remove redo files, but retain original samples, including repeated MLLs.
# The IBS filter on genclone determines which samples share an MLL.

keep_original <- sample_names[!sample_names %in% TR2s]
gensub <- genclone[keep_original, ]
gensub

population_order <- c(
  "GBH", "GCH2", "SKH",
  "BR1", "SKR", "SBR",
  "OR1", "OR2", "OR3"
)

pop(gensub) <- factor(
  as.character(pop(gensub)),
  levels = population_order
)

poppr_results <- poppr(gensub, sample = 189)
print(poppr_results)

# Check whether MLLs occur in more than one population.
mlg.crosspop(gensub, strata = ~POP)
mlg.table(gensub, strata = ~POP)

poppr_results$R <- ifelse(
  poppr_results$N > 1,
  (poppr_results$MLG - 1) / (poppr_results$N - 1),
  NA_real_
)

print(poppr_results[, c("Pop", "N", "MLG", "R")])

write.csv(
  poppr_results,
  "poppr_IBS_diversity_211bamspt12.csv",
  row.names = FALSE
)

# -------------------- SHARED MLL TABLE --------------------

mlg_data <- data.frame(
  sample = indNames(gensub),
  POP = as.character(pop(gensub)),
  MLL = as.character(mll(gensub)),
  stringsAsFactors = FALSE
)

shared_MLLs <- mlg_data %>%
  group_by(MLL) %>%
  summarise(
    number_of_samples = n(),
    number_of_populations = n_distinct(POP),
    populations = paste(sort(unique(POP)), collapse = ", "),
    samples = paste(sample, collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(number_of_populations > 1) %>%
  arrange(desc(number_of_populations), MLL)

print(shared_MLLs)

write.csv(
  shared_MLLs,
  "shared_MLLs_IBS.csv",
  row.names = FALSE
)

# -------------------- MLL FREQUENCY PLOT --------------------

mlg_plot_data <- mlg_data %>%
  mutate(
    POP = factor(POP, levels = population_order),
    MLL = paste0("MLL.", MLL)
  )

mlg_counts <- mlg_plot_data %>%
  count(POP, MLL, name = "number_of_samples")

mlg_plot <- ggplot(mlg_counts, aes(x = MLL, y = number_of_samples)) +
  geom_col(fill = "gray40", color = "white", linewidth = 0.25) +
  facet_wrap(~POP, ncol = 3, scales = "free_x", drop = TRUE) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "MLL frequencies by site",
    subtitle = paste0(
      "N = ", nrow(mlg_data),
      "; unique MLLs = ", n_distinct(mlg_data$MLL),
      "; IBS threshold = ", ibs_threshold
    ),
    x = "Multilocus lineage (MLL)",
    y = "Number of samples"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "gray85", color = "gray50"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(mlg_plot)

ggsave(
  filename = "MLL_frequency_IBS_211bams_pt12.png",
  plot = mlg_plot,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300
)

# -------------------- SELECT ONE REPRESENTATIVE PER SITE AND MLL --------------------
#
# If the same MLL occurs at multiple sites, this keeps one representative
# within each site. This preserves that MLL's occurrence at each site.

cov <- read.table(
  coverage_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(cov) < 2) {
  stop("Coverage file must have sample names in column 1 and coverage in column 2.")
}

cov <- cov[, 1:2]
colnames(cov) <- c("sample", "sites5x")
cov$sample <- basename(as.character(cov$sample))
cov$sites5x <- as.numeric(cov$sites5x)

MLLwCov <- MLLtab %>%
  left_join(sample_info, by = "sample") %>%
  left_join(cov, by = "sample")

if (any(is.na(MLLwCov$sites5x))) {
  warning("Some samples have no matching coverage value; check the coverage file.")
}

samples4popgen <- MLLwCov %>%
  filter(!is.na(sites5x)) %>%
  group_by(POP, MLL) %>%
  slice_max(order_by = sites5x, n = 1, with_ties = FALSE) %>%
  ungroup()

writeLines(
  samples4popgen$sample,
  "representative_samples_per_site_MLL_IBS.txt"
)

write.csv(
  samples4popgen,
  "representative_samples_per_site_MLL_IBS.csv",
  row.names = FALSE
)


# -------------------- HIERFSTAT Ho, He, AND Fis --------------------
#
# Export site-level genotype diversity for:
# 1) all original samples (redo files excluded)
# 2) one coverage-selected representative per site and MLL
#
# These values are calculated from VCF genotypes, not from the IBS matrix.

summarize_by_site <- function(x, stat_name) {
  x <- as.matrix(x)
  
  n_loci <- colSums(!is.na(x))
  means <- colMeans(x, na.rm = TRUE)
  
  se <- vapply(seq_len(ncol(x)), function(j) {
    if (n_loci[j] < 2) {
      NA_real_
    } else {
      stats::sd(x[, j], na.rm = TRUE) / sqrt(n_loci[j])
    }
  }, numeric(1))
  
  data.frame(
    Site = colnames(x),
    setNames(data.frame(means), stat_name),
    setNames(data.frame(se), paste0(stat_name, "_SE")),
    setNames(data.frame(n_loci), paste0("n_loci_", stat_name)),
    row.names = NULL
  )
}

make_hierfstat_summary <- function(dat) {
  basicstat <- basic.stats(dat, diploid = TRUE, digits = 10)
  
  Ho_summary <- summarize_by_site(basicstat$Ho, "Ho")
  He_summary <- summarize_by_site(basicstat$Hs, "He")
  Fis_summary <- summarize_by_site(basicstat$Fis, "Fis")
  
  Reduce(
    function(x, y) merge(x, y, by = "Site", all = TRUE, sort = FALSE),
    list(Ho_summary, He_summary, Fis_summary)
  )
}

# All original samples, with redo files removed
hierfstat_all_original <- make_hierfstat_summary(gensub)
print(hierfstat_all_original)

write.csv(
  hierfstat_all_original,
  "hierfstat_IBS_211bamspt12.csv",
  row.names = FALSE
)

# Clone-corrected representatives, one per site and MLL
hierfstat_representatives <- make_hierfstat_summary(
  genclone_representatives
)
print(hierfstat_representatives)

write.csv(
  hierfstat_representatives,
  "hierfstat_IBS_211bamspt12_with_SE.csv",
  row.names = FALSE
)
