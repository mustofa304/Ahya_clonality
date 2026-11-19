##############
## UPDATED ##
#############

library(dplyr)
library(tidyverse)
library(janitor)
library(ggdendro)

setwd("/Users/mustofa/Desktop/ibsdist211bams")

raw <- read_tsv("ahya_bams211_noout.ibsMat",col_names = FALSE,show_col_types = FALSE) %>%select(-X212)
samples <- read_tsv("bams2nooutliers.txt",col_names = FALSE,show_col_types = FALSE) %>%select(-X2) %>%rename(sample = X1)
meta <- read_tsv("bams2nooutliers_meta.txt") %>%clean_names() %>%mutate(sample = as.factor(sample))

colnames(raw) <- as.vector(samples$sample)
rownames(raw) <- samples$sample

nrow(raw)
ncol(raw)
nrow(samples)

###################### DETERMINE THRESHOLD FROM REPLICATES ######################

replicate_names <- colnames(raw)[grep("redo", colnames(raw))]
replicate_names

line_numbers <- which(colnames(raw) %in% replicate_names)
line_numbers

indices <- c(8,26,36,172,174,199,208,
             7,25,35,171,173,198,207)

replicates <- raw[indices, indices]

replicate_matrix <- as.matrix(replicates)
replicate_matrix[upper.tri(replicate_matrix)] <- NA
diag(replicate_matrix) <- NA

replicate_distances <- as.data.frame(replicate_matrix) %>%
  gather(sample, dist) %>%
  drop_na()

replist <- as.data.frame(colnames(replicates))

range(replicate_distances$dist)

ggplot(replicate_distances) +geom_histogram(aes(dist))

cut <- quantile(replicate_distances %>%filter(dist < 0.125) %>% pull(dist),c(.95))
cut

plotdata <- replicate_distances %>%filter(dist < 0.32) %>%mutate(group = case_when(dist > 0.115 ~ "Outlier",TRUE ~ "yes"))

plotdata

dist <- ggplot(plotdata, aes(x = dist)) +
  geom_density(aes(y = after_stat(density / 5), fill = group), alpha = 0.25) +
  scale_fill_manual(values = c("red", "gray")) +
  theme_classic(base_size = 17) +
  xlab("Genetic Distance (1-IBS)") +
  theme(legend.position = "none",axis.title.x = element_blank()) +
  scale_x_continuous(limits = c(0.00, 0.42),breaks = seq(0.05, 0.42, 0.05), expand = expansion(mult = c(0,0))) +
  coord_flip() +
  annotate("text", x = 0.25, y = 3.2, label = "outliers",size = 4, color = "red", fontface = "italic") +
  annotate("text", x = 0.14, y = 3.2, label = "threshold",size = 4, color = "gray", fontface = "italic") +
  geom_vline(xintercept = cut,linetype = "dashed",color = "gray")

dist

ggsave(filename = "/Users/mustofa/Desktop/ibsdist211bams/densityplot_thres_pt12.png",
       plot = dist,
       width = 3,
       height = 5,
       dpi = 300)

###################### CLUSTER SAMPLES ######################

ibs_matrix <- as.matrix(raw)
ibs_matrix[upper.tri(ibs_matrix)] <- NA
diag(ibs_matrix) <- NA

cluster <- hclust(as.dist(ibs_matrix), method = "complete")
plot(cluster)

clones <- rect.hclust(cluster, h = 0.12)
length(clones)

cloneslist <- do.call(rbind,lapply(seq_along(clones), function(i) {data.frame(CLUSTER = i, clones[[i]])}))

###################### SUMMARIZE CLONALITY ######################

summary_clonality <- cloneslist %>%
  rownames_to_column(var = "sample") %>%
  select(-3) %>%
  left_join(meta, by = "sample") %>%
  select(sample, CLUSTER, site) %>%
  filter(!str_detect(sample, "redo")) %>%
  mutate(CLUSTER = as.factor(CLUSTER)) %>%
  group_by(CLUSTER) %>%
  mutate(sites = n_distinct(site))

summary <- cloneslist %>%
  rownames_to_column(var = "sample") %>%
  select(-3) %>%
  clean_names() %>%
  inner_join(meta, by = "sample") %>%
  filter(!str_detect(sample, "redo")) %>%
  group_by(site) %>%
  summarise(
    clusters = n_distinct(cluster),
    samples = n_distinct(sample),
    .groups = "drop") %>%
  rename(genotypes = clusters) %>%
  mutate(G_R = genotypes / samples) %>%
  separate(site, into = c("block", "trash"), sep = "_", remove = FALSE)

clustercounts <- summary_clonality %>%
  select(CLUSTER) %>%
  group_by(CLUSTER) %>%
  tally()

summary_clonality %>%
  filter(CLUSTER == 211)

###################### CLASSIFY SAMPLES ######################

working <- cloneslist %>%
  rownames_to_column(var = "sample") %>%
  left_join(cloneslist %>%group_by(CLUSTER) %>%tally() %>%filter(n > 1),by = "CLUSTER") %>%
  mutate(group = case_when(grepl("redo|^BR1-71|^GBH-146|^GBH-157|^GBH-71|^SKH-110|^SKH-111|^SKR-14|^SKR-27",sample) ~ "Technical Replicates",
                           n > 1 ~ "Inferred Clones",TRUE ~ "Sample")) %>%
  rename(label = sample) %>%
  mutate(group = case_when(label == "211" ~ "Sample",TRUE ~ as.character(group)))

###################### PREPARE DENDROGRAM ######################

hcdata <- dendro_data(cluster, type = "rectangle")

# Match each leaf label to its group
label_group <- hcdata$labels %>%
  left_join( working %>% select(label, group),by = "label")

# Get x positions for technical replicate and inferred clone tips
technical_replicate_x <- label_group %>%
  filter(group == "Technical Replicates") %>%
  pull(x)

inferred_clone_x <- label_group %>%
  filter(group == "Inferred Clones") %>%
  pull(x)

# Color technical replicate tips red, inferred clone tips steel blue,
# and all other dendrogram lines black
segdata <- hcdata$segments %>%
  mutate(line_color = case_when(
    x == xend & yend == 0 & xend %in% technical_replicate_x ~ "red",
    x == xend & yend == 0 & xend %in% inferred_clone_x ~ "steelblue",
    TRUE ~ "black"),
    yend = case_when(
      yend == 0 ~ y - 0.01,
      TRUE ~ yend))

###################### PLOT DENDROGRAM ######################

tree <- ggplot() +
  geom_hline(yintercept = 0.12,linetype = "dashed",color = "gray") +
  geom_segment(data = segdata,
               aes(x = x, y = y,
                   xend = xend, yend = yend,
                   color = line_color),size = 0.5) +
  scale_color_identity() +
  geom_text(data = hcdata$labels,aes(x = x, y = y, label = label),angle = 90, hjust = 1,size = 0) +
  theme_classic(base_size = 17) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.position = c(0.86, 0.86)) +
  scale_color_identity(
    name = NULL,
    breaks = c("red", "steelblue"),
    labels = c("Technical Replicates", "Inferred Clones"),
    guide = "legend") +
  ylab("Genetic Distance (1-IBS)") +
  scale_y_continuous(
    limits = c(0, 0.42),
    breaks = seq(0.05, 0.4, 0.05),
    expand = expansion(mult = c(0, 0)))

tree

ggsave(filename = "/Users/Mustofa/Desktop/ibsdist211bams/IBS_bam211pt125.png",
       plot = tree,
       width = 9,
       height = 5,
       dpi = 300)

##TEST - Coloured the dendrogram###

hcdata <- dendro_data(
  hclust(as.dist(matrix), method = "complete"),
  type = "rectangle"
)

# Identify terminal branches and position labels at their tips
tip_labels <- hcdata$segments %>%
  filter(yend == 0, x == xend) %>%
  left_join(
    hcdata$labels %>% select(x, label),
    by = "x"
  ) %>%
  left_join(working, by = "label") %>%
  mutate(label_y = y - 0.01)

hcdata$segments <- hcdata$segments %>%
  mutate(
    yend = if_else(yend == 0, y - 0.01, yend)
  ) %>%
  left_join(
    hcdata$labels %>% select(x, label),
    by = "x"
  ) %>%
  left_join(working, by = "label")


tree <- ggplot() +
  geom_hline( yintercept = 0.12, linetype = "dashed", color = "gray") +
  geom_segment(data = hcdata$segments,
               aes(x = x,y = y,xend = xend, yend = yend,color = group),
               linewidth = 0.5) +
  # Labels attached to individual terminal branches
  geom_text(data = tip_labels,
            aes(x = x, y = label_y, label = label, color = group),
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 1,
            show.legend = FALSE) +
  
  scale_color_discrete(
    limits = c(
      "Technical Replicates",
      "Inferred Clones",
      "Sample"
    ),
    name = "group",
    na.value = "gray50"
  ) +
  
  scale_y_continuous(
    limits = c(0, 0.4),
    breaks = seq(0.05, 0.4, 0.05),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  ylab("Genetic Distance (1-IBS)") +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.position = c(0.86, 0.86)
  )

tree

ggsave(filename = "/Users/Mustofa/Desktop/ibsdist211bams/IBS_bam211pt12_label.png",
       plot = tree,
       width = 9,
       height = 5,
       dpi = 300)
