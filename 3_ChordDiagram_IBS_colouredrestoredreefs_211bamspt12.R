setwd("/Users/mustofa/Desktop/ibsdist211bams")

library(tidyverse)
library(circlize)

# Input CSV columns:
# MLG, number_of_samples, number_of_populations, populations, samples
shared <- read_csv("MLG_IBSsummary_cutoff_0.12_211bams_noTR.csv", show_col_types = FALSE) %>%
  mutate(
    MLG = as.character(MLG),
    samples = as.character(samples)
  )

# Assign each site to one of the three groups.
site_groups <- c(
  OR1 = "Group 1",
  OR2 = "Group 1",
  OR3 = "Group 1",
  BR1 = "Group 2",
  SBR = "Group 2",
  SKR = "Group 2",
  GBH = "Group 3",
  GCH2 = "Group 3",
  SKH = "Group 3"
)

site_order <- names(site_groups)

group_colors <- c(
  "Group 1" = "#D55E00",
  "Group 2" = "#0072B2",
  "Group 3" = "#009E73"
)

site_colors <- setNames(
  unname(group_colors[site_groups]),
  names(site_groups)
)

# Small gaps within groups and larger gaps between groups.
gap_values <- c(
  OR1 = 2, OR2 = 2, OR3 = 12,
  BR1 = 2, SBR = 2, SKR = 12,
  GBH = 2, GCH2 = 2, SKH = 12
)

# Expand the sample lists and identify the site from each sample name.
# Assumes sample names start with a site code followed by a hyphen,
# such as BR1-63.
mlg_samples <- shared %>%
  select(MLG, samples) %>%
  separate_rows(samples, sep = ",\\s*") %>%
  mutate(
    samples = str_trim(samples),
    site = str_extract(samples, "^[^-]+")
  ) %>%
  filter(site %in% site_order) %>%
  distinct(MLG, site, samples)

if (nrow(mlg_samples) == 0) {
  stop(
    "No samples matched the site names. Check the 'samples' column ",
    "and the site prefix in your sample names."
  )
}

# Count ramets (samples) for each MLL at each site.
mlg_site_counts <- mlg_samples %>%
  group_by(MLG, site) %>%
  summarise(
    number_of_samples = n_distinct(samples),
    sample_names = paste(sort(unique(samples)), collapse = ", "),
    .groups = "drop"
  )

# MLLs found in more than one ramet at the same site.
# These will be shown as curved within-site loops.
within_site_table <- mlg_site_counts %>%
  filter(number_of_samples > 1) %>%
  mutate(
    within_site_ramet_pairs = choose(number_of_samples, 2)
  ) %>%
  arrange(site, MLG)

# One row per MLL and site for between-site comparisons.
mlg_sites <- mlg_site_counts %>%
  distinct(MLG, site)

# Create every unique pair of sites for each MLL.
# An MLL found at three sites contributes to all three site pairs.
pair_links <- mlg_sites %>%
  group_by(MLG) %>%
  group_modify(function(.x, .y) {
    sites_for_mll <- site_order[site_order %in% unique(.x$site)]
    
    if (length(sites_for_mll) < 2) {
      return(tibble(from = character(), to = character()))
    }
    
    site_pairs <- combn(sites_for_mll, 2)
    
    tibble(
      from = site_pairs[1, ],
      to = site_pairs[2, ]
    )
  }) %>%
  ungroup() %>%
  count(from, to, name = "shared_MLLs") %>%
  arrange(from, to) %>%
  mutate(
    from_group = unname(site_groups[from]),
    to_group = unname(site_groups[to]),
    
    # Links among the restored reefs BR1, SBR, and SKR are blue.
    is_restored_pair =
      from_group == "Group 2" & to_group == "Group 2",
    
    # Links between Group 1 and Group 2 remain orange.
    is_group1_group2 =
      (from_group == "Group 1" & to_group == "Group 2") |
      (from_group == "Group 2" & to_group == "Group 1"),
    
    # Restored-to-restored links take priority and are colored blue.
    link_color = case_when(
      is_restored_pair ~ adjustcolor("#0072B2", alpha.f = 0.80),
      is_group1_group2 ~ adjustcolor("#D55E00", alpha.f = 0.80),
      TRUE ~ adjustcolor("#D9D9D9", alpha.f = 0.80)
    )
  )

# Count unique MLLs found in at least one site from both Groups 1 and 2.
group1_mll <- mlg_sites %>%
  filter(unname(site_groups[site]) == "Group 1") %>%
  distinct(MLG) %>%
  pull(MLG)

group2_mll <- mlg_sites %>%
  filter(unname(site_groups[site]) == "Group 2") %>%
  distinct(MLG) %>%
  pull(MLG)

total_shared_group1_group2 <- length(intersect(group1_mll, group2_mll))

# Total between-site link count for each site.
site_totals <- bind_rows(
  pair_links %>% transmute(site = from, shared_MLLs),
  pair_links %>% transmute(site = to, shared_MLLs)
) %>%
  group_by(site) %>%
  summarise(
    between_site_shared_MLLs = sum(shared_MLLs),
    .groups = "drop"
  )

# Number of distinct MLLs producing within-site loops at each site.
within_site_totals <- within_site_table %>%
  count(site, name = "within_site_MLLs")

# Include every site in the labels, including zero-link sites.
site_labels <- tibble(site = site_order) %>%
  left_join(site_totals, by = "site") %>%
  left_join(within_site_totals, by = "site") %>%
  mutate(
    between_site_shared_MLLs =
      replace_na(between_site_shared_MLLs, 0L),
    within_site_MLLs =
      replace_na(within_site_MLLs, 0L),
    label = paste0(
      site,
      "\nN=", between_site_shared_MLLs + within_site_MLLs
    )
  ) %>%
  select(site, label) %>%
  deframe()

# Print and save the between-site pair count table.
pair_count_table <- pair_links %>%
  select(from_group, from, to_group, to, shared_MLLs) %>%
  arrange(desc(shared_MLLs), from, to)

print(pair_count_table)

write_csv(
  pair_count_table,
  "shared_MLL_counts_by_site_pair.csv"
)

# Save the MLLs and samples responsible for within-site loops.
write_csv(
  within_site_table,
  "within_site_shared_MLLs.csv"
)

# Between-site chords.
between_edges <- pair_links %>%
  transmute(
    from,
    to,
    value = as.numeric(shared_MLLs),
    link_color
  )

# Within-site curved loops.
# Each loop's width represents the number of distinct shared MLLs at that site.
within_edges <- within_site_table %>%
  count(site, name = "value") %>%
  transmute(
    from = site,
    to = site,
    value = as.numeric(value),
    link_color = adjustcolor(site_colors[site], alpha.f = 0.75)
  )

# Add invisible placeholder self-links for sites that have no real links.
# This keeps every site sector, including SBR, in the diagram.
sites_in_edges <- unique(c(
  between_edges$from,
  between_edges$to,
  within_edges$from,
  within_edges$to
))

sites_needing_placeholder <- setdiff(site_order, sites_in_edges)

placeholder_edges <- tibble(
  from = sites_needing_placeholder,
  to = sites_needing_placeholder,
  value = rep(0.000001, length(sites_needing_placeholder)),
  link_color = rep("#FFFFFF00", length(sites_needing_placeholder))
)

all_edges <- bind_rows(
  between_edges,
  within_edges,
  placeholder_edges
)

# Calculate enough sector width for all links.
# A self-link uses two ends within its sector.
sector_widths <- setNames(rep(0, length(site_order)), site_order)

if (nrow(all_edges) > 0) {
  for (i in seq_len(nrow(all_edges))) {
    from_site <- all_edges$from[i]
    to_site <- all_edges$to[i]
    edge_value <- all_edges$value[i]
    
    if (from_site == to_site) {
      sector_widths[from_site] <-
        sector_widths[from_site] + (2 * edge_value)
    } else {
      sector_widths[from_site] <-
        sector_widths[from_site] + edge_value
      sector_widths[to_site] <-
        sector_widths[to_site] + edge_value
    }
  }
}

# Keep zero-link site sectors visible.
sector_widths <- pmax(sector_widths, 1)
sector_widths <- sector_widths[site_order]

draw_chord <- function() {
  circos.clear()
  on.exit(circos.clear(), add = TRUE)
  
  par(mar = c(0, 0, 0, 0))
  
  circos.par(
    start.degree = 90,
    gap.after = gap_values,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01),
    canvas.xlim = c(-1.2, 1.2),
    canvas.ylim = c(-1.2, 1.2)
  )
  
  chordDiagram(
    x = all_edges %>% select(from, to, value),
    order = site_order,
    grid.col = site_colors,
    col = all_edges$link_color,
    transparency = 0,
    self.link = 2,
    reduce = -1,
    xmax = sector_widths,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.15),
    link.border = NA
  )
  
  # Add site labels and counts outside the colored site sectors.
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      site <- CELL_META$sector.index
      
      circos.text(
        x = CELL_META$xcenter,
        y = CELL_META$ycenter,
        labels = site_labels[site],
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        col = site_colors[site],
        cex = 0.75,
        font = 2
      )
    },
    bg.border = NA
  )
}

# Save the PNG in the working directory.
png_file <- file.path(
  getwd(),
  "IBS_chordDiagram_pt12v2.png"
)

png(
  filename = png_file,
  width = 1800,
  height = 1800,
  res = 300,
  bg = "white"
)

draw_chord()
dev.off()

# Draw again so the plot appears in RStudio's Plots tab.
draw_chord()

cat(
  "\nPNG saved to:\n",
  normalizePath(png_file),
  "\n"
)

cat(
  "Unique MLLs shared between Group 1 and Group 2:",
  total_shared_group1_group2,
  "\n"
)

