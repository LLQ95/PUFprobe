#!/usr/bin/env Rscript
#
# plot_network.R — cross-species hypothetical-protein family network figures.
#
# Produces three publication-ready figures from the tables written by
# cross_species_network.py:
#   1. community_network.*  macro-cluster graph (node size = protein count),
#   2. family_network.*     core family graph coloured by community,
#   3. community_bars.*     ranked bar chart showing which cluster dominates.
#
# Figures are written as both PDF and SVG (the SVG is editable in Adobe
# Illustrator). Designed to run on Windows R >= 4.6 and Linux R >= 4.6.
#
# Usage:
#   Rscript plot_network.R <network_dir> [output_dir]
#
# Defaults: network_dir = ../results/cross_species/network relative to this
# script's run; output_dir = network_dir.

options(stringsAsFactors = FALSE)

# ---- required packages -------------------------------------------------------
.cran_pkgs <- c("ggplot2", "dplyr", "readr", "ggrepel", "patchwork")
.missing <- .cran_pkgs[!vapply(.cran_pkgs, requireNamespace, logical(1),
                               quietly = TRUE)]
if (length(.missing) > 0) {
  message("Installing missing packages: ", paste(.missing, collapse = ", "))
  install.packages(.missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
network_dir <- if (length(args) >= 1) args[[1]] else
  file.path("..", "results", "cross_species", "network")
out_dir <- if (length(args) >= 2) args[[2]] else network_dir
network_dir <- normalizePath(network_dir, mustWork = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load --------------------------------------------------------------------
nodes <- read_tsv(file.path(network_dir, "nodes.tsv"), show_col_types = FALSE,
                  col_types = cols())
communities <- read_tsv(file.path(network_dir, "communities.tsv"),
                        show_col_types = FALSE)
comm_edges <- read_tsv(file.path(network_dir, "community_edges.tsv"),
                       show_col_types = FALSE)
core_edges <- read_tsv(file.path(network_dir, "core_edges.tsv"),
                       show_col_types = FALSE)

# numeric coords may be empty for non-core rows
nodes$x <- suppressWarnings(as.numeric(nodes$x))
nodes$y <- suppressWarnings(as.numeric(nodes$y))
communities$x <- suppressWarnings(as.numeric(communities$x))
communities$y <- suppressWarnings(as.numeric(communities$y))

# ---- palette -----------------------------------------------------------------
# A crisp, colour-blind-safe categorical palette; the anchor colours follow the
# project's warm/cyan scheme but the full set keeps clusters distinguishable.
cluster_palette <- c(
  "#E54E21", "#2F6DB5", "#8DD2C6", "#F2A20C", "#7B5EA7",
  "#3A9D54", "#D95F8A", "#5B9BD5", "#C77F2B", "#4575B4",
  "#FC8D59", "#91BF77", "#B07AA1", "#8C8C8C", "#E15759",
  "#4E9F9E")
community_ids <- sort(unique(communities$community))
comm_levels <- paste0("C", community_ids)
community_colors <- setNames(
  cluster_palette[((seq_along(community_ids) - 1L) %% length(cluster_palette)) + 1L],
  comm_levels)
nodes$community <- factor(paste0("C", nodes$community), levels = comm_levels)
communities$community <- factor(paste0("C", communities$community),
                                levels = comm_levels)
communities$cid <- factor(communities$community, levels = comm_levels)

# ---- shared theme ------------------------------------------------------------
clean_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey30"))
theme_set(clean_theme)

# ---- 1. community-level network ---------------------------------------------
ce <- comm_edges %>%
  mutate(community_a = paste0("C", community_a),
         community_b = paste0("C", community_b)) %>%
  left_join(communities %>% select(community, xa = x, ya = y),
            by = c("community_a" = "community")) %>%
  left_join(communities %>% select(community, xb = x, yb = y),
            by = c("community_b" = "community"))

communities <- communities %>%
  mutate(short_label = sub(";.*", "", label),
         node_label = paste0(community, "\n", n_proteins, " proteins"))

p_comm <- ggplot(communities, aes(x = x, y = y)) +
  geom_segment(data = ce, aes(x = xa, y = ya, xend = xb, yend = yb,
                              linewidth = shared_edges),
               colour = "grey70", alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(size = n_proteins, fill = community), shape = 21,
             colour = "white", stroke = 1.1) +
  geom_text(aes(label = community), size = 3.4, fontface = "bold",
            colour = "white") +
  geom_text_repel(aes(label = short_label), size = 2.7, colour = "grey20",
                  max.overlaps = 30, seed = 42, force = 6,
                  box.padding = 1.0, min.segment.length = 0) +
  scale_fill_manual(values = community_colors, name = "Community") +
  scale_size_continuous(range = c(12, 30), name = "Protein count",
                        guide = "none") +
  scale_linewidth_continuous(range = c(0.3, 3), guide = "none") +
  labs(title = "Cross-species hypothetical-protein communities",
       subtitle = "Macro-clusters detected by Leiden on family embeddings") +
  coord_fixed()

# ---- 2. family-level core network -------------------------------------------
core <- nodes %>% filter(core == 1, !is.na(x), !is.na(y))
fe <- core_edges %>%
  left_join(core %>% select(family_id, xa = x, ya = y),
            by = c("family_a" = "family_id")) %>%
  left_join(core %>% select(family_id, xb = x, yb = y),
            by = c("family_b" = "family_id")) %>%
  filter(!is.na(xa), !is.na(xb))

# label only the largest family in each community to avoid clutter
label_families <- core %>%
  group_by(community) %>%
  arrange(desc(size), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

p_fam <- ggplot(core, aes(x = x, y = y)) +
  geom_segment(data = fe, aes(x = xa, y = ya, xend = xb, yend = yb),
               colour = "grey85", alpha = 0.35, show.legend = FALSE) +
  geom_point(aes(size = size, colour = community), alpha = 0.75) +
  geom_text_repel(data = label_families,
                  aes(label = representative), size = 2.4,
                  colour = "grey15", max.overlaps = 40, seed = 42,
                  force = 8, box.padding = 0.6, min.segment.length = 0) +
  scale_colour_manual(values = community_colors, name = "Community") +
  scale_size_continuous(range = c(0.8, 4.5), name = "Proteins in family") +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1),
                               nrow = 4),
         size = guide_legend(nrow = 1)) +
  labs(title = "Core hypothetical-protein families",
       subtitle = "Multi-member families coloured by community; labels show the largest per community") +
  coord_fixed()

# ---- 3. ranked community bar chart ------------------------------------------
# Order by protein count ascending so the largest community sits at the top.
comm_bar <- communities %>%
  arrange(n_proteins) %>%
  mutate(ylabel = paste0(as.character(cid), ": ", sub(";.*", "", label))) %>%
  mutate(ylabel = factor(ylabel, levels = ylabel))

p_bar <- ggplot(comm_bar, aes(x = n_proteins, y = ylabel,
                              fill = cid)) +
  geom_col(width = 0.72, colour = "white") +
  geom_text(aes(label = n_proteins), hjust = -0.15, size = 3.1) +
  scale_fill_manual(values = community_colors, guide = "none", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Community sizes",
       subtitle = "Number of hypothetical proteins assigned to each cluster",
       x = "Number of proteins", y = NULL) +
  theme(axis.text.y = element_text(size = 9),
        plot.margin = margin(10, 30, 10, 10))

# ---- write -------------------------------------------------------------------
save_fig <- function(plot, name, width, height) {
  pdf(file.path(out_dir, paste0(name, ".pdf")),
      width = width, height = height)
  print(plot); dev.off()
  png(file.path(out_dir, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  print(plot); dev.off()
  # grDevices::svg is built in (no svglite system dependency) and produces a
  # vector SVG editable in Adobe Illustrator.
  svg(file.path(out_dir, paste0(name, ".svg")),
      width = width, height = height)
  print(plot); dev.off()
}
save_fig(p_comm, "community_network", 9, 7.5)
save_fig(p_fam, "family_network", 9, 8)
save_fig(p_bar, "community_bars", 8.5, 0.5 * max(3, nrow(comm_bar)) + 1.5)

combined <- p_comm + p_bar + plot_layout(ncol = 1, heights = c(1.5, 1))
pdf(file.path(out_dir, "network_overview.pdf"), width = 10, height = 11)
print(combined); dev.off()
png(file.path(out_dir, "network_overview.png"), width = 10, height = 11,
    units = "in", res = 300)
print(combined); dev.off()

message("Network figures written to ", out_dir)
