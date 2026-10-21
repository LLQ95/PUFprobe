#!/usr/bin/env Rscript
# plot_benchmark.R — render Swiss-Prot hold-out benchmark figures.
#
# Usage: Rscript plot_benchmark.R <benchmark_dir> [out_dir]
#
# Reads benchmark_summary.tsv, metrics_by_identity.tsv and metrics_by_phylum.tsv
# produced by baseline_holdout.py and writes PDF, PNG and SVG versions of:
#   benchmark_channels   coverage and any-correct rate per channel
#   benchmark_identity   any-correct rate across sequence-identity bins
#   benchmark_phylum     per-phylum recovery by channel
#
# Base graphics devices are used so the script runs without svglite/cairo.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript plot_benchmark.R <benchmark_dir> [out_dir]")
bench_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- if (length(args) >= 2) args[2] else bench_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

channel_colors <- c(seq = "#2F6DB5", struct = "#E54E21",
                    union = "#3A9D54", consensus = "#7B5EA7")
channel_labels <- c(seq = "Sequence", struct = "Structure",
                    union = "Fusion (union)", consensus = "Consensus")
channel_levels <- c("seq", "struct", "union", "consensus")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5),
    axis.title = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_blank()
  )

render_all <- function(plot, base) {
  for (ext in c("pdf", "png", "svg")) {
    dev_path <- file.path(out_dir, paste0(base, ".", ext))
    if (ext == "pdf") pdf(dev_path, width = 8, height = 5.2)
    if (ext == "png") png(dev_path, width = 2400, height = 1560, res = 300)
    if (ext == "svg") svg(dev_path, width = 8, height = 5.2)
    print(plot)
    invisible(dev.off())
  }
}

# ---- 1. channel summary ----------------------------------------------------
summary_path <- file.path(bench_dir, "benchmark_summary.tsv")
if (file.exists(summary_path)) {
  summary <- read_tsv(summary_path, show_col_types = FALSE) %>%
    mutate(channel = factor(channel, levels = channel_levels, ordered = TRUE))

  long <- summary %>%
    select(channel, coverage, any_correct = any_correct_rate) %>%
    pivot_longer(-channel, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
                           levels = c("coverage", "any_correct"),
                           labels = c("Coverage (has a prediction)",
                                      "Any-correct rate")))

  p_channels <- ggplot(long, aes(x = channel, y = value, fill = channel)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.35, size = 3.4) +
    facet_wrap(~ metric) +
    scale_x_discrete(labels = channel_labels) +
    scale_fill_manual(values = channel_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    labs(title = "Hold-out performance by evidence channel",
         subtitle = "Swiss-Prot curated Pfam labels hidden during prediction",
         x = NULL, y = "Fraction of proteins") +
    base_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  render_all(p_channels, "benchmark_channels")
}

# ---- 2. identity breakdown (the key structural-gain figure) ----------------
identity_path <- file.path(bench_dir, "metrics_by_identity.tsv")
if (file.exists(identity_path)) {
  identity <- read_tsv(identity_path, show_col_types = FALSE)
  id_levels <- c("<20%", "20-30%", "30-50%", ">=50%", "no hit")
  identity <- identity %>%
    filter(identity_bin %in% id_levels) %>%
    mutate(identity_bin = factor(identity_bin,
                                 levels = id_levels, ordered = TRUE)) %>%
    pivot_longer(c(seq_any_correct, struct_any_correct, union_any_correct),
                 names_to = "channel", values_to = "any_correct") %>%
    mutate(channel = recode(channel,
                            seq_any_correct = "seq",
                            struct_any_correct = "struct",
                            union_any_correct = "union"),
           channel = factor(channel,
                            levels = c("seq", "struct", "union"),
                            ordered = TRUE))

  p_identity <- ggplot(identity,
                       aes(x = identity_bin, y = any_correct,
                           fill = channel)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.72) +
    geom_text(aes(label = sprintf("%.2f", any_correct), group = channel),
              position = position_dodge(width = 0.82),
              vjust = -0.35, size = 2.8, show.legend = FALSE) +
    scale_fill_manual(
      name = "Channel",
      values = channel_colors,
      labels = c(seq = "Sequence", struct = "Structure",
                 union = "Fusion (union)"),
      drop = FALSE) +
    scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    labs(title = "Structural recovery in the remote-homology regime",
         subtitle = "Any-correct rate by identity to the nearest non-self hit",
         x = "Sequence identity to nearest non-self homologue",
         y = "Any-correct Pfam recovery") +
    base_theme
  render_all(p_identity, "benchmark_identity")
}

# ---- 3. phylum breakdown ---------------------------------------------------
phylum_path <- file.path(bench_dir, "metrics_by_phylum.tsv")
if (file.exists(phylum_path)) {
  phylum <- read_tsv(phylum_path, show_col_types = FALSE) %>%
    arrange(desc(n)) %>%
    head(12) %>%
    mutate(phylum = factor(phylum, levels = rev(phylum), ordered = TRUE)) %>%
    pivot_longer(c(seq_correct, struct_correct, union_correct),
                 names_to = "channel", values_to = "any_correct") %>%
    mutate(channel = recode(channel,
                            seq_correct = "seq",
                            struct_correct = "struct",
                            union_correct = "union"),
           channel = factor(channel,
                            levels = c("seq", "struct", "union"),
                            ordered = TRUE))

  p_phylum <- ggplot(phylum,
                     aes(x = phylum, y = any_correct, fill = channel)) +
    geom_col(position = position_dodge(width = 0.82), width = 0.72) +
    scale_fill_manual(
      name = "Channel",
      values = channel_colors,
      labels = c(seq = "Sequence", struct = "Structure",
                 union = "Fusion (union)"),
      drop = FALSE) +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    coord_flip() +
    labs(title = "Hold-out recovery across phyla",
         subtitle = "Top 12 phyla by sample size",
         x = NULL, y = "Any-correct Pfam recovery") +
    base_theme +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey90"))
  render_all(p_phylum, "benchmark_phylum")
}

message("Benchmark figures written to ", out_dir)
