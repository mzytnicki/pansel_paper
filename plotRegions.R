library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
conservedFileName <- args[[1]]
divergentFileName <- args[[2]]
outputFileName <- args[[3]]

types <- c("Conserved", "Divergent")
chrs <- paste0("chr", seq(22))

readFile <- function(fileName, typeId) {
  fileName |>
    read_tsv(col_names = c("chr", "start", "end", "id", "score", "strand", "binSize")) |>
    mutate(type = typeId)
}

d <- c(conservedFileName, divergentFileName) |>
  map2_dfr(types, readFile) |>
  mutate(chr = factor(chr, levels = chrs)) |>
  mutate(type = factor(type))

p <- ggplot(d, aes(x = start, color = type, fill = type)) +
  geom_freqpoly(binwidth = 1000000) +
  facet_grid(rows = vars(chr)) +
  scale_color_discrete() +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0)) +
  theme(legend.title = element_blank()) +
  xlab("position")

ggplot2::ggsave(outputFileName, p, height = 10)
