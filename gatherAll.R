library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
inputFileNames <- head(args, -3)
types <- args[length(args) - 2] |>
  stringr::str_split(pattern = ",") |>
  purrr::pluck(1)
nTypes <- length(types)
breaksFileName <- args[length(args)-1]
outputFileName <- args[length(args)]

xlabels <- readLines(breaksFileName)
xlabels <- setNames(xlabels, 0:7)

p <- purrr::map_dfr(inputFileNames, readr::read_tsv, col_types = "fdc") |>
  dplyr::filter(type %in% types) |>
  dplyr::mutate(type = factor(type, levels = types)) |>
  ggplot2::ggplot(ggplot2::aes(x = category, y = score)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(rows = vars(type), scales = "free") +
  ggplot2::xlab("") +
  ggplot2::ylab("") +
  ggplot2::scale_x_discrete(labels = xlabels) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

ggplot2::ggsave(outputFileName, p, height = nTypes + 1.5)
