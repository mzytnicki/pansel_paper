library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
inputFileNames <- head(args, -2)
types <- args[length(args) - 1] |>
  stringr::str_split(pattern = ",") |>
  purrr::pluck(1)
nTypes <- length(types)
outputFileName <- args[length(args)]

p <- purrr::map_dfr(inputFileNames, readr::read_tsv, col_types = "fdc") |>
  dplyr::filter(type %in% types) |>
  dplyr::mutate(type = factor(type, levels = types)) |>
  ggplot2::ggplot(ggplot2::aes(x = category, y = score)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(rows = vars(type), scales = "free") +
  ggplot2::theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  ggplot2::ylab("")

ggplot2::ggsave(outputFileName, p, height = nTypes)
