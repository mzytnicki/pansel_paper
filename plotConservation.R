library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
inputFileName <- args[1]
sizeThreshold <- as.numeric(args[2])
outputFileName <- args[3]

p <- inputFileName |>
    readr::read_tsv(col_names = c("chr", "start", "end", "id", "conservation", "strand"), show_col_types = FALSE) |>
    dplyr::mutate(size = end - start) |>
    dplyr::mutate(large_region = (size > sizeThreshold)) |>
    dplyr::select(conservation, large_region) |>
    ggplot2::ggplot(ggplot2::aes(x = large_region, y = conservation)) +
    ggplot2::geom_jitter(alpha = 0.1) +
    ggplot2::geom_violin(draw_quantiles = c(0.5)) +
    ggplot2::scale_y_log10() +
    ggplot2::ylim(0.9, 1)
ggplot2::ggsave(outputFileName, p)
