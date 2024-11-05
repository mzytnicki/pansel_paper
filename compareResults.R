library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]
outputFileName <- args[2]

readFile <- function(toolName, res, chr) {
  file.path(directory, toolName, res, paste0("chr", chr, ".tsv")) |>
    readr::read_tsv(col_names = c("id", "start", "end", "index", "n", "m", "s1", "e1", "z1", "s2", "s3", "z2", "e2", "e3"), show_col_types = FALSE) |>
    dplyr::select("start", "end", "index") |>
    dplyr::mutate(tool = toolName, resolution = res, chromosome = chr)
}

tools <- list.dirs(path = directory, full.names = FALSE, recursive = FALSE)
resolutions <- list.dirs(path = file.path(directory, tools[[1]]), full.names = FALSE, recursive = FALSE)
chrs <- seq(22)

table <- tidyr::expand_grid(toolName = tools, res = resolutions, chr = chrs) |>
  purrr::pmap(readFile) |>
  purrr::list_rbind() |>
  tidyr::pivot_wider(names_from = "tool", values_from = "index") |>
  tidyr::drop_na()
p <- table |>
  ggplot2::ggplot(ggplot2::aes(x = MGC, y = PGGB)) +
    ggplot2::geom_point(alpha = 0.1) +
    ggplot2::facet_grid(rows = vars(resolution))
ggplot2::ggsave(outputFileName, p, height = 8, width = 4)

for (res in c("1000", "10000", "100000")) {
  t <- table |> dplyr::filter(resolution == res)
  pearson <- cor(t$MGC, t$PGGB, method = "pearson")
  message(paste0("Pearson for ", res, ": ", pearson))
}
