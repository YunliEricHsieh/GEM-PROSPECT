# Run the existing R postprocessing from this checkout, with explicit stages.
# Usage: Rscript scripts/reproduce_chlamydomonas.R associations|figures|all
# No model parameters, candidate rules, or statistical calculations are changed.
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[[1]] else "associations"
if (length(args) > 1L || !stage %in% c("associations", "figures", "all")) {
  stop("Usage: Rscript scripts/reproduce_chlamydomonas.R associations|figures|all")
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this wrapper with Rscript.")
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
setwd(repo_root)

stage_scripts <- list(
  associations = "identify gene reaction associations.R",
  figures = c("Fig2.R", "Fig3.R")
)
selected <- if (stage == "all") unlist(stage_scripts, use.names = FALSE) else stage_scripts[[stage]]
packages <- c("dplyr", "tidyr")
if (stage %in% c("associations", "all")) packages <- c(packages, "purrr")
if (stage %in% c("figures", "all")) {
  packages <- c(packages, "ggplot2", "patchwork", "scales", "ggbeeswarm", "gridExtra", "svglite")
}
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required R packages first: ", paste(missing_packages, collapse = ", "))
}
dir.create("Results/figures", recursive = TRUE, showWarnings = FALSE)
for (script_name in selected) {
  message("Running existing script: ", script_name)
  source(file.path("Code", "R", script_name), local = new.env(parent = globalenv()))
}
