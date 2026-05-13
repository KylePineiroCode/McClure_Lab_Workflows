#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: ensure_r_environment.R <library_dir>", call. = FALSE)
}

library_dir <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(library_dir, .libPaths()))

required_packages <- c(
  "BiocManager",
  "Gviz",
  "rtracklayer",
  "GenomicRanges",
  "GenomeInfoDb"
)

is_installed_in_library <- function(pkg) {
  pkg_path <- suppressWarnings(
    find.package(pkg, lib.loc = library_dir, quiet = TRUE)
  )
  length(pkg_path) > 0 && nzchar(pkg_path)
}

if (!is_installed_in_library("BiocManager")) {
  message("[INFO] Installing BiocManager into the local R library...")
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org",
    lib = library_dir
  )
}

options(repos = BiocManager::repositories())

missing_bioc_packages <- required_packages[required_packages != "BiocManager"][
  !vapply(
    required_packages[required_packages != "BiocManager"],
    is_installed_in_library,
    logical(1)
  )
]

if (length(missing_bioc_packages) > 0) {
  message(sprintf(
    "[INFO] Installing missing R packages into the local R library: %s",
    paste(missing_bioc_packages, collapse = ", ")
  ))
  BiocManager::install(
    missing_bioc_packages,
    ask = FALSE,
    update = FALSE,
    lib = library_dir,
    force = TRUE
  )
} else {
  message("[INFO] Required R packages already installed in the local R library.")
}

message(sprintf("[INFO] Local R library ready: %s", library_dir))