#!/usr/bin/env Rscript

# This script prepares a local R package library for the rain plot annotation
# workflow. It installs the required CRAN/Bioconductor packages into a specific
# user-provided directory instead of relying on the system-wide R library.
#
# This is useful on HPC systems because users often do not have permission to
# install packages globally. Keeping packages in a local workflow-specific
# library also makes the pipeline more reproducible and easier to rerun later.

args <- commandArgs(trailingOnly = TRUE)

# The script expects exactly one command-line argument:
# the directory where R packages should be installed.
#
# Example:
#   Rscript ensure_r_environment.R src/rainplot_workflow/R_library
if (length(args) != 1) {
  stop("Usage: ensure_r_environment.R <library_dir>", call. = FALSE)
}

# Convert the provided library path into a normalized path.
# winslash = "/" keeps paths consistent across systems, especially if the
# script is ever run on Windows.
#
# mustWork = FALSE allows the path to be normalized even if the directory does
# not exist yet, because we create it in the next step.
library_dir <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)

# Create the local R library directory if it does not already exist.
# recursive = TRUE allows nested directories to be created in one step.
# showWarnings = FALSE avoids unnecessary warnings if the directory already exists.
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)

# Add the local library directory to the front of R's library search paths.
#
# This is important because R checks .libPaths() in order. By placing the local
# library first, the script ensures that packages are installed into and loaded
# from the workflow-specific library before checking system/global libraries.
.libPaths(c(library_dir, .libPaths()))

# List of packages required for the genomic annotation plotting workflow.
#
# BiocManager is needed first because most of the other packages come from
# Bioconductor rather than regular CRAN.
#
# Gviz is used to build genomic feature tracks.
# rtracklayer helps import/export genomic annotation formats.
# GenomicRanges is used for working with genomic intervals.
# GenomeInfoDb helps manage chromosome/genome metadata.
required_packages <- c(
  "BiocManager",
  "Gviz",
  "rtracklayer",
  "GenomicRanges",
  "GenomeInfoDb"
)

# Check whether a package is installed specifically in the local workflow
# library directory.
#
# We check the local library instead of checking all of R's library paths because
# the goal is to make sure this workflow has its own required packages available.
# A package may exist in a system library, but that does not guarantee it will be
# available later in the same HPC/job environment.
is_installed_in_library <- function(pkg) {
  pkg_path <- suppressWarnings(
    find.package(pkg, lib.loc = library_dir, quiet = TRUE)
  )

  # find.package() returns a path if the package is found.
  # length(pkg_path) > 0 confirms something was returned.
  # nzchar(pkg_path) confirms the returned path is not an empty string.
  length(pkg_path) > 0 && nzchar(pkg_path)
}

# Install BiocManager first if it is not already available in the local library.
#
# BiocManager is required before installing Bioconductor packages like Gviz,
# rtracklayer, GenomicRanges, and GenomeInfoDb.
if (!is_installed_in_library("BiocManager")) {
  message("[INFO] Installing BiocManager into the local R library...")

  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org",
    lib = library_dir
  )
}

# Set R's package repositories using BiocManager.
#
# This configures both CRAN and Bioconductor repositories correctly, which is
# important because the required packages may come from different sources.
options(repos = BiocManager::repositories())

# Identify which required Bioconductor packages are missing from the local
# library.
#
# BiocManager is excluded here because it was handled separately above.
# vapply() checks each package and returns a logical vector showing whether the
# package is installed in the local library.
missing_bioc_packages <- required_packages[required_packages != "BiocManager"][
  !vapply(
    required_packages[required_packages != "BiocManager"],
    is_installed_in_library,
    logical(1)
  )
]

# Install only the missing packages.
#
# This avoids reinstalling packages every time the workflow runs, which saves
# time and reduces unnecessary package changes.
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

# Print the final library location so the user or workflow log clearly shows
# where the R environment was prepared.
message(sprintf("[INFO] Local R library ready: %s", library_dir))