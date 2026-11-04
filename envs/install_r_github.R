# Run inside the dbl-r env:  Rscript envs/install_r_github.R
# Log any failure in docs/usability_log.md (it is a usability result).
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Pinned to exact commits so the benchmark is reproducible.
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d", upgrade = "never")  # 2025-03-21
remotes::install_github("13308204545/Chord@00ac9a101f", upgrade = "never")                  # 2022-06-03

for (p in c("Seurat", "DoubletFinder", "scDblFinder", "scds", "Chord")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("%-14s %s  %s\n", p, if (ok) "OK  " else "FAIL",
              if (ok) as.character(packageVersion(p)) else ""))
}
