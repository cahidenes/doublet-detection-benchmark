# scDblFinder (default parameters except dbr) and its computeDoubletDensity variant
# (the updated scran::doubletCells).   tool = scdblfinder | dbldensity
local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        source(file.path(dirname(normalizePath(sub("^--file=", "", a))), "common.R")) })
suppressPackageStartupMessages({ library(SingleCellExperiment); library(scDblFinder) })

args <- parse_cli(); d <- load_counts(args)
sce <- SingleCellExperiment(list(counts = d$m))
set.seed(args$seed)

if (args$tool == "scdblfinder") {
  sce <- scDblFinder(sce, dbr = d$rate)
  score <- sce$scDblFinder.score
  call  <- as.integer(sce$scDblFinder.class == "doublet")
} else if (args$tool == "dbldensity") {
  sce <- scuttle::logNormCounts(sce)
  score <- computeDoubletDensity(sce)
  call  <- NULL  # score only, no default threshold
} else stop("unknown tool ", args$tool)

write_out(args, colnames(sce), score, call, d$rate)
