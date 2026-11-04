# scds: cxds (co-expression, no artificial doublets), bcds (boosting on artificial doublets),
# hybrid (cxds_bcds_hybrid). Defaults; estNdbl = TRUE gives the package's own doublet calls.
# None of them take an expected doublet rate (excluded from the rate-robustness experiment).
#   tool = cxds | bcds | hybrid
local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        source(file.path(dirname(normalizePath(sub("^--file=", "", a))), "common.R")) })
suppressPackageStartupMessages({ library(SingleCellExperiment); library(scds) })
patch_scds_xgboost()  # see common.R

args <- parse_cli(); d <- load_counts(args)
sce <- SingleCellExperiment(list(counts = d$m))
set.seed(args$seed)

sce <- switch(args$tool,
  cxds   = cxds(sce, estNdbl = TRUE),
  bcds   = bcds(sce, estNdbl = TRUE),
  hybrid = cxds_bcds_hybrid(sce, estNdbl = TRUE),
  stop("unknown tool ", args$tool))

score <- colData(sce)[[paste0(args$tool, "_score")]]
callc <- colData(sce)[[paste0(args$tool, "_call")]]
write_out(args, colnames(sce), score, if (is.null(callc)) NULL else as.integer(callc), d$rate)
