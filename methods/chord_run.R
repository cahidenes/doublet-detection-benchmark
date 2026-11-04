# Chord 2.0.1 (ensemble of cxds, bcds, DoubletFinder via gbm), default parameters.
# chord() needs the expected doublet rate; its call = top round(rate * n) by Chord score.
#
# Compatibility shims (Chord was written for DoubletFinder < 2.0.4 and Seurat v4):
#  1. it calls paramSweep_v3 / doubletFinder_v3, renamed to paramSweep / doubletFinder;
#  2. it reads seu@assays$RNA@counts, which Seurat v5 "Assay5" objects lack -> force v3 assays;
#  3. it passes reuse.pANN = FALSE, which DoubletFinder 2.0.6 treats as "reuse" (!is.null(FALSE)).
# The shims are defined in the global env, which R searches after Chord's namespace/imports.
# Chord's code itself is unmodified.
local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        source(file.path(dirname(normalizePath(sub("^--file=", "", a))), "common.R")) })
options(Seurat.object.assay.version = "v3")
suppressPackageStartupMessages({ library(Seurat); library(Chord) })
patch_scds_xgboost()  # 4. Chord runs scds::bcds internally (see common.R)

paramSweep_v3 <- function(...) DoubletFinder::paramSweep(...)
doubletFinder_v3 <- function(seu, ..., reuse.pANN = FALSE)
  DoubletFinder::doubletFinder(seu, ..., reuse.pANN = if (isFALSE(reuse.pANN)) NULL else reuse.pANN)

args <- parse_cli(); d <- load_counts(args)
enter_tmpdir()  # chord() writes CSVs, PDFs and .robj files into the working directory

seu <- CreateSeuratObject(counts = d$m)
called <- chord(seu = seu, doubletrate = d$rate, seed = args$seed)  # default method = "gbm"

sc <- read.csv("outreal_score.csv", row.names = 1)   # columns bcds_s, cxds_s, dbf_s, chord
stopifnot(nrow(sc) == ncol(seu))
barcodes <- colnames(seu)
# chord() returns rownames of its score table; these may be row indices rather than barcodes
idx <- if (all(called %in% barcodes)) match(called, barcodes) else as.integer(called)
stopifnot(!anyNA(idx))
call <- integer(length(barcodes)); call[idx] <- 1L
write_out(args, barcodes, sc$chord, call, d$rate)
