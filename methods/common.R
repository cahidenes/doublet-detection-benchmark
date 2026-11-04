# Shared I/O for R tool wrappers. Same output schema as methods/common.py:
#   barcode  score  call  tool  dataset  seed  expected_rate
# CLI (positional): <input_dir> <output.tsv> <dataset> <tool> <seed> [expected_rate|NA]
suppressPackageStartupMessages(library(Matrix))

parse_cli <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(a) >= 5)
  list(input = normalizePath(a[1]),
       output = file.path(normalizePath(dirname(a[2]), mustWork = FALSE), basename(a[2])),
       dataset = a[3], tool = a[4], seed = as.integer(a[5]),
       rate = if (length(a) >= 6 && a[6] != "NA") as.numeric(a[6]) else NA_real_)
}

# genes x droplets dgCMatrix with names; fills args$rate with the labeled rate if NA
load_counts <- function(args) {
  m <- as(readMM(file.path(args$input, "matrix.mtx.gz")), "CsparseMatrix")
  rownames(m) <- make.unique(readLines(file.path(args$input, "genes.tsv")))
  colnames(m) <- readLines(file.path(args$input, "barcodes.tsv"))
  lab <- read.delim(file.path(args$input, "labels.tsv"))
  stopifnot(identical(lab$barcode, colnames(m)))
  if (is.na(args$rate)) args$rate <- mean(lab$is_doublet)
  list(m = m, rate = args$rate)
}

# scds 1.22 bcds() reads two xgb.cv() result fields that xgboost >= 2.0 moved:
#   $best_iteration -> $early_stop$best_iteration   (NULL -> nrounds = Inf -> crash)
#   $pred           -> $cv_predict$pred (n x 1 matrix) (NULL -> empty bcds_score)
# (xgboost 1.x is not built for R 4.4; 3.x also breaks bcds' xgboost() call.) This wraps xgb.cv
# only inside scds' imports so the fields are where scds expects them; scds' own code and all
# numerics are unchanged.
patch_scds_xgboost <- function() {
  imp <- parent.env(asNamespace("scds"))
  orig <- get("xgb.cv", envir = imp)
  wrapped <- function(...) {
    r <- orig(...)
    if (is.null(r$best_iteration) && !is.null(r$early_stop$best_iteration))
      r$best_iteration <- r$early_stop$best_iteration
    if (is.null(r$pred) && !is.null(r$cv_predict$pred))
      r$pred <- as.vector(r$cv_predict$pred)
    r
  }
  unlockBinding("xgb.cv", imp); assign("xgb.cv", wrapped, envir = imp); lockBinding("xgb.cv", imp)
}

# Some tools write plots/intermediate files into the working directory: isolate them.
enter_tmpdir <- function() {
  d <- tempfile("dbl_"); dir.create(d); setwd(d); invisible(d)
}

write_out <- function(args, barcodes, score, call, rate) {
  if (is.null(call)) call <- rep(NA_integer_, length(barcodes))
  stopifnot(length(score) == length(barcodes), length(call) == length(barcodes))
  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  write.table(data.frame(barcode = barcodes, score = as.numeric(score), call = as.integer(call),
                         tool = args$tool, dataset = args$dataset, seed = args$seed,
                         expected_rate = rate),
              args$output, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
