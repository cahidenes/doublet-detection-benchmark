# Convert one Xi & Li .rds (list: [[1]] genes x droplets count matrix, [[2]] labels
# "doublet"/"singlet") into a language-neutral folder that R and Python both read:
#   matrix.mtx.gz (genes x droplets), genes.tsv, barcodes.tsv, labels.tsv
# Usage: Rscript scripts/convert_rds.R <in.rds> <out_dir>
# NOTE: verify the list structure on the first file (str(readRDS(...))) after download.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
d <- readRDS(args[1]); out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

m <- as(d[[1]], "CsparseMatrix")
lab <- as.character(d[[2]])
stopifnot(ncol(m) == length(lab), all(lab %in% c("doublet", "singlet")))
if (is.null(colnames(m))) colnames(m) <- sprintf("cell%06d", seq_len(ncol(m)))
if (is.null(rownames(m))) rownames(m) <- sprintf("gene%06d", seq_len(nrow(m)))

mtx <- file.path(out, "matrix.mtx")
writeMM(m, mtx); system2("gzip", c("-f", mtx))
writeLines(rownames(m), file.path(out, "genes.tsv"))
writeLines(colnames(m), file.path(out, "barcodes.tsv"))
write.table(data.frame(barcode = colnames(m), is_doublet = as.integer(lab == "doublet")),
            file.path(out, "labels.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d genes x %d droplets, %.1f%% doublets\n",
            basename(args[1]), nrow(m), ncol(m), 100 * mean(lab == "doublet")))
