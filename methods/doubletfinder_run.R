# DoubletFinder, following the package README workflow:
#   Seurat preprocessing -> pK chosen by max BCmetric (paramSweep, no ground truth)
#   -> doubletFinder(PCs = 1:10, pN = 0.25, nExp = round(rate * n)).
# No homotypic adjustment of nExp (needs cluster annotations); score = pANN.
local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        source(file.path(dirname(normalizePath(sub("^--file=", "", a))), "common.R")) })
suppressPackageStartupMessages({ library(Seurat); library(DoubletFinder) })

args <- parse_cli(); d <- load_counts(args)
enter_tmpdir()  # find.pK draws a plot -> Rplots.pdf
set.seed(args$seed)
PCs <- 1:10

seu <- CreateSeuratObject(counts = d$m)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE)

sweep <- paramSweep(seu, PCs = PCs, sct = FALSE)
bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
nExp <- round(d$rate * ncol(seu))
message(sprintf("DoubletFinder: pK = %g, nExp = %d", pK, nExp))

seu <- doubletFinder(seu, PCs = PCs, pN = 0.25, pK = pK, nExp = nExp, sct = FALSE)
md <- seu@meta.data
score <- md[[grep("^pANN_", colnames(md), value = TRUE)]]
call  <- as.integer(md[[grep("^DF.classifications_", colnames(md), value = TRUE)]] == "Doublet")
write_out(args, colnames(seu), score, call, d$rate)
