# =============================================================================
# snp_cnv_pipeline.R
# -----------------------------------------------------------------------------
# SNP calling and Copy Number Variation (CNV) analysis pipeline for
# genomic integrity assessment of reprogrammed cells (hiPSCs vs parental).
#
# Workflow:
#   1. Load and preprocess Illumina SNP array data (.idat files)
#   2. Genotype calling via CRLMM (krlmm method)
#   3. BAF (B-Allele Frequency) and CNV computation
#   4. Per-chromosome anomaly detection (duplication / deletion)
#   5. Export per-sample BAF/CNV tables (.tsv)
#   6. Export per-family PDF reports with chromosome ideograms
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/snp-calling-ipsc-integrity
# Project : Genomic integrity of human iPSC lines
# Paper   : Gaignerie A. et al., Scientific Reports, 2018
#           DOI: 10.1038/s41598-018-32645-2
# Date    : 2017 (MSc M1 internship, CRTI UMR 1064, Nantes Universite)
# =============================================================================


# =============================================================================
# I. DEPENDENCIES
# =============================================================================

.install_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' not found. Install with install.packages('", pkg, "')")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Bioconductor packages
for (pkg in c("illuminaio", "crlmm", "oligoClasses", "SNPchip")) .install_bioc(pkg)

# CRAN packages
for (pkg in c("png", "stringr", "matrixStats", "ff", "tools", "grid")) .load_pkg(pkg)


# =============================================================================
# II. PARAMETERS — edit this section before running
# =============================================================================

# Root directory containing .idat files, manifest, and samplesheet
DATA_DIR  <- "path/to/your/data"        # <-- set your data directory

# Samplesheet filename (CSV, Illumina format)
SAMPLESHEET <- "SampleSheet.csv"        # <-- set your samplesheet filename

# Family / cell-type annotation file (TSV with Patient, CellType, Sample_ID)
FAMILY_FILE <- "FamilyAnnotation.txt"   # <-- set your family annotation file

# Manifest file (CSV with ";" separator, Illumina format)
MANIFEST_FILE <- "ManifestFile.csv"     # <-- set your manifest filename

# Reference genome build
GENOME <- "hg19"

# CRLMM genotyping confidence threshold (SNPs below this are excluded)
CONF_THRESHOLD <- 0.90

# BAF anomaly detection parameters
BAF_DUP_WINDOW  <- 0.08    # sliding window size (Mb) for duplication detection
BAF_DUP_MIN_PTS <- 5       # min points in window to flag duplication
BAF_DEL_RATIO   <- 0.17    # max BAF heterozygous ratio below which deletion is flagged
BAF_MIN_MARKERS <- 2000    # min markers required for anomaly detection

# Chromosome numeric encoding for sex chromosomes
CHR_MAP <- c(X = 23L, Y = 24L, XY = 25L, MT = 26L)

# Characters to sanitise in filenames (Windows-incompatible)
ILLEGAL_CHARS <- c("/", ":", "<", ">", "|", '"')


# =============================================================================
# III. UTILITY FUNCTIONS
# =============================================================================

#' Create a directory if it does not exist
#' @param path Character. Directory path.
create_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


#' Sanitise a string for use as a filename
#' Replaces Windows-illegal characters with underscores.
#' @param x Character.
#' @return Character.
sanitise_filename <- function(x) {
  for (ch in ILLEGAL_CHARS) x <- gsub(ch, "_", x, fixed = TRUE)
  x
}


#' Read a samplesheet in Illumina CSV format
#' Illumina samplesheets have 8 metadata rows before the data rows.
#' @param filepath Character. Path to samplesheet CSV.
#' @return data.frame. Samplesheet with SentrixBarcode_A and SentrixPosition_A.
read_samplesheet <- function(filepath) {
  header_row  <- read.csv(filepath, header = FALSE, as.is = TRUE)[9, ]
  samplesheet <- read.csv(filepath, header = FALSE,
                          col.names = header_row, as.is = TRUE)
  samplesheet <- samplesheet[-c(1:9), 1:3]
  return(samplesheet)
}


#' Read and preprocess an Illumina manifest CSV
#' Encodes sex/MT chromosomes as integers and adds positional annotation columns.
#' @param filepath Character. Path to manifest CSV (semicolon-separated).
#' @return data.frame. Annotated manifest.
read_manifest <- function(filepath) {
  manifest_cols <- c(
    "IlmnID", "Name", "IlmnStrand", "SNP",
    "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq",
    "GenomeBuild", "chromosome", "MapInfo", "Ploidy", "Species",
    "Source", "SourceVersion", "SourceStrand", "SourceSeq",
    "TopGenomicSeq", "BeadSetID", "Exp_Clusters", "RefStrand"
  )
  manifest <- read.csv2(filepath, header = FALSE, fill = TRUE,
                        stringsAsFactors = FALSE)
  manifest <- manifest[-c(1:8), ]
  colnames(manifest) <- manifest_cols

  # Encode sex / mitochondrial chromosomes as integers
  for (sym in names(CHR_MAP)) {
    manifest$chromosome[manifest$chromosome == sym] <- as.character(CHR_MAP[sym])
  }
  manifest$chromosome <- as.integer(manifest$chromosome)
  manifest$isSnp      <- TRUE
  manifest$position   <- as.integer(manifest$MapInfo)

  # Add feature annotation columns required by genotype.Illumina
  annot <- manifest[, c("Name", "chromosome", "position")]
  colnames(annot) <- c("featureNames", "chr", "position")
  manifest <- cbind(manifest, annot)
  return(manifest)
}


# =============================================================================
# IV. BAF / CNV PLOTTING FUNCTION
# =============================================================================

#' Plot BAF and CNV for one chromosome and detect anomalies
#'
#' Generates a PNG with two panels (BAF + CNV + idiogram) and returns the
#' chromosome number if a genomic anomaly is detected (duplication or deletion).
#'
#' @param cn_reset    crlmmSet object. Output of genotype.Illumina / crlmmCopynumber.
#' @param chr         Integer. Chromosome number (1-23).
#' @param sample_idx  Integer. Sample column index in cn_reset.
#' @param sample_id   Character. Sample identifier (for plot title).
#' @param out_dir     Character. Directory for PNG output.
#' @param label_bands Logical. Label cytoband names on idiogram. Default FALSE.
#'
#' @return Integer (chr) if anomaly detected, else NULL (invisible).
plot_baf_cnv <- function(cn_reset,
                          chr,
                          sample_idx,
                          sample_id,
                          out_dir,
                          label_bands = FALSE) {
  markers <- which(chromosome(cn_reset) == chr)

  # Extract BAF and CNV
  A   <- CA(cn_reset, i = markers, j = sample_idx)
  B   <- CB(cn_reset, i = markers, j = sample_idx)
  BAF <- B / (A + B)
  cnv <- totalCopynumber(cn_reset, i = markers, j = sample_idx)
  pos <- position(cn_reset)[markers]

  # Build output filename
  png_file <- file.path(out_dir, paste0(sanitise_filename(sample_id),
                                        "_Chr", chr, ".png"))
  png(png_file)
  on.exit(dev.off())

  # Chromosomes 13/14/15/21/22 are acrocentric — use fixed xlim
  is_acrocentric <- chr %in% c(13, 14, 15, 21, 22)
  xlim_arg       <- if (is_acrocentric) c(0, max(pos)) else NULL

  par(mfcol = c(2, 1), mex = 0.5, cex = 0.5, mar = c(4.1, 0, 6.6, 8.0))
  chr_label <- paste("/ chr.", chr)

  # BAF panel
  plot(pos, BAF,
       cex.main = 4, pch = 16, las = 1, yaxt = "n", xlab = "",
       cex = 1.4, xaxt = "n", col = "chartreuse4",
       ylim = c(0, 1), xlim = xlim_arg, ylab = "BAF",
       main = paste(sample_id, paste(chr_label, "BAF", sep = "   :  "), sep = " /"))
  axis(side = 4, at = seq(0, 1, by = 0.5), labels = TRUE, las = 1, cex.axis = 3)

  # CNV panel
  plot(pos, cnv,
       cex.main = 4, pch = 16, las = 1, yaxt = "n", xlab = "",
       cex = 1.4, xaxt = "n", col = "blue",
       ylim = c(-2, 4), xlim = xlim_arg, ylab = "CNV",
       main = paste(sample_id, paste(chr_label, "CNV", sep = "   :  "), sep = " /"))
  axis(side = 4, at = seq(0, 4, by = 2), labels = TRUE, las = 1, cex.axis = 3)

  # Chromosome idiogram
  plotIdiogram(chr, GENOME,
               unit = c("bp", "Mb"), new = FALSE,
               label.cytoband = label_bands,
               cytoband.ycoords = c(-2, -1), label.y = c(-3), verbose = FALSE)

  # ── Anomaly detection ─────────────────────────────────────────────────────
  res <- data.frame(BAF = BAF, CNV = cnv, x = pos / 1e6)
  res <- res[complete.cases(res), ]
  res <- res[order(res$x), ]

  anomaly <- NULL

  if (nrow(res) > BAF_MIN_MARKERS) {
    n_below_half <- sum(res$BAF <= 0.50)
    balance_ratio <- n_below_half / nrow(res)

    if (balance_ratio > 0.45 && balance_ratio < 0.55) {
      bar_threshold <- max(7L, as.integer(nrow(res) / 700))

      # Check duplication (BAF bands at ~0.25 and ~0.75)
      for (band in list(
        res[res$BAF > 0.19 & res$BAF < 0.405, ],
        res[res$BAF > 0.595 & res$BAF < 0.81, ]
      )) {
        if (nrow(band) > BAF_DUP_MIN_PTS) {
          count <- 0L
          prev  <- -Inf
          for (p in band$x) {
            if (abs(p - prev) < BAF_DUP_WINDOW) {
              count <- count + 1L
            } else {
              count <- 0L
            }
            prev <- p
            if (count > bar_threshold && is.null(anomaly)) {
              warning(sprintf(
                "Possible duplication — Sample: %s | Chr: %d | BAF banding detected",
                sample_id, chr))
              anomaly <- chr
              break
            }
          }
        }
        if (!is.null(anomaly)) break
      }

      # Check deletion (very few heterozygous SNPs)
      n_het <- nrow(res[res$BAF > 0.20 & res$BAF < 0.80, ])
      if ((n_het / nrow(res)) < BAF_DEL_RATIO) {
        warning(sprintf(
          "Possible deletion — Sample: %s | Chr: %d | Low heterozygosity",
          sample_id, chr))
        anomaly <- chr
      }
    }
  }

  return(anomaly)
}


# =============================================================================
# V. PER-SAMPLE TSV EXPORT
# =============================================================================

#' Export BAF and CNV values for all autosomes to a TSV file
#'
#' @param cn_reset   crlmmSet object.
#' @param sample_idx Integer. Sample column index.
#' @param sample_id  Character. Sample identifier (used in filename).
#' @param out_dir    Character. Output directory.
#' @return Invisible NULL.
export_baf_cnv_tsv <- function(cn_reset, sample_idx, sample_id, out_dir) {
  create_dir(out_dir)
  autosomes <- 1:22

  markers <- unlist(lapply(autosomes, function(c)
    which(chromosome(cn_reset) == c)))

  A   <- CA(cn_reset, i = markers, j = sample_idx)
  B   <- CB(cn_reset, i = markers, j = sample_idx)
  BAF <- B / (A + B)
  cnv <- totalCopynumber(cn_reset, i = markers, j = sample_idx)

  out_table          <- data.frame(BAF = BAF, CNV = cnv)
  out_file           <- file.path(out_dir,
                                  paste0(sanitise_filename(sample_id), ".tsv"))
  write.table(out_table, file = out_file,
              sep = "\t", col.names = NA, dec = ".", quote = FALSE)
  message("Exported: ", out_file)
  invisible(NULL)
}


# =============================================================================
# VI. MAIN PIPELINE
# =============================================================================

message("=== SNP / CNV Genomic Integrity Pipeline ===")
message("[1/5] Loading samplesheet and manifest...")

samplesheet <- read_samplesheet(file.path(DATA_DIR, SAMPLESHEET))
manifest    <- read_manifest(file.path(DATA_DIR, MANIFEST_FILE))

# Build array file paths
array_names <- file.path(DATA_DIR,
                          paste(samplesheet$SentrixBarcode_A,
                                samplesheet$SentrixPosition_A, sep = "_"))
array_info  <- list(barcode = "SentrixBarcode_A", position = "SentrixPosition_A")

# Create ff scratch directory
ff_dir <- file.path(DATA_DIR, paste0("ff_files_", tools::file_path_sans_ext(SAMPLESHEET)))
create_dir(ff_dir)
ldPath(ff_dir)


message("[2/5] Running CRLMM genotyping (krlmm method)...")
cn_reset <- genotype.Illumina(
  sampleSheet        = samplesheet,
  arrayNames         = array_names,
  ids                = NULL,
  path               = "",
  arrayInfoColNames  = array_info,
  highDensity        = FALSE,
  sep                = "_",
  fileExt            = list(green = "Grn.idat", red = "Red.idat"),
  XY                 = NULL,
  anno               = manifest,
  genome             = GENOME,
  call.method        = "krlmm",
  trueCalls          = NULL,
  cdfName            = "nopackage",
  copynumber         = TRUE,
  batch              = NULL,
  saveDate           = TRUE,
  stripNorm          = TRUE,
  useTarget          = TRUE,
  quantile.method    = "between",
  nopackage.norm     = "loess",
  mixtureSampleSize  = 100,
  fitMixture         = TRUE,
  eps                = 0.1,
  verbose            = TRUE,
  seed               = 10,
  probs              = rep(1 / 3, 3),
  DF                 = 6,
  SNRMin             = 5,
  recallMin          = 2,
  recallRegMin       = 100,
  gender             = NULL,
  returnParams       = TRUE,
  badSNP             = 0.7
)
open(cn_reset)


message("[3/5] Computing BAF and CNV (confidence threshold: ", CONF_THRESHOLD, ")...")
crlmmCopynumber(
  cn_reset,
  MIN.SAMPLES   = 10,
  SNRMin        = 5,
  MIN.OBS       = 1,
  DF.PRIOR      = 50,
  bias.adj      = FALSE,
  prior.prob    = rep(1 / 4, 4),
  seed          = 1,
  verbose       = TRUE,
  GT.CONF.THR   = CONF_THRESHOLD,
  MIN.NU        = 8,
  MIN.PHI       = 8,
  THR.NU.PHI    = TRUE,
  type          = c("SNP", "NP", "X.SNP", "X.NP"),
  fit.linearModel = TRUE
)


message("[4/5] Generating per-family PDF reports...")
group      <- read.table(file.path(DATA_DIR, FAMILY_FILE),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
images_dir <- file.path(DATA_DIR, "images")
create_dir(images_dir)
results_dir <- file.path(DATA_DIR, "results")
create_dir(results_dir)

families <- unique(group$Patient)

for (patient in families) {
  message("  Family: ", patient)
  family_data <- group[group$Patient == patient, ]
  parental    <- family_data[family_data$CellType == "Parental" &
                               !is.na(family_data$Sample_ID), ]
  hipsc       <- family_data[family_data$CellType == "hiPSC" &
                               !is.na(family_data$Sample_ID), ]
  sample_rows <- c(as.integer(rownames(parental)), as.integer(rownames(hipsc)))

  family_name  <- sanitise_filename(as.character(patient))
  family_dir   <- file.path(images_dir, family_name)
  create_dir(family_dir)
  pdf_path     <- file.path(family_dir, paste0(family_name, ".pdf"))

  pdf(pdf_path, width = 10, height = 10)

  for (sample_idx in sample_rows) {
    sample_id    <- samplesheet$Sample_ID[sample_idx]
    anomalies    <- c()
    png_files    <- c()

    # ── Per-chromosome overview plots ──────────────────────────────────────
    for (chr in 1:23) {
      anomaly  <- plot_baf_cnv(cn_reset, chr, sample_idx,
                                sample_id, family_dir)
      png_files <- c(png_files,
                     file.path(family_dir,
                               paste0(sanitise_filename(sample_id),
                                      "_Chr", chr, ".png")))
      if (!is.null(anomaly)) anomalies <- c(anomalies, anomaly)
    }

    # ── Export BAF/CNV TSV for this sample ─────────────────────────────────
    export_baf_cnv_tsv(cn_reset, sample_idx, sample_id, results_dir)

    # ── Overview page: 23 chromosomes as raster tiles ──────────────────────
    par(mfcol = c(4, 6), mex = 0.5, cex = 0.5, mar = c(4.1, 0, 3.1, 2.1))
    for (chr in 1:23) {
      pn <- tryCatch(readPNG(png_files[chr]), error = function(e) NULL)
      unlink(png_files[chr])   # clean up temp PNG
      plot.new()
      if (!is.null(pn)) rasterImage(pn, 0, 0, 1, 1)
    }

    # ── Zoom page: one page per anomalous chromosome ───────────────────────
    for (anom_chr in anomalies) {
      par(mfrow = c(2, 1), mar = c(5.1, 4.1, 4.1, 2.1))
      markers  <- which(chromosome(cn_reset) == anom_chr)
      A        <- CA(cn_reset, i = markers, j = sample_idx)
      B        <- CB(cn_reset, i = markers, j = sample_idx)
      BAF      <- B / (A + B)
      cnv      <- totalCopynumber(cn_reset, i = markers, j = sample_idx)
      pos      <- position(cn_reset)[markers]
      chr_lab  <- paste("/ chr.", anom_chr)
      is_acro  <- anom_chr %in% c(13, 14, 15, 21, 22)
      xlim_arg <- if (is_acro) c(0, max(pos)) else NULL

      # BAF zoom
      plot(pos, BAF,
           cex.main = 0.9, pch = 16, las = 1, yaxt = "n", xlab = "",
           cex = 1, xaxt = "n", col = "chartreuse4",
           ylim = c(0, 1), xlim = xlim_arg, ylab = "BAF",
           main = paste(sample_id, paste(chr_lab, "BAF", sep = "   :  "), sep = " /"))
      axis(1, at = pretty(pos), labels = pretty(pos) / 1e6)
      axis(side = 2, at = seq(0, 1, by = 0.25), labels = TRUE, las = 1)

      # CNV zoom
      plot(pos, cnv,
           cex.main = 0.9, pch = 16, las = 1, yaxt = "n", xlab = "",
           cex = 1, xaxt = "n", col = "blue",
           ylim = c(-2, 4), xlim = xlim_arg, ylab = "CNV",
           main = paste(sample_id, paste(chr_lab, "CNV", sep = "   :  "), sep = " /"))
      axis(1, at = pretty(pos), labels = pretty(pos) / 1e6)
      axis(side = 2, at = seq(0, 4, by = 2), labels = TRUE, las = 1)

      plotIdiogram(anom_chr, GENOME,
                   unit = c("bp", "Mb"), new = FALSE,
                   label.cytoband = FALSE,
                   cytoband.ycoords = c(-2, -1), label.y = c(-3), verbose = FALSE)
    }
  }

  dev.off()
  message("  PDF saved: ", pdf_path)
}

message("[5/5] Done. Results in: ", results_dir)
