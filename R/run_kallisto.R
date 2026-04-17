#' Run kallisto-based quantification and collect expression matrices
#'
#' `MyMuse()` is a lightweight constructor that can:
#' 1) run `kallisto index` and `kallisto quant`, and
#' 2) load abundance tables into TPM / count / effective-length matrices.
#'
#' @param cdsfile Path to transcript FASTA file used for kallisto indexing.
#' @param fastqfiles_dir Directory containing FASTQ(.gz) files.
#' @param replicates Number of replicates per condition (used for averaging).
#' @param paired_end Logical; `TRUE` for paired-end, `FALSE` for single-end.
#' @param fragment_length Fragment length (single-end only).
#' @param sd Fragment length standard deviation (single-end only).
#' @param bootstrap_samples Number of kallisto bootstrap samples.
#' @param kallisto_bin Path to kallisto binary.
#' @param index_file Path for temporary kallisto index.
#' @param run_quant Logical; if `FALSE`, skip quant and only parse existing outputs.
#' @param cleanup_index Logical; remove index file after run.
#'
#' @return A `MyMuse` object (named list) with quantification matrices.
#' @export
MyMuse <- function(
  cdsfile,
  fastqfiles_dir,
  replicates,
  paired_end = TRUE,
  fragment_length = NULL,
  sd = NULL,
  bootstrap_samples = 100,
  kallisto_bin = "kallisto",
  index_file = "index.idx",
  run_quant = TRUE,
  cleanup_index = TRUE
) {
  stopifnot(file.exists(cdsfile), dir.exists(fastqfiles_dir), replicates >= 1)

  value <- list(
    cdsfile = cdsfile,
    fastqfiles_dir = fastqfiles_dir,
    reps = replicates,
    paired_end = paired_end,
    fragment_length = fragment_length,
    sd = sd,
    bootstrap_samples = bootstrap_samples,
    kallisto_bin = kallisto_bin,
    index_file = index_file,
    names = NULL,
    tpms = NULL,
    eff_lengths = NULL,
    est_counts = NULL,
    rpkms = NULL,
    average_rpkms = NULL,
    fpkms = NULL,
    average_fpkms = NULL
  )
  class(value) <- "MyMuse"

  if (run_quant) {
    value$names <- RunKallisto(value, cleanup_index = cleanup_index)
  } else {
    value$names <- infer_sample_names(fastqfiles_dir, paired_end)
  }

  value$tpms <- GetTPMs(value)
  value$eff_lengths <- Getefflengths(value)
  value$est_counts <- GetFinalCounts(value)

  if (isTRUE(value$paired_end)) {
    value$fpkms <- GetFPKMs(value)
    value$average_fpkms <- GetAverageFPKMS(value)
  } else {
    value$rpkms <- GetRPKMs(value)
    value$average_rpkms <- GetAverageRPKMS(value)
  }

  value
}

#' Run kallisto index + quant
#' @param MuseObject A `MyMuse` object.
#' @param cleanup_index Logical; delete index file after quantification.
#' @return Character vector of sample names.
#' @export
RunKallisto <- function(MuseObject, cleanup_index = TRUE) {
  build_kallisto_index(
    kallisto_bin = MuseObject$kallisto_bin,
    index_file = MuseObject$index_file,
    cdsfile = MuseObject$cdsfile
  )

  if (isTRUE(MuseObject$paired_end)) {
    sample_pairs <- discover_paired_fastqs(MuseObject$fastqfiles_dir)
    sample_names <- names(sample_pairs)
    run_kallisto_paired(MuseObject, sample_pairs)
  } else {
    sample_files <- discover_single_fastqs(MuseObject$fastqfiles_dir)
    sample_names <- remove_fastq_extensions(basename(sample_files))
    run_kallisto_single(MuseObject, sample_files, sample_names)
  }

  if (isTRUE(cleanup_index) && file.exists(MuseObject$index_file)) {
    file.remove(MuseObject$index_file)
  }

  sample_names
}

#' Read TSV with stable defaults
#' @param x File path.
#' @param sep Delimiter.
#' @return data.frame
#' @export
read_delim <- function(x, sep = "\t") {
  read.delim(x, sep = sep, header = TRUE, stringsAsFactors = FALSE)
}

#' Collect TPM matrix from kallisto output directories
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and sample columns.
#' @export
GetTPMs <- function(MuseObject) {
  collect_abundance_field(MuseObject, "tpm")
}

#' Collect effective length matrix from kallisto output directories
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and sample columns.
#' @export
Getefflengths <- function(MuseObject) {
  collect_abundance_field(MuseObject, "eff_length")
}

#' Collect estimated count matrix from kallisto output directories
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and sample columns.
#' @export
GetFinalCounts <- function(MuseObject) {
  collect_abundance_field(MuseObject, "est_counts")
}

#' Compute RPKM matrix (single-end)
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and RPKM columns.
#' @export
GetRPKMs <- function(MuseObject) {
  counts <- as.matrix(MuseObject$est_counts[, setdiff(names(MuseObject$est_counts), "target_id")])
  eff_lengths <- MuseObject$eff_lengths[[2]]
  col_sums <- colSums(counts)
  rpkm <- sweep(counts, 2, col_sums / 1e6, "/")
  rpkm <- sweep(rpkm, 1, eff_lengths / 1e3, "/")
  data.frame(target_id = MuseObject$est_counts$target_id, as.data.frame(rpkm), check.names = FALSE)
}

#' Compute per-condition average RPKM matrix
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and averaged columns.
#' @export
GetAverageRPKMS <- function(MuseObject) {
  average_by_replicates(MuseObject$rpkms, MuseObject$reps)
}

#' Compute FPKM matrix (paired-end)
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and FPKM columns.
#' @export
GetFPKMs <- function(MuseObject) {
  counts <- as.matrix(MuseObject$est_counts[, setdiff(names(MuseObject$est_counts), "target_id")])
  eff_lengths <- MuseObject$eff_lengths[[2]]
  col_sums <- colSums(counts)
  fpkm <- sweep(counts, 2, col_sums / 1e6, "/")
  fpkm <- sweep(fpkm, 1, eff_lengths / 1e3, "/")
  data.frame(target_id = MuseObject$est_counts$target_id, as.data.frame(fpkm), check.names = FALSE)
}

#' Compute per-condition average FPKM matrix
#' @param MuseObject A `MyMuse` object.
#' @return data.frame with `target_id` and averaged columns.
#' @export
GetAverageFPKMS <- function(MuseObject) {
  average_by_replicates(MuseObject$fpkms, MuseObject$reps)
}

#' Reorder sample names
#' @param MuseObject A `MyMuse` object.
#' @param order Integer index order.
#' @return Character vector.
#' @export
reorder_names <- function(MuseObject, order) {
  MuseObject$names[order]
}

#' Location where users can place example FASTA/FASTQ files
#' @return Character path inside installed package.
#' @export
kallisto_test_data_dir <- function() {
  system.file("extdata", "kallisto", package = "NatTuretzek")
}

infer_sample_names <- function(fastq_dir, paired_end) {
  if (isTRUE(paired_end)) {
    names(discover_paired_fastqs(fastq_dir))
  } else {
    remove_fastq_extensions(basename(discover_single_fastqs(fastq_dir)))
  }
}

build_kallisto_index <- function(kallisto_bin, index_file, cdsfile) {
  cmd <- paste(shQuote(kallisto_bin), "index", "-i", shQuote(index_file), shQuote(cdsfile))
  status <- system(cmd)
  if (status != 0) stop("kallisto index failed", call. = FALSE)
}

run_kallisto_single <- function(MuseObject, sample_files, sample_names) {
  if (is.null(MuseObject$fragment_length) || is.null(MuseObject$sd)) {
    stop("fragment_length and sd are required for single-end runs.", call. = FALSE)
  }

  for (i in seq_along(sample_files)) {
    out_dir <- sample_names[[i]]
    cmd <- paste(
      shQuote(MuseObject$kallisto_bin), "quant",
      "-i", shQuote(MuseObject$index_file),
      "-o", shQuote(out_dir),
      "--single",
      "-l", MuseObject$fragment_length,
      "-s", MuseObject$sd,
      "-b", MuseObject$bootstrap_samples,
      shQuote(sample_files[[i]])
    )
    status <- system(cmd)
    if (status != 0) stop(sprintf("kallisto quant failed for sample %s", out_dir), call. = FALSE)
  }
}

run_kallisto_paired <- function(MuseObject, sample_pairs) {
  for (sample_name in names(sample_pairs)) {
    pair <- sample_pairs[[sample_name]]
    cmd <- paste(
      shQuote(MuseObject$kallisto_bin), "quant",
      "-i", shQuote(MuseObject$index_file),
      "-o", shQuote(sample_name),
      "-b", MuseObject$bootstrap_samples,
      shQuote(pair[[1]]), shQuote(pair[[2]])
    )
    status <- system(cmd)
    if (status != 0) stop(sprintf("kallisto quant failed for sample %s", sample_name), call. = FALSE)
  }
}

discover_single_fastqs <- function(fastq_dir) {
  files <- list.files(fastq_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No FASTQ files found.", call. = FALSE)
  sort(files)
}

discover_paired_fastqs <- function(fastq_dir) {
  files <- discover_single_fastqs(fastq_dir)
  base <- basename(files)

  cleaned <- gsub("(_R?[12]|_[12])(?=\\.(fastq|fq)(\\.gz)?$)", "", base, perl = TRUE, ignore.case = TRUE)
  groups <- split(files, cleaned)
  groups <- groups[lengths(groups) == 2]

  if (length(groups) == 0) {
    stop("Could not infer paired-end FASTQ pairs from filenames.", call. = FALSE)
  }

  lapply(groups, sort)
}

remove_fastq_extensions <- function(x) {
  gsub("\\.(fastq|fq)(\\.gz)?$", "", x, ignore.case = TRUE)
}

collect_abundance_field <- function(MuseObject, field) {
  tables <- lapply(MuseObject$names, function(sample_name) {
    file <- file.path(sample_name, "abundance.tsv")
    if (!file.exists(file)) stop(sprintf("Missing abundance file: %s", file), call. = FALSE)
    read_delim(file)
  })

  target_id <- tables[[1]][["target_id"]]
  values <- lapply(tables, function(df) df[[field]])
  names(values) <- MuseObject$names
  values[["target_id"]] <- target_id
  data.frame(values, check.names = FALSE)
}

average_by_replicates <- function(expression_df, replicates) {
  stopifnot(replicates >= 1)
  ids <- expression_df$target_id
  mat <- as.matrix(expression_df[, setdiff(names(expression_df), "target_id")])

  group_starts <- seq(1, ncol(mat), by = replicates)
  grouped <- lapply(group_starts, function(start_col) {
    end_col <- min(start_col + replicates - 1, ncol(mat))
    rowMeans(mat[, start_col:end_col, drop = FALSE])
  })

  out <- as.data.frame(grouped)
  names(out) <- paste0("group_", seq_along(grouped))
  data.frame(target_id = ids, out, check.names = FALSE)
}
