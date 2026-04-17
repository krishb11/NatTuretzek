library("dplyr")
library("readr")

#' Parse exonerate GFF and extract transcript FASTA sequences
#'
#' @param file Path to GFF file.
#' @param species Species label for naming output FASTA.
#' @param style GFF style; currently supports `"exonerate"`.
#' @param genomefile Path to genome FASTA used for coordinate extraction.
#'
#' @return `MyDeep` object.
#' @export
MyDeep <- function(file, species, style = "exonerate", genomefile) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for MyDeep().", call. = FALSE)
  }
  stopifnot(file.exists(file), file.exists(genomefile))

  value <- list(
    pathtofile = file,
    species = species,
    data = NULL,
    edited_gff_object = list(),
    style = style,
    fasta = NULL,
    genomefile = genomefile
  )
  class(value) <- "MyDeep"

  value$data <- read_gff(value)
  value$edited_gff_object <- gffeditor(value)
  value$fasta <- gff2fastaprinter(DeepObject = value)
  value
}

#' Read and normalize GFF file
#' @param DeepObject A `MyDeep` object.
#' @param pathtofile Optional direct path; overrides `DeepObject$pathtofile`.
#' @param filter_exonerate Logical; filter to `exonerate:est2genome` lines.
#' @return data.frame with canonical GFF columns.
#' @export
read_gff <- function(DeepObject, pathtofile = NULL, filter_exonerate = TRUE) {
  file <- if (is.null(pathtofile)) DeepObject$pathtofile else pathtofile

  if (isTRUE(filter_exonerate)) {
    lines <- readLines(file, warn = FALSE)
    lines <- lines[grepl("exonerate:est2genome", lines, ignore.case = TRUE)]
    lines <- lines[!grepl("source-version", lines, fixed = TRUE)]
    if (length(lines) == 0) {
      stop("No exonerate:est2genome records found in GFF.", call. = FALSE)
    }
    tf <- tempfile(fileext = ".gff")
    writeLines(lines, tf)
    file <- tf
  }

  gff_input <- suppressWarnings(readr::read_delim(file, delim = "\t", col_names = FALSE, comment = "#", show_col_types = FALSE))
  if (ncol(gff_input) > 9) {
    stop("GFF format supports up to 9 columns.", call. = FALSE)
  }

  gff_names <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
  names(gff_input)[seq_len(ncol(gff_input))] <- gff_names[seq_len(ncol(gff_input))]
  as.data.frame(gff_input)
}

#' Split GFF into per-gene exon tables
#' @param DeepObject A `MyDeep` object.
#' @return Named list of exon-only data.frames grouped by gene.
#' @export
gffeditor <- function(DeepObject) {
  if (!identical(DeepObject$style, "exonerate")) {
    stop("Only style = 'exonerate' is currently supported.", call. = FALSE)
  }

  gene_rows <- which(DeepObject$data$type == "gene")
  if (length(gene_rows) == 0) stop("No gene rows found in GFF.", call. = FALSE)

  chunks <- lapply(seq_along(gene_rows), function(i) {
    start <- gene_rows[[i]]
    end <- if (i == length(gene_rows)) nrow(DeepObject$data) else gene_rows[[i + 1]] - 1
    DeepObject$data[start:end, , drop = FALSE]
  })

  gene_table <- DeepObject$data %>% filter(.data$type == "gene")
  gene_names <- vapply(gene_table$attribute, function(att) {
    parts <- strsplit(att, ";")[[1]]
    if (length(parts) < 2) return(att)
    sub(".*\\s", "", trimws(parts[[2]]))
  }, character(1))

  names(chunks) <- make.unique(gene_names)
  lapply(chunks, function(df) df %>% filter(.data$type == "exon"))
}

#' Extract transcript FASTA sequences from exon coordinates
#' @param DeepObject A `MyDeep` object.
#' @param edited_gff_object Optional list returned by [gffeditor()].
#' @param output_file Optional output FASTA path.
#' @return DNAStringSet of assembled transcript sequences.
#' @export
gff2fastaprinter <- function(DeepObject, edited_gff_object = NULL, output_file = NULL) {
  exon_list <- if (is.null(edited_gff_object)) DeepObject$edited_gff_object else edited_gff_object

  genome <- Biostrings::readDNAStringSet(DeepObject$genomefile)
  names(genome) <- vapply(strsplit(names(genome), " "), `[`, character(1), 1)

  seqs <- lapply(exon_list, function(exons) {
    if (nrow(exons) == 0) return(NA_character_)
    contig <- exons$seqid[[1]]
    ref <- genome[names(genome) == contig]
    if (length(ref) == 0) return(NA_character_)

    pieces <- vapply(seq_len(nrow(exons)), function(i) {
      start <- exons$start[[i]]
      end <- exons$end[[i]]
      if (start > Biostrings::width(ref) || end > Biostrings::width(ref)) return("")
      as.character(Biostrings::subseq(ref, start = start, end = end))
    }, character(1))

    paste0(pieces, collapse = "")
  })

  seqs <- unlist(seqs, use.names = TRUE)
  seqs <- seqs[!is.na(seqs)]
  dat <- Biostrings::DNAStringSet(seqs)

  if (is.null(output_file)) {
    output_file <- paste(DeepObject$species, "transcripts", "fasta", sep = ".")
  }
  Biostrings::writeXStringSet(dat, output_file)
  dat
}
