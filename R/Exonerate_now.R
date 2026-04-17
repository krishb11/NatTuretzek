library("dplyr")
library("readr")

# Resolve executable path with optional user override.
resolve_binary <- function(binary_name, user_path = NULL, arg_name = NULL) {
  if (!is.null(user_path) && nzchar(user_path)) {
    if (file.exists(user_path) || nzchar(Sys.which(user_path))) {
      return(user_path)
    }
    stop(sprintf("Provided path for '%s' is not valid: %s", binary_name, user_path), call. = FALSE)
  }

  found <- Sys.which(binary_name)
  if (nzchar(found)) {
    return(found)
  }

  stop(
    sprintf(
      "Could not find '%s' on PATH. Please provide `%s` explicitly.",
      binary_name,
      if (is.null(arg_name)) paste0(binary_name, "_path") else arg_name
    ),
    call. = FALSE
  )
}

#' Prepare BLAST and Exonerate inputs
#'
#' Constructor helper that runs BLAST, splits transcript/contig mappings,
#' generates per-sequence FASTA files, and writes shell scripts for exonerate.
#'
#' @param blasttype BLAST mode, e.g. `"blastn"`, `"megablast"`, `"tblastn"`.
#' @param pathtogenome Path to target genome FASTA.
#' @param pathtoqueryfile Path to BLAST query FASTA.
#' @param species Species label used in output filenames.
#' @param pathtocdsfile Optional CDS FASTA path for exonerate query.
#' @param out_dir Output directory for generated files.
#' @param blast_path Optional BLAST executable path override.
#' @param exonerate_path Optional exonerate executable path override.
#' @param blast_params Named list of BLAST options.
#' @param exonerate_params Named list of exonerate options.
#'
#' @return `MyWitness` object.
#' @export
MyWitness <- function(
  blasttype,
  pathtogenome,
  pathtoqueryfile,
  species,
  pathtocdsfile = NULL,
  out_dir = ".",
  blast_path = NULL,
  exonerate_path = NULL,
  blast_params = list(),
  exonerate_params = list()
) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for MyWitness().", call. = FALSE)
  }
  stopifnot(file.exists(pathtogenome), file.exists(pathtoqueryfile))

  cds_path <- if (is.null(pathtocdsfile)) pathtoqueryfile else pathtocdsfile

  blast_defaults <- list(
    outfmt = 10,
    max_target_seqs = 1,
    num_threads = 10,
    split = 2000,
    additional_args = NULL
  )
  exonerate_defaults <- list(
    model = "e2g",
    softmaskquery = "yes",
    softmasktarget = "yes",
    bestn = 1,
    minintron = 20,
    maxintron = 20000,
    showalignment = "false",
    showtargetgff = TRUE,
    gfffile = "exonerate_output.gff",
    additional_args = NULL
  )

  blast_cfg <- modifyList(blast_defaults, blast_params)
  exonerate_cfg <- modifyList(exonerate_defaults, exonerate_params)

  blast_exec_name <- if (blasttype %in% c("tblastx", "tblastn")) blasttype else "blastn"

  value <- list(
    blast = blasttype,
    genome = pathtogenome,
    cdsfile = cds_path,
    species = species,
    query = pathtoqueryfile,
    out_dir = out_dir,
    csv_dir = file.path(out_dir, "csvfiles"),
    fasta_dir = file.path(out_dir, "fastafiles"),
    blast_bin = resolve_binary(blast_exec_name, blast_path, "blast_path"),
    exonerate_bin = resolve_binary("exonerate", exonerate_path, "exonerate_path"),
    blast_params = blast_cfg,
    exonerate_params = exonerate_cfg
  )
  class(value) <- "MyWitness"

  dir.create(value$csv_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(value$fasta_dir, showWarnings = FALSE, recursive = TRUE)

  runBlast(
    value,
    split = value$blast_params$split,
    additional_args = value$blast_params$additional_args
  )
  GetDataForExonerate(value)
  ExonerateNow(value)
  value
}

#' Run BLAST and split transcript-contig mappings
#' @param WitnessObject A `MyWitness` object.
#' @param split Number of rows per mapping chunk.
#' @param additional_args Additional raw CLI arguments appended to BLAST command.
#' @return Invisibly returns list of generated CSV chunk files.
#' @export
runBlast <- function(WitnessObject, split = 2000, additional_args = NULL) {
  index_csv <- file.path(WitnessObject$csv_dir, "index.csv")
  blast_cfg <- WitnessObject$blast_params

  additional_args <- if (is.null(additional_args)) blast_cfg$additional_args else additional_args

  cmd_parts <- c()
  if (WitnessObject$blast %in% c("tblastx", "tblastn")) {
    cmd_parts <- c(
      shQuote(WitnessObject$blast_bin),
      "-subject", shQuote(WitnessObject$genome),
      "-query", shQuote(WitnessObject$query),
      "-outfmt", blast_cfg$outfmt,
      "-out", shQuote(index_csv),
      "-max_target_seqs", blast_cfg$max_target_seqs,
      "-num_threads", blast_cfg$num_threads
    )
  } else {
    cmd_parts <- c(
      shQuote(WitnessObject$blast_bin),
      "-task", WitnessObject$blast,
      "-subject", shQuote(WitnessObject$genome),
      "-query", shQuote(WitnessObject$query),
      "-outfmt", blast_cfg$outfmt,
      "-out", shQuote(index_csv),
      "-max_target_seqs", blast_cfg$max_target_seqs,
      "-num_threads", blast_cfg$num_threads
    )
  }

  if (!is.null(additional_args) && nzchar(additional_args)) {
    cmd_parts <- c(cmd_parts, additional_args)
  }

  cmd <- paste(cmd_parts, collapse = " ")
  status <- system(cmd)
  if (status != 0) stop("BLAST command failed.", call. = FALSE)

  index <- utils::read.csv(index_csv, header = FALSE, stringsAsFactors = FALSE)
  index <- unique(index[, c(1, 2)])
  names(index) <- c("transcript", "contig")

  index$chunk <- ceiling(seq_len(nrow(index)) / split)
  chunk_files <- index %>%
    group_by(.data$chunk) %>%
    group_map(~ {
      path <- file.path(WitnessObject$csv_dir, paste0(WitnessObject$species, .y$chunk, ".csv"))
      readr::write_csv(.x[, c("transcript", "contig")], path)
      path
    }) %>%
    unlist(use.names = FALSE)

  invisible(chunk_files)
}

#' Create sequence-level FASTA files for exonerate
#' @param WitnessObject A `MyWitness` object.
#' @return Invisibly returns output FASTA file paths.
#' @export
GetDataForExonerate <- function(WitnessObject) {
  genome <- Biostrings::readDNAStringSet(WitnessObject$genome)
  cdsfiles <- Biostrings::readDNAStringSet(WitnessObject$cdsfile)

  write_one_fasta <- function(seq_set) {
    lapply(seq_along(seq_set), function(i) {
      name <- names(seq_set)[[i]]
      out <- file.path(WitnessObject$fasta_dir, paste0(name, ".fasta"))
      Biostrings::writeXStringSet(seq_set[i], out)
      out
    })
  }

  files <- c(write_one_fasta(genome), write_one_fasta(cdsfiles))
  invisible(unlist(files, use.names = FALSE))
}

#' List CSV mapping chunks created for exonerate
#' @param WitnessObject Optional `MyWitness` object.
#' @param csv_dir Optional fallback directory.
#' @return Character vector of CSV filenames.
#' @export
GatherCSV <- function(WitnessObject = NULL, csv_dir = "./csvfiles") {
  target_dir <- if (!is.null(WitnessObject)) WitnessObject$csv_dir else csv_dir
  list.files(path = target_dir, pattern = "[0-9]\\.csv$", full.names = TRUE)
}

#' Generate exonerate shell scripts
#' @param WitnessObject A `MyWitness` object.
#' @param gfffile Optional output GFF filename override.
#' @param additional_args Optional raw exonerate arguments appended to command.
#' @return Invisibly returns generated shell script paths.
#' @export
ExonerateNow <- function(WitnessObject, gfffile = NULL, additional_args = NULL) {
  csv_files <- GatherCSV(WitnessObject)
  ex_cfg <- WitnessObject$exonerate_params

  if (is.null(gfffile)) {
    gfffile <- ex_cfg$gfffile
  }
  if (is.null(additional_args)) {
    additional_args <- ex_cfg$additional_args
  }

  scripts <- lapply(csv_files, function(csv_path) {
    data <- utils::read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

    commands <- vapply(seq_len(nrow(data)), function(i) {
      gfile <- file.path(WitnessObject$fasta_dir, paste0(data$contig[[i]], ".fasta"))
      cfile <- file.path(WitnessObject$fasta_dir, paste0(data$transcript[[i]], ".fasta"))

      cmd_parts <- c(
        shQuote(WitnessObject$exonerate_bin),
        "--model", ex_cfg$model,
        "-q", shQuote(cfile),
        "-t", shQuote(gfile),
        "--softmaskquery", ex_cfg$softmaskquery,
        "--softmasktarget", ex_cfg$softmasktarget,
        "--bestn", ex_cfg$bestn,
        "--minintron", ex_cfg$minintron,
        "--maxintron", ex_cfg$maxintron,
        "--showalignment", ex_cfg$showalignment
      )

      if (isTRUE(ex_cfg$showtargetgff)) {
        cmd_parts <- c(cmd_parts, "--showtargetgff")
      }

      if (!is.null(additional_args) && nzchar(additional_args)) {
        cmd_parts <- c(cmd_parts, additional_args)
      }

      cmd_parts <- c(cmd_parts, ">>", shQuote(file.path(WitnessObject$out_dir, gfffile)))
      paste(cmd_parts, collapse = " ")
    }, character(1))

    script_name <- sub("\\.csv$", ".sh", basename(csv_path))
    script_path <- file.path(WitnessObject$out_dir, script_name)
    writeLines(commands, script_path)
    script_path
  })

  invisible(unlist(scripts, use.names = FALSE))
}

#' Remove generated temporary files
#' @param WitnessObject Optional `MyWitness` object.
#' @param remove_csv Remove csv directory.
#' @param remove_fasta Remove fasta directory.
#' @return Invisibly returns deleted paths.
#' @export
RemoveFiles <- function(WitnessObject = NULL, remove_csv = TRUE, remove_fasta = TRUE) {
  csv_dir <- if (!is.null(WitnessObject)) WitnessObject$csv_dir else "csvfiles"
  fasta_dir <- if (!is.null(WitnessObject)) WitnessObject$fasta_dir else "fastafiles"

  deleted <- character(0)
  if (remove_fasta && dir.exists(fasta_dir)) {
    unlink(fasta_dir, recursive = TRUE, force = TRUE)
    deleted <- c(deleted, fasta_dir)
  }
  if (remove_csv && dir.exists(csv_dir)) {
    unlink(csv_dir, recursive = TRUE, force = TRUE)
    deleted <- c(deleted, csv_dir)
  }

  invisible(deleted)
}
