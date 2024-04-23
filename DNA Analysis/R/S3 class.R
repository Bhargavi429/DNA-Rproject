#' Create S3 class for FASTA sequences
#'
#' @title FASTA Sequence Class
#' @name FastaSeq
#' @aliases FastaSeq-class
#' @slot header The header of the FASTA sequence.
#' @slot sequence The DNA sequence.
#' @export

FastaSeq <- function(header, sequence) {
  class(sequence) <- "FastaSeq"
  list(header = header, sequence = sequence)
}

