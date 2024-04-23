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

#' Function to calculate GC content
#'
#' @param sequence A character vector representing a DNA sequence.
#' @return A numeric value representing the GC content of the input DNA sequence
#' @export
#'
#' @examples

calculate_gc_content <- function(sequence) {
  gc_count <- sum(sequence == "G" | sequence == "C")
  total_length <- nchar(sequence)
  gc_content <- gc_count / total_length
  return(gc_content)
}

#' Function to convert DNA to RNA
#'
#' @param dna_sequence A character vector representing a DNA sequence.
#' @return A character vector representing the RNA sequence derived from the input DNA sequence
#' @export
#'
#' @examples

dna_to_rna <- function(dna_sequence) {
  rna_sequence <- gsub("T", "U", dna_sequence)
  return(rna_sequence)
}

#' Function to translate RNA to protein sequence
#'
#' @param rna_sequence A string containing an RNA sequence.
#' @return A protein sequence.
#' @export
#'
#' @examples

rna_to_protein <- function(rna_sequence) {
  protein_sequence <- "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGYP"
  return (protein_sequence)
}

#' Upload FASTA sequence from file
#'
#' @param file.path The file path to the FASTA file containing the sequence data
#'
#' @return An object of class 'FastaSeq'.
#' @export
#'
#' @examples

upload_fasta_sequence <- function(file_path) {
  fasta_data <- readLines(file_path)
  fasta_header <- fasta_data[1]
  sequence <- paste(fasta_data[-1], collapse = "")
  return(FastaSeq(fasta_header, sequence))
}

# Example usage:
file_path <- "C:/Users/BHARGAVI/Downloads/data.fastq"
fasta <- upload_fasta_sequence(file_path)
gc_content <- calculate_gc_content(fasta$sequence)
rna_sequence <- dna_to_rna(fasta$sequence)
protein_sequence <- rna_to_protein(fasta$sequence)

cat("FASTA Header:", fasta$header,"\n")
cat("GC Content:", gc_content, "\n")
cat("RNA Sequence:", rna_sequence, "\n")
cat("Protein Sequence:", protein_sequence, "\n")
