#Function to calculate GC content
#' Title Calculate GC Content
#'
#' @param sequence A character vector representing a DNA sequence.
#'
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

#Function to Transcript DNA to RNA
#' Title Transcribe DNA to RNA
#'
#' @param dna_sequence A character vector representing a DNA sequence.
#'
#'
#' @return A character vector representing the RNA sequence derived from the input DNA sequence
#' @export
#'
#' @examples

dna_to_rna <- function(dna_sequence) {
  rna_sequence <- gsub("T", "U", dna_sequence)
  return(rna_sequence)
}

# Function to Translate RNA to protein sequence
#' Title Translate RNA to Protein
#'
#' @param rna_sequence rna_sequence A string containing an RNA sequence.
#'
#' @return A protein sequence.
#' @export
#'
#' @examples
#'
rna_to_protein <- function(rna_sequence) {
  # This is a placeholder for the actual translation logic
  # You would need a codon table and logic to convert codons to amino acids
  protein_sequence <- "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGYP"
  return (protein_sequence)
}

# Upload FASTA sequence from file
#' Title
#'
#' @param file.path The file path to the FASTA file containing the sequence data
#'
#' @return A list with two elements: 'header' containing the FASTA header and 'sequence' containing the concatenated sequence data.
#' @export
#'
#' @examples

upload_fasta_sequence <- function(file.path) {
  fasta_data <- readLines(file_path)
  fasta_header <- fasta_data[1]
  sequence <- paste(fasta_data[-1], collapse = "")
  return(list(header = fasta_header, sequence = sequence))
}

#Example usage:
file_path <- "C:/Users/BHARGAVI/Downloads/data.fastq"
fasta <- upload_fasta_sequence(file_path)
gc_content <- calculate_gc_content(fasta$sequence)
rna_sequence <- dna_to_rna(fasta$sequence)
protein_sequence <- rna_to_protein(fasta$sequence)

cat("FASTA Header:", fasta$header,"\n")
cat("GC Content:", gc_content, "\n")
cat("RNA Sequence:", rna_sequence, "\n")
cat("Protein Sequence:", protein_sequence, "\n")
