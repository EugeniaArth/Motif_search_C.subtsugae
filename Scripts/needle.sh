#!/bin/bash

# alignment with neddle


export PATH="/User/bedops/bin:$PATH"


INPUT_FILE=""  # Path to txt file with  data about GENE_ID,	STRAND,	SIZE
MOTIF_FASTA=""  # Path to motif FASTA file which we align with promoter DNA of every gene
OUTPUT_DIR=""  # Directory to save alignment results

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process each line of the input file
while IFS=$'\t' read -r GENE_ID STRAND SIZE; do
  if [[ "$GENE_ID" == "GENE_ID" ]]; then
    continue  # Skip header
  fi

  # Adjust search length based on size using arithmetic expansion $(( ))
  if (( SIZE > 200 )); then
    SEARCH_LENGTH=200
  else
    SEARCH_LENGTH=$SIZE
  fi

  # Define input protein FASTA file and output alignment result file
  GENE_FASTA="${GENE_ID}_${SEARCH_LENGTH}.fasta"  # Adjust this to your actual fasta file names
  OUTPUT_FILE="${OUTPUT_DIR}/${GENE_ID}.txt"

  # Perform alignment using needle
  needle -asequence "$MOTIF_FASTA" -bsequence "$GENE_FASTA" -outfile "$OUTPUT_FILE" -gapopen 100.0 -gapextend 0.5

  echo "Aligned motif.fasta with $GENE_FASTA, result saved to $OUTPUT_FILE"

done < "$INPUT_FILE"
