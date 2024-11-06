#!/bin/bash

#this script extract 200 or exact nt size based on for_search.txt file
#for "+" it convert it to reverse sequence
#for "-" it convert it to compliment sequence for normal alignment with neddle


export PATH="/User/bedops/bin:$PATH"


# Define input files and variables
GFF_FILE=""   # Path to the GFF file
GENOME_FASTA=""  # Path to the genome FASTA file
CHROM_SIZES="/Users/chrom.sizes"    # Path to chromosome sizes file
# Path to txt file containing data about PROTEIN_ID,	STRAND,	SIZE to extract:
INPUT_FILE="" #in tab separated .txt


# Process each line of the input file
while IFS=$'\t' read -r GENE_ID STRAND SIZE; do
  if [[ "$GENE_ID" == "GENE_ID" ]]; then
    continue  # Skip header
  fi

  # Create the BED file for the current protein
  BED_FILE="${GENE_ID}.bed"

  # Adjust search length based on size using arithmetic expansion $(( ))
  if (( SIZE > 200 )); then
    SEARCH_LENGTH=200
  else
    SEARCH_LENGTH=$SIZE
  fi

  # Extract upstream sequence based on strand
  if [[ "$STRAND" == "+" ]]; then #for forward sequence from NCBI fasta
    bedtools flank -i "$BED_FILE" -g "$CHROM_SIZES" -l "$SEARCH_LENGTH" -r 0 | \
    bedtools getfasta -fi "$GENOME_FASTA" -bed - -fo "${GENE_ID}_${SEARCH_LENGTH}.fasta"


  elif [[ "$STRAND" == "-" ]]; then #for reverse sequense extract and convert to complement
    bedtools flank -i "$BED_FILE" -g "$CHROM_SIZES" -l "$SEARCH_LENGTH" -r 0 -s | \
    bedtools getfasta -fi "$GENOME_FASTA" -bed - -fo "${GENE_ID}_${SEARCH_LENGTH}_tmp.fasta"

    # Convert the sequence to complement (A <-> T, C <-> G)
    # Complement the sequence (without complementing headers)
    awk '/^>/ {print $0} !/^>/ {print}' "${GENE_ID}_${SEARCH_LENGTH}_tmp.fasta" | \
    tr "ATGCatgc" "TACGtacg" > "${GENE_ID}_${SEARCH_LENGTH}.fasta"

    # Clean up temporary file
    rm "${GENE_ID}_${SEARCH_LENGTH}_tmp.fasta"
  else
    echo "Unknown strand: $STRAND for protein ID: $GENE_ID"
    continue
  fi

  echo "Processed $GENE_ID with strand $STRAND and size $SIZE"

done < "$INPUT_FILE"
