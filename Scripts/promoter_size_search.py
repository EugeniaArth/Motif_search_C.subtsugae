#This script performs search of direction of gene transcription, calculate the possible size of  upstream zone
#where the promoter can be located and data about possible overlaps. For search we need annotation in gff and
#gene_data.csv with the list of genes (ORF)
#It will add additional columns to left, resulted file name is "Gene_Promoter_Table_with_Overlap.csv"
#It takes into account the strand at which the gene is located - if the gene on "+", it searches the promoter before the gene,
#on the left side from TSS
#If the gene is located on "-" strans, it searches promoter at the rigth side from TSS

import pandas as pd
from collections import defaultdict

# Load the CSV and GFF files
csv_file_path = ''  # CSV file path with gene names in column named ORF (genes_data_for_promoter.csv)
gff_file_path = ''  # GFF file path

# Load the CSV into a DataFrame
df = pd.read_csv(csv_file_path, delimiter=";")

# Initialize columns for transcription direction, promoter size, and overlapping gene
df["Direction of gene transcription"] = ""
df["Size of upstream zone"] = ""
df["Overlap with gene in upstream zone"] = ""

# Load and parse the GFF annotation file to extract gene information
# We will organize genes by sequence ID (chromosome/contig) and strand
genes_by_seqid_strand = defaultdict(list)

with open(gff_file_path, 'r') as gff_file:
    for line in gff_file:
        if line.startswith("#"):
            continue  # Skip comments
        parts = line.strip().split("\t")
        if len(parts) == 9:  # GFF format expects 9 fields
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            attributes_dict = {k: v for k, v in [field.split("=") for field in attributes.strip().split(";") if "=" in field]}
            if feature == "gene":
                gene_info = {
                    'seqid': seqid,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'gene_id': attributes_dict.get("ID", "").replace("gene-", "")
                }
                # Organize genes by seqid and strand
                genes_by_seqid_strand[(seqid, strand)].append(gene_info)

# Sort the genes for each seqid and strand based on their start position
for key in genes_by_seqid_strand:
    genes_by_seqid_strand[key].sort(key=lambda x: x['start'])

# Function to calculate promoter size with correct neighbor identification
def calculate_promoter_size(gene_list, current_gene):
    idx = gene_list.index(current_gene)
    strand = current_gene['strand']
    start = current_gene['start']
    end = current_gene['end']

    if strand == '+':  # Forward strand
        # Check if there is a previous gene
        if idx == 0:
            return "more than 50", None  # No previous gene on this strand
        prev_gene = gene_list[idx - 1]
        prev_end = prev_gene['end']
        # Overlap check
        if start <= prev_end:
            return "0", prev_gene['gene_id']  # Overlap detected, return gene ID
        distance = start - prev_end - 1
    elif strand == '-':  # Reverse strand
        # Check if there is a next gene
        if idx == len(gene_list) - 1:
            return "more than 50", None  # No next gene on this strand
        next_gene = gene_list[idx + 1]
        next_start = next_gene['start']
        # Overlap check
        if end >= next_start:
            return "0", next_gene['gene_id']  # Overlap detected, return gene ID
        distance = next_start - end - 1
    # Check distance thresholds
    if distance < 50:
        return "less than 50", None
    else:
        return distance, None  # Return exact promoter size and no overlap

# Create a mapping from gene IDs to gene info for quick access
gene_id_to_gene_info = {}
for gene_list in genes_by_seqid_strand.values():
    for gene in gene_list:
        gene_id_to_gene_info[gene['gene_id']] = gene

# Iterate over each gene in the CSV file to calculate promoter size and direction
for index, row in df.iterrows():
    orf = row["ORF"]
    matching_gene = gene_id_to_gene_info.get(orf)
    if matching_gene:
        seqid = matching_gene['seqid']
        strand = matching_gene['strand']
        # Get the sorted gene list for this seqid and strand
        gene_list = genes_by_seqid_strand[(seqid, strand)]
        # Fill in transcription direction
        df.at[index, "Direction of gene transcription"] = strand
        # Calculate and fill in promoter size and overlapping gene
        promoter_size, overlapping_gene_id = calculate_promoter_size(gene_list, matching_gene)
        df.at[index, "Size of upstream zone"] = promoter_size
        df.at[index, "Overlap with gene in upstream zone"] = overlapping_gene_id if overlapping_gene_id else ""
    else:
        # If the gene is not found in the GFF data
        df.at[index, "Direction of gene transcription"] = "Not found"
        df.at[index, "Size of upstream zone"] = "Not found"
        df.at[index, "Overlap with gene in upstream zone"] = ""

# Save the updated table to a new CSV file
df.to_csv("Gene_Promoter_Table_with_Overlap.csv", index=False)

print("Analysis complete. Updated table saved as 'Gene_Promoter_Table_with_Overlap.csv'.")
