#!/usr/bin/env bash
set -euo pipefail

# SPORTS upstream workflow template.
# Run SPORTS for each FASTQ sample first, then use sports_summary_to_counts.R
# and general_count_matrix_de.R for R differential expression.

# Edit these settings before running.
sports_bin="sports.pl"
fastq_dir="your_trimmed_fastq_directory"
output_dir="your_sports_output_directory"
threads="your_thread_number"
mismatches="your_mismatch_number"

genome_index="your_genome_bowtie_index"
mirna_index="your_miRNA_bowtie_index"
rrna_index="your_rRNA_bowtie_index"
trna_index="your_tRNA_bowtie_index"
pirna_index="your_piRNA_bowtie_index"
ensembl_index="your_ensembl_ncRNA_bowtie_index"
rfam_index="your_rfam_bowtie_index"

mkdir -p "$output_dir"

for fastq in "$fastq_dir"/*.fastq "$fastq_dir"/*.fastq.gz "$fastq_dir"/*.fastqsanger.gz; do
  [[ -e "$fastq" ]] || continue

  sample="$(basename "$fastq")"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fastqsanger.gz}"
  sample="${sample%.fastq}"

  sample_out="$output_dir/$sample"
  mkdir -p "$sample_out"

  "$sports_bin" \
    -i "$fastq" \
    -p "$threads" \
    -M "$mismatches" \
    -k \
    -g "$genome_index" \
    -m "$mirna_index" \
    -r "$rrna_index" \
    -t "$trna_index" \
    -w "$pirna_index" \
    -e "$ensembl_index" \
    -f "$rfam_index" \
    -o "$sample_out" \
    > "$sample_out/${sample}_SPORTS.log" 2>&1
done

# After SPORTS finishes:
# 1. Copy or point summary_dir in scripts/sports_summary_to_counts.R to the
#    directory containing all *_summary.txt files.
# 2. Run the R count-matrix script.
# 3. Run scripts/general_count_matrix_de.R separately for miRNA, tsRNA, and rsRNA.
#
# Example:
# Rscript scripts/sports_summary_to_counts.R
# Rscript scripts/general_count_matrix_de.R
