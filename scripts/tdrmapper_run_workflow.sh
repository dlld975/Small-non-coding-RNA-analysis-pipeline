#!/usr/bin/env bash
set -euo pipefail

# tDRMapper upstream workflow template.
# Run tDRMapper for each FASTQ sample first, build a tRNA count matrix, then
# use tdrmapper_de.R for R differential expression.

# Edit these settings before running.
tdrmapper_script="TdrMappingScripts.pl"
trna_reference="your_tRNA_reference.fa"
fastq_dir="your_trimmed_fastq_directory"
output_dir="your_tdrmapper_output_directory"
count_matrix="your_tRNA_count_matrix.txt"
tdr_feature_column="your_tRNA_ID_column_number"
tdr_count_column="your_read_count_column_number"

mkdir -p "$output_dir"
cd "$output_dir"

for fastq in "$fastq_dir"/*.fastq "$fastq_dir"/*.fastq.gz "$fastq_dir"/*.fastqsanger.gz; do
  [[ -e "$fastq" ]] || continue

  sample="$(basename "$fastq")"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fastqsanger.gz}"
  sample="${sample%.fastq}"

  mkdir -p "$sample"
  (
    cd "$sample"
    perl "$tdrmapper_script" "$trna_reference" "$fastq" > "${sample}_tDRMapper.log" 2>&1
  )
done

# Build a count matrix from tDRMapper speciesInfo outputs.
# Set tdr_feature_column to the tRNA ID column and tdr_count_column to the
# read-count column in your tDRMapper speciesInfo files.
find "$output_dir" -name "*.hq_cs.mapped.speciesInfo.txt" -print > speciesInfo_files.txt

awk -v feature_col="$tdr_feature_column" -v count_col="$tdr_count_column" 'BEGIN{FS=OFS="\t"}
FNR==1 {
  sample=FILENAME
  sub(/.*\//, "", sample)
  sub(/\.hq_cs\.mapped\.speciesInfo\.txt$/, "", sample)
  samples[sample]=1
}
{
  feature=$feature_col
  count=$count_col
  if (feature != "" && count ~ /^[0-9.]+$/) {
    counts[feature,sample] += count
    features[feature]=1
  }
}
END{
  printf "feature"
  for (s in samples) printf OFS s
  printf ORS

  for (f in features) {
    printf f
    for (s in samples) printf OFS (counts[f,s] == "" ? 0 : counts[f,s])
    printf ORS
  }
}' $(cat speciesInfo_files.txt) > "$count_matrix"

# After the count matrix is created, edit scripts/tdrmapper_de.R so count_file
# points to this matrix, then run:
# Rscript scripts/tdrmapper_de.R
