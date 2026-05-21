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
# This template assumes tRNA ID is in column 2 and read count is in column 4,
# matching the workflow notes. Check your tDRMapper output before running.
find "$output_dir" -name "*.hq_cs.mapped.speciesInfo.txt" -print > speciesInfo_files.txt

awk 'BEGIN{FS=OFS="\t"}
FNR==1 {
  sample=FILENAME
  sub(/.*\//, "", sample)
  sub(/\.hq_cs\.mapped\.speciesInfo\.txt$/, "", sample)
  samples[sample]=1
}
{
  feature=$2
  count=$4
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
