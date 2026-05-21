#!/usr/bin/env bash
set -euo pipefail

# proTRAC upstream workflow template.
# Run proTRAC per sample, collect cluster GTF files, build a shared cluster
# count matrix, then use protrac_cluster_de.R for R differential expression.

# Edit these settings before running.
protrac_bin="proTRAC.pl"
sam_dir="your_sam_or_bam_directory"
output_dir="your_protrac_output_directory"
reference_genome="your_reference_genome.fa"
threads="your_thread_number"
pirna_min_length="your_piRNA_min_length"
pirna_max_length="your_piRNA_max_length"
first_base_or_tenth_base_bias="your_1T_or_10A_bias_threshold"
cluster_strand_bias="your_cluster_strand_bias_threshold"

merged_clusters_gtf="your_combined_piRNA_clusters.gtf"
merged_clusters_bed="your_combined_piRNA_clusters.bed"
count_matrix="your_piRNA_cluster_count_matrix.txt"

mkdir -p "$output_dir"

for alignment in "$sam_dir"/*.sam "$sam_dir"/*.bam; do
  [[ -e "$alignment" ]] || continue

  sample="$(basename "$alignment")"
  sample="${sample%.sam}"
  sample="${sample%.bam}"
  sample_out="$output_dir/proTRAC_${sample}"
  mkdir -p "$sample_out"

  # Adjust arguments to match your proTRAC installation and documentation.
  "$protrac_bin" \
    -map "$alignment" \
    -genome "$reference_genome" \
    -pimin "$pirna_min_length" \
    -pimax "$pirna_max_length" \
    -1Tor10A "$first_base_or_tenth_base_bias" \
    -clstrand "$cluster_strand_bias" \
    -o "$sample_out" \
    -p "$threads" \
    > "$sample_out/${sample}_proTRAC.log" 2>&1
done

# Merge all proTRAC cluster GTF files into one cluster set.
find "$output_dir" -name "clusters.gtf" -print0 | xargs -0 cat > "$merged_clusters_gtf"

# Add gene_id to proTRAC cluster attributes.
awk 'BEGIN{FS=OFS="\t"} {
  gsub(/.*piRNA cluster no: ([^;]+);.*/, "gene_id \"cluster_"$9"\";", $9)
  print
}' "$merged_clusters_gtf" > "${merged_clusters_gtf%.gtf}.modified.gtf"

# Convert merged GTF to BED6.
awk 'BEGIN{FS=OFS="\t"} $3=="piRNA_cluster" {print $1, $4-1, $5, $9, ".", $7}' \
  "${merged_clusters_gtf%.gtf}.modified.gtf" \
  | sort -k1,1 -k2,2n \
  > "$merged_clusters_bed"

# Count overlaps between each sample alignment file and the shared cluster BED.
# For BAM files, bedtools multicov can create a cluster-by-sample count table.
# If your inputs are SAM, convert to sorted/indexed BAM first.
bam_files=("$sam_dir"/*.bam)
bedtools multicov -bams "${bam_files[@]}" -bed "$merged_clusters_bed" > "${count_matrix%.txt}.raw.tsv"

# Convert BED + count columns to a DESeq2-friendly count matrix.
sample_names=""
for bam in "${bam_files[@]}"; do
  sample="$(basename "$bam")"
  sample="${sample%.bam}"
  sample_names="${sample_names}${sample_names:+,}${sample}"
done

awk -v sample_names="$sample_names" 'BEGIN{FS=OFS="\t";
  sample_count=split(sample_names, samples, ",")
  printf "cluster"
  for (i=1; i<=sample_count; i++) printf OFS samples[i]
  printf ORS
}
{
  cluster=$1":"$2"-"$3":"$6
  printf cluster
  for (i=7; i<=NF; i++) printf OFS $i
  printf ORS
}' "${count_matrix%.txt}.raw.tsv" > "$count_matrix"

# After the count matrix is created, edit scripts/protrac_cluster_de.R so
# count_file points to this matrix, then run:
# Rscript scripts/protrac_cluster_de.R
