#!/usr/bin/env bash
set -euo pipefail

# proTRAC non-DE unique cluster workflow.
# This merges treatment and control proTRAC cluster GTF files separately,
# identifies clusters with no overlap in the opposite group, then intersects
# unique clusters with gene annotations for functional study.

# Edit these settings before running.
workdir="your_protrac_project_directory"
treatment_label="your_TREATMENT"
control_label="CONTROL"
gene_gtf="your_gene_annotation.gtf"
genes_bed="genes.bed"
intersect_gene_name_column="your_gene_name_column_after_bedtools_intersect"

cd "$workdir"

# Replace these glob patterns with the cluster GTF files for each group.
treatment_gtfs=(proTRAC_your_treatment_sample*/clusters.gtf)
control_gtfs=(proTRAC_your_control_sample*/clusters.gtf)

cat "${treatment_gtfs[@]}" > "${treatment_label}combined.gtf"
cat "${control_gtfs[@]}" > "${control_label}combined.gtf"

# Add a simple gene_id attribute to each proTRAC cluster.
awk 'BEGIN{FS=OFS="\t"} {
  gsub(/.*piRNA cluster no: ([^;]+);.*/, "gene_id \"cluster_"$9"\";", $9)
  print
}' "${treatment_label}combined.gtf" > "${treatment_label}modified.gtf"

awk 'BEGIN{FS=OFS="\t"} {
  gsub(/.*piRNA cluster no: ([^;]+);.*/, "gene_id \"cluster_"$9"\";", $9)
  print
}' "${control_label}combined.gtf" > "${control_label}modified.gtf"

# Convert proTRAC piRNA clusters from GTF to BED6.
awk 'BEGIN{FS=OFS="\t"} $3=="piRNA_cluster" {print $1, $4-1, $5, $9, ".", $7}' \
  "${treatment_label}modified.gtf" > "${treatment_label}modified.bed"

awk 'BEGIN{FS=OFS="\t"} $3=="piRNA_cluster" {print $1, $4-1, $5, $9, ".", $7}' \
  "${control_label}modified.gtf" > "${control_label}modified.bed"

sort -k1,1 -k2,2n "${treatment_label}modified.bed" > "${treatment_label}modified.sorted.bed"
sort -k1,1 -k2,2n "${control_label}modified.bed" > "${control_label}modified.sorted.bed"

echo "${treatment_label} chroms:"
cut -f1 "${treatment_label}modified.sorted.bed" | sort -u
echo "${control_label} chroms:"
cut -f1 "${control_label}modified.sorted.bed" | sort -u

echo "Chroms in ${treatment_label} but not in ${control_label}:"
comm -23 \
  <(cut -f1 "${treatment_label}modified.sorted.bed" | sort -u) \
  <(cut -f1 "${control_label}modified.sorted.bed" | sort -u)

echo "Chroms in ${control_label} but not in ${treatment_label}:"
comm -13 \
  <(cut -f1 "${treatment_label}modified.sorted.bed" | sort -u) \
  <(cut -f1 "${control_label}modified.sorted.bed" | sort -u)

# Any overlap between treatment and control clusters.
bedtools intersect \
  -a "${treatment_label}modified.sorted.bed" \
  -b "${control_label}modified.sorted.bed" \
  -wa -wb \
  > "${treatment_label}_vs_${control_label}.overlap.tsv"

# Classify overlapping treatment clusters as WITHIN, INCLUDES, or OVERLAP_ONLY.
awk -v treatment="$treatment_label" -v control="$control_label" 'BEGIN{FS=OFS="\t"}
{
  aS=$2; aE=$3
  bS=$8; bE=$9
  rel="OVERLAP_ONLY"
  if (aS>=bS && aE<=bE) rel=treatment"_WITHIN_"control
  else if (aS<=bS && aE>=bE) rel=treatment"_INCLUDES_"control
  print $1,aS,aE,rel,$4,bS,bE,$10
}' "${treatment_label}_vs_${control_label}.overlap.tsv" \
  > "${treatment_label}_vs_${control_label}.classified.tsv"

# Core non-DE outputs:
# Treatment unique no-overlap clusters can be used for functional study
# as treatment-specific/upregulated-like clusters.
bedtools intersect \
  -a "${treatment_label}modified.sorted.bed" \
  -b "${control_label}modified.sorted.bed" \
  -v \
  > "${treatment_label}_unique_noOverlap.bed"

# Control unique no-overlap clusters can be used for functional study
# as control-specific/downregulated-like clusters.
bedtools intersect \
  -a "${control_label}modified.sorted.bed" \
  -b "${treatment_label}modified.sorted.bed" \
  -v \
  > "${control_label}_unique_noOverlap.bed"

awk -v rel="${treatment_label}_WITHIN_${control_label}" '$4==rel' \
  "${treatment_label}_vs_${control_label}.classified.tsv" \
  > "${treatment_label}_within_${control_label}.tsv"

awk -v rel="${treatment_label}_INCLUDES_${control_label}" '$4==rel' \
  "${treatment_label}_vs_${control_label}.classified.tsv" \
  > "${treatment_label}_includes_${control_label}.tsv"

# Build BED6 gene annotation from GTF.
awk 'BEGIN{FS=OFS="\t"} $3=="gene" {
  g=$9
  sub(/.*gene_name "/, "", g)
  sub(/".*/, "", g)
  print $1, $4-1, $5, g, ".", $7
}' "$gene_gtf" | sort -k1,1 -k2,2n > "$genes_bed"

annotate_unique_clusters() {
  local unique_bed="$1"
  local prefix="$2"
  local stranded="${3:-no}"
  local intersect_file="${prefix}_intersected_genes.bed"
  local collapsed_file="${prefix}_genes_collapsed.tsv"
  local output_file="${prefix}_clusters_withGenes.tsv"

  if [[ "$stranded" == "yes" ]]; then
    bedtools intersect -a "$unique_bed" -b "$genes_bed" -s -wa -wb > "$intersect_file"
  else
    bedtools intersect -a "$unique_bed" -b "$genes_bed" -wa -wb > "$intersect_file"
  fi

  # Set intersect_gene_name_column to the gene-name field in the bedtools -wa -wb output.
  awk -v gene_col="$intersect_gene_name_column" 'BEGIN{FS=OFS="\t"}
  {
    key=$1":"$2"-"$3":"$6
    gene=$gene_col
    if (gene=="") gene="NA"
    if (!seen[key,gene]++) {
      genes[key] = (genes[key]=="" ? gene : genes[key]","gene)
    }
  }
  END{
    for (k in genes) print k, genes[k]
  }' "$intersect_file" > "$collapsed_file"

  # ARGIND keeps clusters with no gene hits as NA, even if collapsed_file is empty.
  awk 'BEGIN{FS=OFS="\t";
    print "Cluster_ID","Chromosome","Start","End","Genes","Score","Strand"
  }
  ARGIND==1 {
    g[$1]=$2
    next
  }
  ARGIND==2 {
    key=$1":"$2"-"$3":"$6
    cid="cluster_" FNR
    gene_list=(key in g ? g[key] : "NA")
    print cid,$1,$2,$3,gene_list,$5,$6
  }' "$collapsed_file" "$unique_bed" > "$output_file"
}

# If proTRAC cluster strand is ".", use "no"; if strand is reliable, use "yes".
annotate_unique_clusters "${treatment_label}_unique_noOverlap.bed" "${treatment_label}" "no"
annotate_unique_clusters "${control_label}_unique_noOverlap.bed" "${control_label}" "no"

echo "Done."
echo "Treatment unique clusters: ${treatment_label}_unique_noOverlap.bed"
echo "Control unique clusters: ${control_label}_unique_noOverlap.bed"
echo "Treatment clusters with genes: ${treatment_label}_clusters_withGenes.tsv"
echo "Control clusters with genes: ${control_label}_clusters_withGenes.tsv"
