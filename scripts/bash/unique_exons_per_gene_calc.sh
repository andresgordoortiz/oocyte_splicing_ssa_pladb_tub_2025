awk '
  BEGIN { FS = "\t"; OFS = "\t" }
  $3 == "exon" &&
  ($9 ~ /gene_biotype "protein_coding"/ || $9 ~ /gene_type "protein_coding"/) {
    # extract gene_id
    match($9, /gene_id "([^"]+)"/, m)
    gid = m[1]

    # unique span: chr:start-end
    span = $1 ":" $4 "-" $5

    # composite key: gene_id SUBSEP span
    key = gid SUBSEP span

    if (!seen[key]++) {
      count[gid]++
    }
  }
  END {
    for (g in count) {
      print g, count[g]
    }
  }
' Mus_musculus.GRCm39.113.gtf \
| sort -k2,2nr -k1,1 \
> coding_unique_exon_counts_per_gene.txt
