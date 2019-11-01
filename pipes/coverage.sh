#Sort bam file
samtools sort {input.bam} -o {output.bam_sorted}
#Filter bam
bamm filter --bamfile {output.bam_sorted} --percentage_id 0.95 --percentage_aln 0.9 -o {params.out_dir}
#Compute tpmean headers(contig, length, tpmean)
bamm parse -c {output.tpmean} -m tpmean -b {output.bam_filtered}
#Compute genomecov 
bedtools genomecov -dz -ibam {input.bam_filtered} > {output.bam_cov}
#Parse genomecov headers(covered_length, contig)
cut -f 1 {output.bam_cov} | sort| uniq -c | sort -nr | sed -e 's/^[[:space:]]*//' > {output.cov_final}