demeter: /media/demeter/Seagate Expansion Drive/20200703_FS10000714_10_BPL20321-3127/Data/Intensities/BaseCalls/0/Wheat_virome/

cd /home/lmf/apps/Mosaic/workflow
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/00_RAW_READS min_len="0" min_cov="0" -k -j 16 -n min_len="1000"
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/newMosaic/00_RAW_READS -k -j 16 -n


#CONTAMINATION
cd /home/lmf/apps/Mosaic/workflow
snakemake --use-conda -p abundance_from_db_contigs --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment/00_RAW_READS representative_contigs=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment/Triticum_aestivum.tot.fasta results_dir=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment


#Mosaic + amplified reads

snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome_3/00_RAW_READS sampling="tot" -k -j 16 -n
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome_5/00_RAW_READS sampling="tot" -k -j 32 -n

#Nanopore
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome_6/00_RAW_READS nanopore_pooled="True" nanopore_pooled_name=rapid_virome sampling="tot" -j 64 -k

#APSE phages

APSE-1  APSE-1
APSE-2  APSE-2
APSE-3  APSE-3
APSE-4  APSE-2
APSE-5  APSE-2
APSE-6  APSE-6
APSE-7  APSE-7
APSE-8  APSE-2

APSE has 5 different types:   1 2 3 6 7
Full sequences for 4:         1 2 3   7

#Identify which APSE are representative
snakemake --use-conda -p abundance_from_db_contigs --config input_dir=/home/lmf/PhylloVir/Wheat_Virome_3/APSE/00_RAW_READS representative_contigs=/home/lmf/PhylloVir/Wheat_Virome_3/APSE/APSE_genomes.tot.fasta -k -j 32 -n
/home/lmf/apps/Mosaic/workflow/tools/weeSAM/weeSAM --bam BamM_*_sorted.tot.bam --html weeSAM_plots


contig_id	contig_length	gene_count	viral_genes	host_genes	checkv_quality	provirus	termini	seqname	group	name_x	rank	taxonomy_mmseqs	name_y	taxonomy_vcontact2	OTU	FL_depth	FL_amp_depth	OL_depth	OL_amp_depth
26	NODE_38_length_45369_cov_157.060004	45369	52	39	0	Complete	No	55-bp-DTR	NODE_38_length_45369_cov_157.060004	dsDNAphage	NODE_38_length_45369_cov_157.060004	genus	Zindervirus			NODE_38_length_45369_cov_157.060004	470.0358061	83.32669639	0	0
5288	NODE_34_length_44979_cov_93.040268	44979	54	32	0	Complete	No	55-bp-DTR	NODE_34_length_44979_cov_93.040268	dsDNAphage	NODE_34_length_44979_cov_93.040268	no rank	unclassified	NODE_34_length_44979_cov_93.040268	Autographiviridae [family]	NODE_34_length_44979_cov_93.040268	0	0	135.8739028	15.00752882
5274	NODE_7_length_90711_cov_145.333028	90711	188	55	2	Complete	No	55-bp-DTR	NODE_7_length_90711_cov_145.333028	dsDNAphage	NODE_7_length_90711_cov_145.333028	no rank	unclassified			NODE_7_length_90711_cov_145.333028	0	0	41.49901658	347.5742026
6874	NODE_4_length_146679_cov_145.852500	146679	263	140	1	Complete	No	55-bp-DTR	NODE_4_length_146679_cov_145.852500	dsDNAphage	NODE_4_length_146679_cov_145.852500	no rank	unclassified			NODE_4_length_146679_cov_145.852500	71.42930548	199.8586936	11.83895211	96.18311964

#Test depth
snakemake -p test_assembly_depth --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome_4/00_RAW_READS sampling="tot" max_norm=10000  -k -j 32 -n

#WGS
snakemake -p assembly --use-conda --config input_dir=/home/lmf/PhylloVir/WGS_2/00_RAW_READS sampling="tot" -k -j 32 -n
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/WGS_2/00_RAW_READS sampling="tot" -k -j 32 -n

#kraken2
kraken2 --db /home/lmf/db/KRAKEN/kraken/ --paired 02_CLEAN_DATA/OL_forward_paired_clean.tot.fastq 02_CLEAN_DATA/OL_reverse_paired_clean.tot.fastq --threads 32 --output OL_full_kraken_out.txt --use-names --report OL_full_kraken_report.txt &
kraken2 --db /home/lmf/db/KRAKEN/kraken/ --paired 02_CLEAN_DATA/FL_forward_paired_clean.tot.fastq 02_CLEAN_DATA/FL_reverse_paired_clean.tot.fastq --threads 32 --output FL_full_kraken_out.txt --use-names --report FL_full_kraken_report.txt

#What the phage
nextflow run replikation/What_the_Phage -r v1.0.0 --cores 8 -profile local,docker --fasta merged_scaffolds.tot.fasta

#CRISPR
conda activate Mosaic
cd ~/PhylloVir/WGS_2/CRISPR
/home/lmf/apps/Mosaic/workflow/tools/pilercr1.06/pilercr -noinfo -in merged_scaffolds.tot.fasta -out piler_cr.txt
minced -spacers merged_scaffolds.tot.fasta minced.txt

/home/lmf/apps/spacepharer/build/bin/spacepharer createsetdb viral_merged_scaffolds.tot.fasta viralTargetDB tmpFolder
/home/lmf/apps/spacepharer/build/bin/spacepharer createsetdb viral_merged_scaffolds.tot.fasta viralTargetDB_rev tmpFolder --reverse-fragments 1

/home/lmf/apps/spacepharer/build/bin/spacepharer parsespacer piler_cr.txt spacers_pilerDB

/home/lmf/apps/spacepharer/build/bin/spacepharer createsetdb minced_spacers.fa spacers_mincedSetDB tmpFolder --extractorf-spacer 1
/home/lmf/apps/spacepharer/build/bin/spacepharer createsetdb spacers_pilerDB spacers_pilerDBSetDB tmpFolder --extractorf-spacer 1
/home/lmf/apps/spacepharer/build/bin/spacepharer downloaddb spacers_shmakov_et_al_2017 shmakovSetDB tmpFolder

/home/lmf/apps/spacepharer/build/bin/spacepharer predictmatch spacers_mincedSetDB viralTargetDB viralTargetDB_rev spacepharer_minced.tsv tmpFolder
/home/lmf/apps/spacepharer/build/bin/spacepharer predictmatch spacers_pilerDBSetDB viralTargetDB viralTargetDB_rev spacepharer_piler.tsv tmpFolder
/home/lmf/apps/spacepharer/build/bin/spacepharer predictmatch shmakovSetDB viralTargetDB viralTargetDB_rev spacepharer_shmakov.tsv tmpFolder
rm *DB*
