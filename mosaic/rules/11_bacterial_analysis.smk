
def input_microbial_merge(wildcards):
	input_list=[]
	input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.tot.fasta", sample=SAMPLES)),
	if len(config['additional_reference_contigs'])>0:
		input_list.append(config['additional_reference_contigs'])
	return input_list
rule merge_microbial:
	input:
		assembled_contigs=input_microbial_merge,
	output:
		combined_positive_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial.tot.fasta",
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
		derreplicated_microbial_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		derreplicated_tmp=directory(dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp"),
	params:
		rep_name="combined_microbial_derreplicated_tot_tmp",
		rep_name_full=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot_tmp_rep_seq.fasta",
		rep_temp="combined_microbial_derreplicated_tot_tmp",
		dir_assembly=dirs_dict["ASSEMBLY_DIR"],
		min_len=1000
	message:
		"Derreplicating assembled contigs with mmseqs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/merge_microbial/tot.tsv"
	threads: 16
	shell:
		"""
		cat {input.assembled_contigs} | seqtk seq -l 0 - > {output.combined_positive_contigs}
		seqtk seq -L 2000 {output.combined_positive_contigs} > {output.combined_positive_contigs_2k}
		cd {params.dir_assembly}
		mmseqs easy-cluster --createdb-mode 1 --min-seq-id 1 -c 1 --cov-mode 1 {output.combined_positive_contigs} {params.rep_name} {params.rep_temp} --threads {threads}
		seqtk seq -L 2000 {params.rep_name_full} > {output.derreplicated_microbial_contigs}
		rm -f -- {params.rep_name:q}_all_seqs.fasta
		"""
		

rule buildBowtieDB_microbial:
	input:
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
	output:
		contigs_bt2=temp(expand(dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"])),
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot",
	message:
		"Creating contig DB with Bowtie2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBowtieDB_microbial/tot.tsv"
	threads: 64
	shell:
		"""
		bowtie2-build {input.combined_positive_contigs_2k} {params.prefix} --threads {threads}
		"""

rule mapReadsToContigs_microbial:
	input:
		contigs_bt2=expand(dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot.{part}.bt2", part=["1", "2", "3", "4", "rev.1", "rev.2"]),
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.{sampling}.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.{sampling}.fastq.gz"),
	output:
		sam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}.sam"),
		bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}.bam"),
		sorted_bam=(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_sorted.bam"),
		sorted_bam_idx=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_sorted.bam.bai"),
		filtered_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_filtered.bam"),
		flagstats=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_{sample}_.{sampling}.txt",
		flagstats_filtered=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_filtered_{sample}.{sampling}.txt",
		flagstats_unique=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_flagstats_unique_{sample}.{sampling}.txt",
		unique_sam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique.sam"),
		unique_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique.bam"),
		unique_sorted_bam=temp(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_sorted.bam"),
		covstats=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_covstats.txt",
		covstats_unique=dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_covstats.txt",
		basecov=(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_basecov.txt.gz") if PER_BASE_COVERAGE else [],
		unique_basecov=(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_{sampling}_unique_basecov.txt.gz") if PER_BASE_COVERAGE else [],
	params:
		prefix=dirs_dict["ASSEMBLY_DIR"] + "/2K_combined_microbial.tot",
	message:
		"Mapping microbial reads to assembly"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_mapping.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/mapReadsToContigs_microbial/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bowtie2 -x {params.prefix} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam} --threads {threads} --no-unal --all
		samtools view  -@ {threads} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {threads} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		coverm filter -b {output.sorted_bam} -o {output.filtered_bam} --min-read-percent-identity 95 --min-read-aligned-percent 85 -t {threads}
		samtools flagstat {output.filtered_bam} > {output.flagstats_filtered}
		samtools view -@ 144 -hf 0x2 {output.filtered_bam} | grep -v "XS:i:" > {output.unique_sam}
		samtools view  -@ 144 -bS {output.unique_sam}> {output.unique_bam}
		samtools sort -@ 144 {output.unique_bam} -o {output.unique_sorted_bam}
		samtools index {output.unique_sorted_bam}
		samtools flagstat {output.unique_bam}> {output.flagstats_unique}
		#genomecov
		for basecov in {output.basecov:q}; do
			bedtools genomecov -dz -ibam {output.filtered_bam} | gzip -c > "$basecov"
		done
		for basecov in {output.unique_basecov:q}; do
			bedtools genomecov -dz -ibam {output.unique_sorted_bam} | gzip -c > "$basecov"
		done
		#covstats
		coverm contig -b {output.filtered_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats}
		coverm contig -b {output.unique_sorted_bam} -m mean length covered_bases count variance trimmed_mean rpkm  -o {output.covstats_unique}
		"""

rule bacterial_binning_VAMB:
	input:
		combined_positive_contigs_2k=dirs_dict["ASSEMBLY_DIR"]+ "/2K_combined_microbial.tot.fasta",
		sorted_bam=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam", sample=SAMPLES_NO_TECHNICAL),
		sorted_bam_index=expand(dirs_dict["MAPPING_DIR"]+ "/MICROBIAL/bowtie2_{sample}_tot_sorted.bam.bai", sample=SAMPLES_NO_TECHNICAL),
	output:
		vamb_bins=directory(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins"),
	params:
		vamb_outdir_temp=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results_temp"),
		vamb_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results"),
		min_votu_len=config['min_votu_length'],
	message:
		"Binning microbial contigs with vamb"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/bacterial_binning_VAMB/tot.tsv"
	threads: 64
	shell:
		"""
		rm -rf {params.vamb_outdir_temp}
		vamb -o "_" --outdir {params.vamb_outdir_temp} --fasta {input.combined_positive_contigs_2k}  \
				--bamfiles {input.sorted_bam} --minfasta {params.min_votu_len} -p {threads}
		mkdir {params.vamb_outdir}
		mv {params.vamb_outdir_temp}/* {params.vamb_outdir}
		rm -rf {params.vamb_outdir_temp}
		for d in {params.vamb_outdir}/bins/*/ ; do cp ${{d}}*fna {output.vamb_bins} & done
		"""

rule predict_spacers:
	input:
		combined_positive_contigs=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_tot.fasta",
		minced_dir=(os.path.join(workflow.basedir, config['minced_dir'])),
	output:
		spacers=(dirs_dict["ANNOTATION"] + "/minced_predicted_spacers.tsv"),
	message:
		"Getting CRISPR spacers with MinCED"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/predict_spacers/tot.tsv"
	threads: 64
	shell:
		"""
		"$CONDA_PREFIX/bin/java" -jar {input.minced_dir:q}/minced.jar -spacers {input.combined_positive_contigs:q} {output.spacers:q}
		"""

rule estimateBinningQuality:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_temp")),
		checkMoutdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM"),
		# checkMoutplots=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_plots"),
	params:
		checkm_table=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/tab_results_checkM.csv"),
		checkm_outfile=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM/output_results_checkM.txt"),
		all_bins=dirs_dict["ASSEMBLY_DIR"] + "/microbial_checkM_temp/all_bins",
	log:
		checkMoutdir=(dirs_dict["vOUT_DIR"] + "/microbial_checkM_log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateBinningQuality/tot.tsv"
	threads: 32
	shell:
		"""
		mkdir -p {output.checkMoutdir_temp}
		cp -r {input.vamb_bins} {output.checkMoutdir_temp}
		cd {params.all_bins}
		checkm lineage_wf --tab_table -t {threads} -f {params.checkm_outfile} -x fna {params.all_bins} {output.checkMoutdir} 1> {log}
		"""

rule taxonomy_binning:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		gtdbtk_db=(config['gtdbtk_db']),
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/microbial_vamb_GTDB-Tk"),
		tempdir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_dir")),
	params:
		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_vamb_GTDB-Tk_mash"),
	message:
		"Assigning microbial taxonomy with GTDB-Tk "
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/taxonomy_binning/tot.tsv"
	threads: 64
	shell:
		"""
		mkdir {output.tempdir}
		cd {output.tempdir}
		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
		gtdbtk classify_wf --genome_dir {input.vamb_bins}/ --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fna
		"""

rule taxonomy_binning_assembly:
	input:
		racoon_assembly=expand(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample_nanopore}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta", sample_nanopore=NANOPORE_SAMPLES, sampling=SAMPLING_TYPE),
		gtdbtk_db=(config['gtdbtk_db']),
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_racon_microbial_GTDB-Tk"),
	params:
		mash_outdir=(dirs_dict["ASSEMBLY_DIR"] + "/microbial_racon_GTDB-Tk_mash"),
	message:
		"Assigning microbial taxonomy with GTDB-Tk "
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/taxonomy_binning_assembly/tot.tsv"
	threads: 64
	shell:
		"""
		mkdir {output.GTDB_outdir}/input_assemblies
		cp {input.racoon_assembly} {output.GTDB_outdir}/input_assemblies
		conda env config vars set GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
		gtdbtk classify_wf --genome_dir {output.GTDB_outdir}/input_assemblies/ --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fasta
		"""
rule DRAM_microbial_annotation:
	input:
		vamb_bins=(dirs_dict["ASSEMBLY_DIR"] + "/vamb_binning_results/all_bins/"),
		DRAM_db=config['DRAM_db'],
	output:
		DRAM_output=directory(dirs_dict["ANNOTATION"]+ "/DRAM_annotate_results_{sampling}"),
		DRAM_summary=directory(dirs_dict["ANNOTATION"]+ "/DRAM_distill_results_{sampling}"),
	params:
		DRAM_annotations=dirs_dict["ANNOTATION"]+ "/DRAM_annotate_results_{sampling}/annotations.tsv",
		# trna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/trnas.tsv"),
		# rrna=directory(dirs_dict["vOUT_DIR"]+ "/DRAM_combined_" + VIRAL_CONTIGS_BASE + "_derreplicated_rep_seq_{sampling}/rrnas.tsv"),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/DRAM_microbial_annotation/sampling={sampling}.tsv"
	message:
		"Annotate contigs with DRAM"
	threads: 32
	shell:
		"""
		DRAM.py annotate -i '{input.DAS_Tool_bins}/*fna' -o {output.DRAM_output} --threads {threads}
		DRAM.py distill -i {params.DRAM_annotations} -o {output.DRAM_summary} 
		"""

def input_taxonomy_gtdbtk_bacteria_all(wildcards):
	input_list=[]
	if NANOPORE & (NANOPORE_ONLY):
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER +"_corrected_scaffolds_pilon.{sampling}.fasta", sample=NANOPORE_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & (PACBIO_ONLY):
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	if PACBIO & (PACBIO_HYBRID):
		input_list.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta", sample=PACBIO_SAMPLES, sampling=wildcards.sampling))
	return(input_list)


rule taxonomy_gtdbtk_bacteria:
	input:
		assemblies=input_taxonomy_gtdbtk_bacteria_all,
		gtdbtk_db=config['gtdbtk_db']
	output:
		GTDB_outdir=directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_bacteria_GTDB-Tk_{sampling}"),
		GTDB_temp=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/assembly_bacteria_singlefasta_GTDB-Tk_{sampling}"))
	params:
		mash_outdir=dirs_dict["ASSEMBLY_DIR"] + "/assembly_bacteria_GTDB-Tk_mash_{sampling}"
	message:
		"Assigning bacterial taxonomy with GTDB-Tk"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/taxonomy_gtdbtk_bacteria/sampling={sampling}.tsv"
	threads: 64
	shell:
		"""
		mkdir -p {output.GTDB_temp}
		for assembly in {input.assemblies}; do sample=$(basename "$assembly" .fasta); cp "$assembly" {output.GTDB_temp}/"$sample".fasta; done
		export GTDBTK_DATA_PATH={input.gtdbtk_db}/release214/
		gtdbtk classify_wf --genome_dir {output.GTDB_temp} --out_dir {output.GTDB_outdir} --cpus {threads} --mash_db {params.mash_outdir} --extension fasta
		"""

rule single_fasta_microbial:
	input:
		derreplicated_microbial_contigs=lambda wildcards: SOURMASH_CONTIG_FASTAS[wildcards.sourmash_catalogue],
	output:
		derreplicated_microbial_contigs_dir=temp(directory(dirs_dict["ASSEMBLY_DIR"] + "/single_{sourmash_catalogue}")),
	wildcard_constraints:
		sourmash_catalogue="|".join(SOURMASH_CONTIG_FASTAS),
	message:
		"Formatting contig representatives into single FASTA files"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/single_fasta_microbial/catalogue={sourmash_catalogue}.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.derreplicated_microbial_contigs_dir:q}
		if [ -s {input.derreplicated_microbial_contigs:q} ]; then
			seqkit split --quiet -w 0 -i {input.derreplicated_microbial_contigs:q} --out-dir {output.derreplicated_microbial_contigs_dir:q}
		fi
	 	"""

rule sourmash_sketch_microbial:
	input:
		derreplicated_microbial_contigs=lambda wildcards: SOURMASH_CONTIG_FASTAS[wildcards.sourmash_catalogue],
		derreplicated_microbial_contigs_dir=dirs_dict["ASSEMBLY_DIR"] + "/single_{sourmash_catalogue}",
	output:
		manysketch_csv=temp(dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_manysketch.csv"),
		sketch=temp(dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_sourmash.sig.zip"),
		query_stats=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_sourmash_queries.tsv",
	wildcard_constraints:
		sourmash_catalogue="|".join(SOURMASH_CONTIG_FASTAS),
	params: 
		min_length=lambda wildcards: 0 if wildcards.sourmash_catalogue == "combined_microbial_derreplicated_tot" else int(config.get("sourmash_contig_catalogue_min_length", 5000)),
	message:
		"Building sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_microbial/catalogue={sourmash_catalogue}.tsv"
	threads: 64
	shell:
		"""
		python - {input.derreplicated_microbial_contigs:q} {input.derreplicated_microbial_contigs_dir:q} {output.manysketch_csv:q} {params.min_length} <<-'PYTHON'
		import csv
		from pathlib import Path
		import sys
		import screed
		fasta, split_dir, csv_path, minimum = sys.argv[1:]
		with open(csv_path, "w") as handle:
		    writer = csv.writer(handle)
		    writer.writerow(["name", "genome_filename", "protein_filename"])
		    if Path(fasta).stat().st_size:
		        for record in screed.open(fasta):
		            if len(record.sequence) >= int(minimum):
		                name = record.name.split()[0]
		                path = Path(split_dir) / (Path(fasta).stem + ".part_" + name + ".fasta")
		                writer.writerow([name, str(path), ""])
		PYTHON
		if [ "$(wc -l < {output.manysketch_csv:q})" -gt 1 ]; then
			sourmash scripts manysketch {output.manysketch_csv:q} -p k=31,abund,DNA -o {output.sketch:q} -c {threads}
		else
			python - {output.sketch:q} <<-'PYTHON'
			import sys
			from sourmash.sourmash_args import SaveSignaturesToLocation
			with SaveSignaturesToLocation(sys.argv[1]):
			    pass
			PYTHON
		fi
		python - {input.derreplicated_microbial_contigs:q} {output.sketch:q} {output.query_stats:q} {params.min_length} <<-'PYTHON'
		import csv
		from pathlib import Path
		import sys
		import screed
		import sourmash
		fasta, sketches, output_path, minimum = sys.argv[1:]
		hashes = {{signature.name: len(signature.minhash) for signature in sourmash.load_file_as_signatures(sketches)}}
		with open(output_path, "w") as handle:
		    writer = csv.writer(handle, delimiter="\t")
		    writer.writerow(["cluster_rep", "sourmash_representative_length_bp", "sourmash_query_hashes", "sourmash_screen_status"])
		    if Path(fasta).stat().st_size:
		        for record in screed.open(fasta):
		            name = record.name.split()[0]
		            length = len(record.sequence)
		            count = hashes.get(name, 0)
		            status = "below_length_cutoff" if length < int(minimum) else ("no_usable_hashes" if count == 0 else "screened")
		            writer.writerow([name, length, count, status])
		PYTHON
		"""

rule sourmash_gather_microbial:
	input:
		sketch=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_sourmash.sig.zip",
		query_stats=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_sourmash_queries.tsv",
		sourmash_rocksdb=config['sourmash_rocksdb'],
	output:
		gather=temp(dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_gather_sourmash.csv"),
	wildcard_constraints:
		sourmash_catalogue="|".join(SOURMASH_CONTIG_FASTAS),
	message:
		"Genome containtment with sourmash gather"
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_gather_microbial/catalogue={sourmash_catalogue}.tsv"
	threads: 64
	shell:
		"""
		printf '' > {output.gather:q}
		if awk -F '\t' 'NR > 1 && $3 > 0 {{found=1}} END {{exit !found}}' {input.query_stats:q}; then
			sourmash scripts fastmultigather {input.sketch:q} {input.sourmash_rocksdb:q} -c {threads} -o {output.gather:q} -t {params.threshold_bp} -s 1000
		fi
		if [ ! -s {output.gather:q} ]; then
			printf 'query_name,name,f_unique_weighted,f_unique_to_query,unique_intersect_bp,remaining_bp,query_md5,query_filename,query_bp,ksize,scaled,query_n_hashes\n' > {output.gather:q}
		fi
		"""

rule sourmash_tax_microbial:
	input:
		gather=dirs_dict["ANNOTATION"] + "/{sourmash_catalogue}_gather_sourmash.csv",
		sourmash_tax=config['sourmash_tax'],
	output:
		csv_report=dirs_dict["ANNOTATION"] + "/sourmash_{sourmash_catalogue}.classifications.csv",
	wildcard_constraints:
		sourmash_catalogue="|".join(SOURMASH_CONTIG_FASTAS),
	params:
		outdir=(dirs_dict["ANNOTATION"]),
		name="sourmash_{sourmash_catalogue}",
	message:
		"Assigning taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax_microbial/catalogue={sourmash_catalogue}.tsv"
	threads: 1
	shell:
		"""
		if [ "$(wc -l < {input.gather:q})" -gt 1 ]; then
			sourmash tax genome --gather-csv {input.gather:q} -t {input.sourmash_tax:q} -o {params.name:q}\
				--output-dir {params.outdir:q} -F csv_summary
		else
			printf 'query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank\n' > {output.csv_report:q}
		fi
		"""
		
rule defense_finder:
	input:
		aa=dirs_dict["ASSEMBLY_DIR"]+ "/combined_microbial_derreplicated_ORFs_tot.faa",
	output:
		defenseFinder_dir=directory(dirs_dict["ANNOTATION"]+ "/DefenseFinder_results_{sampling}/"),
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/defense_finder/sampling={sampling}.tsv"
	message:
		"Detecting defense systems with DefenseFinder"
	threads: 32
	shell:
		"""
		python scripts/process_fasta_satellite_finder.py {params.faa} {output.faa_temp}
		defense-finder run -–db-type gembase -w {threads} --out-dir {output.defenseFinder_dir} {input.aa}
		"""

def input_estimateBacterialGenomeCompletness(wildcards):
	input_list=[]
	if NANOPORE & (NANOPORE_ONLY):
		return(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta")
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER +"_corrected_scaffolds_pilon.{sampling}.fasta")
	if PACBIO & (PACBIO_ONLY):
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
	if PACBIO & (PACBIO_HYBRID):
		return(dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
	if ISOLATES:
		return(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.{sampling}.fasta")


rule estimateBacterialGenomeCompletness:
	input:
		fasta=input_estimateBacterialGenomeCompletness,
		checkm_db=(config['checkm_db']),
	output:
		checkMoutdir_temp=temp(directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}_temp")),
		checkMoutdir=(directory(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}")),
	params:
		checkv_db=dirs_dict["vOUT_DIR"] + "/{sample}_checkV_{sampling}",
	log:
		checkMoutdir=temp(dirs_dict["vOUT_DIR"] + "/{sample}_checkM_{sampling}.log"),
	message:
		"Estimating genome completeness with CheckM "
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/estimateBacterialGenomeCompletness/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		mkdir -p {output.checkMoutdir_temp}
		cp {input.fasta} {output.checkMoutdir_temp}
		cd {output.checkMoutdir_temp}
		checkm lineage_wf -t {threads} -x fasta {output.checkMoutdir_temp} {output.checkMoutdir} 1> {log}
		"""

rule combine_logs_to_csv:
	input:
		logs = expand((dirs_dict["vOUT_DIR"] + "/{sample}_checkM_tot.log"), sample=SAMPLES)
	output:
		csv = dirs_dict["PLOTS_DIR"] + "/checkM_summary.csv"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/combine_logs_to_csv/tot.tsv"
	run:
		import pandas as pd
		from io import StringIO
		import glob
		dfs = []
		for log_file in input.logs:
			with open(log_file, "r") as f:
					lines = f.readlines()

			# Find the table start and end
			start_idx = None
			end_idx = None
			for i, line in enumerate(lines):
					if line.strip().startswith("Bin Id"):
						start_idx = i
					elif start_idx and line.strip().startswith("[") and "INFO" in line:
						end_idx = i
						break

			if start_idx is not None:
					# Extract the table lines
					table_lines = lines[start_idx:end_idx]
					# Remove separator lines (-----)
					table_lines = [l for l in table_lines if not l.strip().startswith('---')]
					# Join lines and read with pandas
					table_str = "".join(table_lines)
					df = pd.read_csv(StringIO(table_str), sep=r'\s{2,}', engine='python')
					dfs.append(df)

		# Combine all tables
		final_df = pd.concat(dfs, ignore_index=True)
		final_df.to_csv(output.csv, index=False)

rule fastani_all_vs_all:
	input:
		fasta=expand(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.tot.fasta", sample=SAMPLES),
	output:
		fastANI=dirs_dict["ASSEMBLY_DIR"]+ "/spades_filtered_scaffolds_ANI.tot.csv",
		query_list=temp(dirs_dict["ASSEMBLY_DIR"]+ "/spades_filtered_scaffolds_ANI_query.tot.csv"),
	threads: 8
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/fastani_all_vs_all/tot.tsv"
	shell:
		"""
		mkdir -p fastani
		printf "%s\n" {input.fasta} > {output.query_list}
		fastANI --matrix -t {threads} \
			--ql {output.query_list} \
			--rl {output.query_list} \
			-o {output.fastANI} 
			
		"""

rule single_fasta_microbial_isolate:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
	output:
		single_contigs_dir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_single_{sampling}")),
	message:
		"formating microbial contigs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/single_fasta_microbial_isolate/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		seqkit split --quiet -w 0 -i {input.fasta} --out-dir {output.single_contigs_dir}
	 	"""

rule sourmash_sketch_microbial_isolate:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_spades_filtered_scaffolds.{sampling}.fasta",
		single_contigs_dir=((dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_single_{sampling}")),
	output:
		manysketch_csv=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_manysketch.csv"),
		sketch=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_sourmash.sig.zip"),
	params: 
		name="{sample}_spades_filtered_scaffolds.{sampling}"
	message:
		"Building sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_microbial_isolate/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		echo name,genome_filename,protein_filename > {output.manysketch_csv}
		grep "^>" {input.fasta} | sed 's/^>//' | awk -v dir="{input.single_contigs_dir}/{params.name}.part_" '{{print $1 "," dir $1 ".fasta,"}}' >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,abund,DNA -o {output.sketch} -c {threads}
		"""

rule sourmash_gather_microbial_isolate:
	input:
		sketch=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_sourmash.sig.zip"),
		sourmash_rocksdb=config['sourmash_rocksdb'],
	output:
		gather=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_gather_sourmash.csv"),
	message:
		"Genome containtment with sourmash gather"
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_gather_microbial_isolate/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
		"""


rule sourmash_tax_microbial_isolate:
	input:
		gather=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_gather_sourmash.csv"),
		sourmash_tax=config['sourmash_tax'],
	output:
		csv_report=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}.classifications.csv"),
	params:
		outdir=(dirs_dict["ASSEMBLY_DIR"]),
		name="{sample}_{sampling}"
	message:
		"Assigning taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax_microbial_isolate/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax}  -o {params.name}\
			--output-dir {params.outdir} -F csv_summary --rank strain
		"""


rule single_fasta_pacbio:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta",
	output:
		single_contigs_dir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_pacbio_single_{sampling}")),
	message:
		"formating PacBio contigs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/single_fasta_pacbio/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		seqkit split --quiet -w 0 -i {input.fasta} --out-dir {output.single_contigs_dir}
	 	"""

rule sourmash_sketch_pacbio:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta",
		single_contigs_dir=((dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_pacbio_single_{sampling}")),
	output:
		manysketch_csv=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_manysketch.csv"),
		sketch=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_sourmash.sig.zip"),
	params: 
		name="{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}"
	message:
		"Building PacBio sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_pacbio/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		echo name,genome_filename,protein_filename > {output.manysketch_csv}
		grep "^>" {input.fasta} | sed 's/^>//' | awk -v dir="{input.single_contigs_dir}/{params.name}.part_" '{{print $1 "," dir $1 ".fasta,"}}' >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,abund,DNA -o {output.sketch} -c {threads}
		"""

rule sourmash_gather_pacbio:
	input:
		sketch=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_sourmash.sig.zip"),
		sourmash_rocksdb=config['sourmash_rocksdb'],
	output:
		gather=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_gather_sourmash.csv"),
	message:
		"Genome containment with sourmash gather"
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_gather_pacbio/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
		"""

rule sourmash_tax_pacbio:
	input:
		gather=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_gather_sourmash.csv"),
		sourmash_tax=config['sourmash_tax'],
	output:
		csv_report=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio.classifications.csv"),
	params:
		outdir=(dirs_dict["ASSEMBLY_DIR"]),
		name="{sample}_{sampling}_pacbio"
	message:
		"Assigning PacBio taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax_pacbio/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax}  -o {params.name}\
			--output-dir {params.outdir} -F csv_summary --rank strain
		"""

rule single_fasta_pacbio_hybrid:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta",
	output:
		single_contigs_dir=temp(directory(dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_pacbio_hybrid_single_{sampling}")),
	message:
		"formating PacBio hybrid contigs into single fasta"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/single_fasta_pacbio_hybrid/sample={sample}__sampling={sampling}.tsv"
	threads: 1
	shell:
		"""
		seqkit split --quiet -w 0 -i {input.fasta} --out-dir {output.single_contigs_dir}
	 	"""

rule sourmash_sketch_pacbio_hybrid:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"]+ "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta",
		single_contigs_dir=((dirs_dict["ASSEMBLY_DIR"]+ "/{sample}_pacbio_hybrid_single_{sampling}")),
	output:
		manysketch_csv=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid_manysketch.csv"),
		sketch=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid_sourmash.sig.zip"),
	params: 
		name="polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}"
	message:
		"Building PacBio hybrid sketches with sourmash"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_pacbio_hybrid/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		echo name,genome_filename,protein_filename > {output.manysketch_csv}
		grep "^>" {input.fasta} | sed 's/^>//' | awk -v dir="{input.single_contigs_dir}/{params.name}.part_" '{{print $1 "," dir $1 ".fasta,"}}' >> {output.manysketch_csv}
		sourmash scripts manysketch {output.manysketch_csv} -p k=31,abund,DNA -o {output.sketch} -c {threads}
		"""

rule sourmash_gather_pacbio_hybrid:
	input:
		sketch=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid_sourmash.sig.zip"),
		sourmash_rocksdb=config['sourmash_rocksdb'],
	output:
		gather=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid_gather_sourmash.csv"),
	message:
		"Genome containment with sourmash gather"
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_gather_pacbio_hybrid/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
		"""

rule sourmash_tax_pacbio_hybrid:
	input:
		gather=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid_gather_sourmash.csv"),
		sourmash_tax=config['sourmash_tax'],
	output:
		csv_report=(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_pacbio_hybrid.classifications.csv"),
	params:
		outdir=(dirs_dict["ASSEMBLY_DIR"]),
		name="{sample}_{sampling}_pacbio_hybrid"
	message:
		"Assigning PacBio hybrid taxonomy with sourmash tax"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax_pacbio_hybrid/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax}  -o {params.name}\
			--output-dir {params.outdir} -F csv_summary --rank strain
		"""

rule sourmash_sketch_nanopore_only_bacteria:
	input:
		fasta=dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta"
	output:
		sketch=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_sourmash.sig.zip")
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_sketch_nanopore_only_bacteria/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		sourmash sketch dna -p k=31,abund {input.fasta} -o {output.sketch}
		"""

rule sourmash_gather_nanopore_only_bacteria:
	input:
		sketch=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_sourmash.sig.zip",
		sourmash_rocksdb=config['sourmash_rocksdb']
	output:
		gather=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_gather_sourmash.csv")
	params:
		threshold_bp="0"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_gather_nanopore_only_bacteria/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
		"""
		
rule sourmash_tax_nanopore_only_bacteria:
	input:
		gather=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_gather_sourmash.csv",
		sourmash_tax=config['sourmash_tax']
	output:
		csv_report=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore.classifications.csv"
	params:
		outdir=dirs_dict["ASSEMBLY_DIR"],
		name="{sample}_{sampling}_nanopore"
	conda:
		dirs_dict["ENVS_DIR"]+ "/sourmash.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/sourmash_tax_nanopore_only_bacteria/sample={sample}__sampling={sampling}.tsv"
	threads: 4
	shell:
		"""
		sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax} -o {params.name} --output-dir {params.outdir} -F csv_summary --rank strain
		"""
def input_bakta_assembly(wildcards):
	if NANOPORE & NANOPORE_ONLY:
		return(dirs_dict["ASSEMBLY_DIR"] + "/racon_{sample}_contigs_2_"+ LONG_ASSEMBLER + ".{sampling}.fasta")
	if NANOPORE & (not NANOPORE_ONLY) & PAIRED:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_"+ LONG_ASSEMBLER +"_corrected_scaffolds_pilon.{sampling}.fasta")
	if PACBIO & PACBIO_ONLY:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
	if PACBIO & PACBIO_HYBRID:
		return(dirs_dict["ASSEMBLY_DIR"] + "/polypolish_{sample}_contigs_"+ LONG_ASSEMBLER_PACBIO + ".{sampling}.fasta")
	if ISOLATES:
		return(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_spades_filtered_scaffolds.{sampling}.fasta")

rule annotate_bakta:
	input:
		assembly=input_bakta_assembly,
		db=(config["bakta_db"]),
	output:
		outdir=directory(dirs_dict["ANNOTATION"] + "/bakta_{sample}_{sampling}"),
		gff=dirs_dict["ANNOTATION"] + "/bakta_{sample}_{sampling}/{sample}.gff3",
		ffn=dirs_dict["ANNOTATION"] + "/bakta_{sample}_{sampling}/{sample}.ffn",
		faa=dirs_dict["ANNOTATION"] + "/bakta_{sample}_{sampling}/{sample}.faa",
		tsv=dirs_dict["ANNOTATION"] + "/bakta_{sample}_{sampling}/{sample}.tsv"
	params:
		prefix="{sample}",
	message:
		"Annotating bacterial/fungal assembly with Bakta"
	conda:
		dirs_dict["ENVS_DIR"] + "/bakta.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/annotate_bakta/sample={sample}__sampling={sampling}.tsv"
	threads: 8
	shell:
		"""
		bakta --db {input.db} --threads {threads} --prefix {params.prefix} --output {output.outdir} --force {input.assembly}
		"""

rule sourmash_sketch_nanopore_hybrid_bacteria:
    benchmark:
        dirs_dict["BENCHMARKS"] + "/sourmash_sketch_nanopore_hybrid_bacteria/sample={sample}__sampling={sampling}.tsv"
    input:
        fasta=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_" + LONG_ASSEMBLER + "_corrected_scaffolds_pilon.{sampling}.fasta"
    output:
        sketch=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid_sourmash.sig.zip")
    conda:
        dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
    threads: 4
    shell:
        """
        sourmash sketch dna -p k=31,abund {input.fasta} -o {output.sketch}
        """

rule sourmash_gather_nanopore_hybrid_bacteria:
    benchmark:
        dirs_dict["BENCHMARKS"] + "/sourmash_gather_nanopore_hybrid_bacteria/sample={sample}__sampling={sampling}.tsv"
    input:
        sketch=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid_sourmash.sig.zip",
        sourmash_rocksdb=config['sourmash_rocksdb']
    output:
        gather=temp(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid_gather_sourmash.csv")
    params:
        threshold_bp="0"
    conda:
        dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
    threads: 8
    shell:
        """
        sourmash scripts fastmultigather {input.sketch} {input.sourmash_rocksdb} -c {threads} -o {output.gather} -t {params.threshold_bp} -s 1000
        """

rule sourmash_tax_nanopore_hybrid_bacteria:
    benchmark:
        dirs_dict["BENCHMARKS"] + "/sourmash_tax_nanopore_hybrid_bacteria/sample={sample}__sampling={sampling}.tsv"
    input:
        gather=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid_gather_sourmash.csv",
        sourmash_tax=config['sourmash_tax']
    output:
        csv_report=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_{sampling}_nanopore_hybrid.classifications.csv"
    params:
        outdir=dirs_dict["ASSEMBLY_DIR"],
        name="{sample}_{sampling}_nanopore_hybrid"
    conda:
        dirs_dict["ENVS_DIR"] + "/sourmash.yaml"
    threads: 4
    shell:
        """
        sourmash tax genome --gather-csv {input.gather} -t {input.sourmash_tax} -o {params.name} --output-dir {params.outdir} -F csv_summary --rank strain
        """
