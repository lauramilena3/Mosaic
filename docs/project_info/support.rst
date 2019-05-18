Cluster Execution
=================

The use of Snakemake allows Mosaic's jobs to be executed on cluster engines with SURFsara-Grid, PBS-TORQUE or Slurm job schedulers.

In order to execute Mosaic in a cluster, run the following code according to your job scheduler::

   snakemake --cluster qsub -j 32		#PBS-TORQUE
   snakemake --cluster sbatch -j 32		#Slurm

For additional information about cluster configuration, please refer to Snakemake examples 
avaliable `here <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration>`_.



Config File
===========

Many values can be modified depending on user needs. Mosaic’s JSON config file allows the user to select types of data and thresholds used for the analysis depending on their needs. This file may be used as supplementary material to replicate the analyses. Using this, plus the conda JSON files describing all dependencies versions, we expect the users to consistently retrieve the same results independent on the machine. 

In order to change this values, you can append the desired values to the main command. For example, in order to change the minimun length of a contig to 2000bp, you can use::

   snakemake -j $nCores --use-conda --config min_len=2000 basecalled_data=$fastqDir results_dir=$reusultDir 




Below, you can find all the variables, and their default values.  

Single End
----------

Mosaic automatically detects the Illumina input file. Please note that Spades does not support metagenomic assembly for SE files. Hence, Mosaic assembles SE reads using Spades single-cell flag.

General:
--------

+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| input_dir                  | "data/00_RAW_DATA"  | Input folder with raw data                  |
+----------------------------+---------------------+---------------------------------------------+
| results_dir                | "data"              | Directory for output files                  |
+----------------------------+---------------------+---------------------------------------------+
| operating_system           | "linux"             | Operating system ([linux], [macOs])         |
+----------------------------+---------------------+---------------------------------------------+
| forward_tag                | "R1"                | Forward tag used for Illumina reads         |
+----------------------------+---------------------+---------------------------------------------+
| reverse_tag                | "R2"                | Reverse tag used for Illumina reads         |
+----------------------------+---------------------+---------------------------------------------+
| nanopore_tag               | "nanopore"          | Tag used for Nanopore files                 |
+----------------------------+---------------------+---------------------------------------------+


Nanopore:
---------
+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| nanopore_pooled            | "False"             | Nanopore run pooling                        |
+----------------------------+---------------------+---------------------------------------------+
| nanopore_pooled_name       | "raw"               | Tag used for pooled Nanopore run            |
+----------------------------+---------------------+---------------------------------------------+



Trimmomatic:
------------

+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| trimmomatic_leading        | "15"                | Leading minimum base quality                |
+----------------------------+---------------------+---------------------------------------------+
| trimmomatic_trailing       | "15"                | Trailing minimum base quality               |
+----------------------------+---------------------+---------------------------------------------+
| trimmomatic_window_size    | "4"                 | Window size                                 |
+----------------------------+---------------------+---------------------------------------------+
| trimmomatic_window_quality | "15"                | Window minimum quality                      |
+----------------------------+---------------------+---------------------------------------------+
| trimmomatic_minlen         | "50"                | Read minimum length                         |
+----------------------------+---------------------+---------------------------------------------+

Cleaning:
---------

+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| contaminants_list          | "GCF_000819615"     | List of contaminants RefSeq accesion IDs    |
|                            |                     | Eg: "GCF_000819615 GCA_000929915.1"         |
+----------------------------+---------------------+---------------------------------------------+
| adapters_file              | "NexteraPE-PE.fa"   | Adapter file for barcode removal            |
|                            |                     | Please note that if you use other adapters  |
|                            |                     | you need to add them to pipes/db/adapters   |
+----------------------------+---------------------+---------------------------------------------+

Normalization and Subsampling:
------------------------------

+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| min_norm                   | "2"                 | Minimum normalization coverage              |
+----------------------------+---------------------+---------------------------------------------+
| max_norm                   | "100"               | Maximum normalization coverage              |
+----------------------------+---------------------+---------------------------------------------+
| max_subsample              | "20000000"          | Maximum number of reads for subsampling     |
+----------------------------+---------------------+---------------------------------------------+


Contigs:
--------

+----------------------------+---------------------+---------------------------------------------+
| Name                       | Default             | Description                                 |
+============================+=====================+=============================================+
| min_len                    | "1000"              | Minimum contig length                       |
+----------------------------+---------------------+---------------------------------------------+
| min_cov                    | "2"                 | Minimum contig coverage                     |
+----------------------------+---------------------+---------------------------------------------+



Outputs
=======


01_QC:
------

MultiQC report of Illumina pre-processing fastq files::
   
   pre_processing_multiqc_report.html

Nanopore report of each nanopore fastq file::

   {sample}_nanopore_report.html

2_CLEAN_DATA:
-------------

Clean Illumina fastq files (not normalized) for subsampled and total reads::

   {sample}_unpaired_clean.sub.fastq
   {sample}_unpaired_clean.tot.fastq

Clean Nanopore fastq files (not normalized) for subsampled and total reads::

   {sample}_clean.sub.fastq
   {sample}_clean.tot.fastq

3_CONTIGS:
----------


Assembled Spades scaffolds::

   {sample}_spades_filtered_scaffolds.sub.fasta
   {sample}_spades_filtered_scaffolds.tot.fasta

Assembled Canu scaffolds, error corrected with Pilon::

   {sample}_canu_filtered_scaffolds.sub.fasta
   {sample}_canu_filtered_scaffolds.tot.fasta

Quast report of (Spades) or (Spades + Canu)  assembly::

   {sample}_quast_report.sub.txt
   {sample}_quast_report.tot.txt

4_vOTUs:
--------

All assembled scaffolds::

   merged_scaffolds.sub.fasta
   merged_scaffolds.tot.fasta

Representatives contigs after clustering::

   merged_scaffolds.sub_95-80.fna
   merged_scaffolds.tot_95-80.fna

5_VIRAL_ID:
-----------

High confidence viral contigs::

   high_confidence.sub.fasta
   high_confidence.tot.fasta

Low confidence viral contigs::

   low_confidence.sub.fasta
   low_confidence.tot.fasta

High confidence open reading frames (ORFs)::

   high_confidence_ORFs.sub.fasta
   high_confidence_ORFs.tot.fasta

Low confidence open reading frames (ORFs)::

   low_confidence_ORFs.sub.fasta
   low_confidence_ORFs.tot.fasta

6_MAPPING:
----------

Sorted bam files of sample clean reads mapped to high confidence contigs::

   {sample}_high_confidence_sorted.sub_filtered.bam
   {sample}_high_confidence_sorted.tot_filtered.bam

Sorted bam files of sample clean reads mapped to low confidence contigs::

   {sample}_low_confidence_sorted.sub_filtered.bam
   {sample}_low_confidence_sorted.tot_filtered.bam

Indexed bam files of sample clean reads mapped to high confidence contigs::

   {sample}_high_confidence_sorted.sub_filtered.bam.bai
   {sample}_high_confidence_sorted.tot_filtered.bam.bai

Indexed bam files of sample clean reads mapped to low confidence contigs::

   {sample}_low_confidence_sorted.sub_filtered.bam.bai
   {sample}_low_confidence_sorted.tot_filtered.bam.bai

Abundance tables in txt format::

   high_confidence_vOTU_abundance_table.sub.txt
   high_confidence_vOTU_abundance_table.tot.txt
   low_confidence_vOTU_abundance_table.sub.txt
   low_confidence_vOTU_abundance_table.tot.txt

Abundance tables in BIOM JSON format::

   high_confidence_vOTU_abundance_table_json.sub.biom
   high_confidence_vOTU_abundance_table_json.tot.biom
   low_confidence_vOTU_abundance_table_json.sub.biom
   low_confidence_vOTU_abundance_table_json.tot.biom

Contigs summary table::

   vOTU_summary.txt



Support
=======

The easiest way to get help with the project is to send an email to lm.forero10@uniandes.edu.co.
The other good way is to open an issue on Mosaic's Github page at https://github.com/lauramilena3/Mosaic/issues.

.. _https://github.com/lauramilena3/Mosaic/issues: https://github.com/lauramilena3/Mosaic/issues


