.. image:: https://user-images.githubusercontent.com/12174880/79860360-d483ce00-83d2-11ea-93ff-8030ff58fc24.png?raw=true
   :width: 200px
   :height: 100px
   :scale: 50 %
   :alt: alternate text
   :align: right

.. _getting_started:

Requirements
============

The only required software for running Mosaic is Anaconda. Please follow the installation guide avaliable `here <https://docs.anaconda.com/anaconda/install/>`_ .

Installation
============

Clone GitHub  repo and enter directory::
   
   git clone https://github.com/lauramilena3/Mosaic
   cd Mosaic/workflow

Create Mosaic virtual environment and activate it::
   
   conda env create -n Mosaic -f Mosaic.yaml
   source activate Mosaic

Running Mosaic
==============

View the number of avaliable cores with::
   
   nproc #Linux

Create variables for the number of cores, your raw data directory and the results directory::
   
   nCores=16
   in_dir="/path/to/your/raw/data"
   out_dir="/path/to/your/desired/results/dir"

Run Mosaic's pipeline with the desired number of cores and choosen directories::
   
   snakemake -j $nCores --use-conda --config input_dir=$in_dir results_dir=$out_dir

NOTE: Every time you run Mosaic you need to: 1) activate the virtual environment and 2) run it from the Mosaic/workflow folder.

Visualize the workflow 
+++++++++++++++++++++++

DAG
***

The directed acyclic graph (DAG) includes all files and rules used in the workflow.

To view the DAG use::

   snakemake --dag | dot -Tpdf > dag.pdf

RULES
*****

If you want to visualize a simpler diagram, only including the rules, use::

   snakemake --rulegraph | dot -Tpdf > rulegraph.pdf

Mosaic modules
==============

Other mosaic submodules::

   snakemake --use-conda --config input_dir=$in_dir/00_RAW_READS/ -j 16 -k
   
   snakemake --use-conda --config input_dir=$in_dir/00_RAW_READS/ nanopore_pooled="True" nanopore_pooled_name=direct_virome nanopore_quality=7 -k -j 64 -n

   snakemake -j 32 --use-conda -p clean_reads --config input_dir=$in_dir/00_RAW_DATA/ contaminants_list="GCF_000001405.39"

   snakemake -j 16 --use-conda -p assembly_vs_reference --config results_dir=$output_dir/ VIRAL_CONTIGS=$ncbi.fasta

   snakemake --use-conda -p annotate_VIGA --config representative_contigs=$contigs_file.tot.fasta -k -j 16 -n

   snakemake --use-conda -p annotate_contigs --config representative_contigs=$contigs_file.tot.fasta -k -j 16 -n

   snakemake --use-conda -p annotate_VIBRANT_contigs --config representative_contigs=$contigs_file.tot.fasta -k -j 16 -n

   snakemake --use-conda -p abundance_from_db_contigs --config input_dir=in_dir/00_RAW_READS/  representative_contigs=$contigs_file.tot.fasta results_dir=$output_dir/

   snakemake --use-conda -k -j 16 -p taxonomyAssignmentvContact --config representative_contigs=$contigs_file.tot.fasta 

   snakemake --use-conda -p taxonomyAssignmenMMseqs --config representative_contigs=$contigs_file.tot.fasta  -k -j 16






