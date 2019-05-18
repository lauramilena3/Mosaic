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

   snakemake -j $nCores --use-conda --config basecalled_data=$fastqDir results_dir=$reusultDir min_len=2000

Below, you can find all the variables, and their default values.  


+-------------------+---------------------+---------------------------------------------+
| Name              | Default             | Description                                 |
+===================+=====================+=============================================+
| input_dir         | "data/00_RAW_DATA"  | Input folder with raw data                  |
+-------------------+---------------------+---------------------------------------------+
| results_dir       | "data"              | Directory for output files                  |
+-------------------+---------------------+---------------------------------------------+
| operating_system  | "linux"             | Operating system ([linux], [macOs])         |
+-------------------+---------------------+---------------------------------------------+
| forward_tag       | "R1"                | Forward tag used for Illumina reads         |
+-------------------+---------------------+---------------------------------------------+
| reverse_tag       | "R2"                | Reverse tag used for Illumina reads         |
+-------------------+---------------------+---------------------------------------------+
| nanopore_tag      | "nanopore"          | Tag used for Nanopore reads                 |
+-------------------+---------------------+---------------------------------------------+

nanopore_pooled: "True"
nanopore_pooled_name: "raw"

trimmomatic_leading: "15"
trimmomatic_trailing: "15"
trimmomatic_window_size: "4"
trimmomatic_window_quality: "15"
trimmomatic_minlen: "50"

#bbnorm
min_norm: "2" 
max_norm: "100"
max_subsample: "20000000"

#Spades and canu assembly
min_len: "1000"
min_cov: "2"

contaminants_list: "GCF_000819615.1"
adapters_file: "NexteraPE-PE.fa"
basecalled_dir: "data/basecalled/"

Outputs
=======




Support
=======

The easiest way to get help with the project is to join the #crawler
channel on Freenode.

We hang out there and you can get real-time help with your projects.
The other good way is to open an issue on Github.

The mailing list at `https://groups.google.com/forum/#!forum/crawler`_ is also available for support.

Freenode: `freenode`_ 

Github: `github_page`_

.. _https://groups.google.com/forum/#!forum/crawler: https://groups.google.com/forum/#!forum/crawler
.. _freenode: irc://freenode.net
.. _github_page: http://github.com/example/crawler/issues

