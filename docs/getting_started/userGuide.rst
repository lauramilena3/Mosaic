.. _getting_started:

Requirements
============

The only required software for running Mosaic is Anaconda. Please follow the installation guide avaliable `here <https://docs.anaconda.com/anaconda/install/>`_ .

Installation
============

Clone GitHub  repo and enter directory::
   
   git clone https://github.com/lauramilena3/Mosaic
   cd Mosaic/pipes

Create Mosaic virtual environment and activate it::
   
   conda env create -n Mosaic -f Mosaic.yaml
   source activate Mosaic

Running Mosaic
==============

View the number of avaliable cores with::
   
   nproc #Linux
   sysctl -n hw.ncpu #MacOs

Go into Mosaic directory and create variables for the number of cores, your raw data directory and the results directory of your choice::
   
   nCores="cores"
   fastqDir="/path/to/your/raw/data"
   reusultDir="/path/to/your/desired/results/dir"

Run Mosaic's pipeline with the desired number of cores and choosen directories::
   
   snakemake -j $nCores --use-conda --config basecalled_data=$fastqDir results_dir=$reusultDir

NOTE: Please notice that every time you run Mosaic: 1) you will need to activate the virtual environment and 2) you need to run it from the Mosaic/pipes folder. If you are using your laptop we suggest you to leave 2 free processors for other system tasks. 

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





