Requirements
============

- Anaconda

You can follow the `installation guide <https://docs.anaconda.com/anaconda/install/>`_ .

Installation
============

Clone GitHub  repo and enter directory::
   
   git clone https://github.com/lauramilena3/Mosaic
   cd On-rep-seq

Create On-rep-seq virtual environment and activate it::
   
   conda env create -n Mosaic -f Mosaic.yaml
   source activate Mosaic


Running Mosaic
==============

View the number of avaliable cores with::
   
   nproc #linux
   sysctl -n hw.ncpu #Os

Go into Mosaic directory and create variables for the number of cores,
your basecalled data directory and the results directory of your choice::
   
   nCores="cores"
   fastqDir="/path/to/your/basecalled/data"
   reusultDir="/path/to/your/desired/results/dir"

Run Mosaic's pipeline with the desired number of cores::
   
   snakemake -j $nCores --use-conda --config basecalled_data=$fastqDir results_dir=$reusultDir

If you are using your laptop we suggest you to leave 2 free processors
for other system tasks. 

View dag of jobs to visualize the workflow 
++++++++++++++++++++++++++++++++++++++++++

To view the dag run::

   snakemake --dag | dot -Tpdf > dag.pdf





