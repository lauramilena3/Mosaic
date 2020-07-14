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





