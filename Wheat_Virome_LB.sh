demeter: /media/demeter/Seagate Expansion Drive/20200703_FS10000714_10_BPL20321-3127/Data/Intensities/BaseCalls/0/Wheat_virome/

cd /home/lmf/apps/Mosaic/workflow
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/00_RAW_READS min_len="0" min_cov="0" -k -j 16 -n min_len="1000"
snakemake --use-conda --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/newMosaic/00_RAW_READS -k -j 16 -n


#CONTAMINATION
cd /home/lmf/apps/Mosaic/workflow
snakemake --use-conda -p abundance_from_db_contigs --config input_dir=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment/00_RAW_READS representative_contigs=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment/Triticum_aestivum.tot.fasta results_dir=/home/lmf/PhylloVir/Wheat_Virome/ContaminationAssesment
