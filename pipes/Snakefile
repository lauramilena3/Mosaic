import os
import re
#import pandas as pd
#======================================================
# Config files
#======================================================
configfile: "config.yaml"

#======================================================
# Global variables
#======================================================

RAW_DATA_DIR =config['input_dir'].rstrip("/") 
RESULTS_DIR=config['results_dir'].rstrip("/")
dir_list = ["RULES_DIR","ENVS_DIR", "ADAPTERS_DIR", "CONTAMINANTS_DIR","RAW_DATA_DIR", "QC_DIR", "CLEAN_DATA_DIR", "ASSEMBLY_DIR", "vOUT_DIR", "VIRAL_DIR", "MAPPING_DIR", "PROFILE_DIR", "MERGE_DIR"]
dir_names = ["rules", "../envs", "db/adapters",  RESULTS_DIR + "/db/contaminants" ,RAW_DATA_DIR, RESULTS_DIR + "/01_QC", RESULTS_DIR + "/02_CLEAN_DATA",RESULTS_DIR + "/03_CONTIGS", RESULTS_DIR + "/04_vOTUs", RESULTS_DIR + "/05_VIRAL_ID",RESULTS_DIR + "/06_MAPPING", "05_ANVIO_PROFILE", "06_MERGED"]
dirs_dict = dict(zip(dir_list, dir_names))
RULES_DIR = 'rules'
CONFIDENCE_TYPES=["high", "low"]
SAMPLING_TYPE=["sub", "tot"]
SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{sample}_" + config['forward_tag'] + ".fastq")

NANOPORE_SAMPLES=SAMPLES
CONTAMINANTS=config['contaminants_list'].split()
NANOPORE=False
PAIRED=False
POOLED=config['nanopore_pooled']



READ_TYPES=[config['forward_tag']]
for fname in os.listdir(RAW_DATA_DIR):
	if fname.endswith(config['reverse_tag'] + '.fastq'):
		PAIRED=True
	elif fname.endswith(config['nanopore_tag'] + '.fastq'):
		NANOPORE=True
if PAIRED:
	READ_TYPES.append(config['reverse_tag'])
if POOLED == "True":
	NANOPORE_SAMPLES=config['nanopore_pooled_name']
if len(SAMPLES)==1:
	SAMPLING_TYPE=["tot"]

print(READ_TYPES)
print(SAMPLES)
print(CONTAMINANTS)

#======================================================
# Rules
#======================================================
 
def inputAll(wildcards):
	inputs=[]
	inputs.append(dirs_dict["QC_DIR"]+ "/preQC_multiqc_report.html")
	#inputs.append(dirs_dict["QC_DIR"]+ "/postQC_multiqc_report.html")
	inputs.extend(expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_quast_report.{sampling}.txt", sample=SAMPLES, sampling=SAMPLING_TYPE))
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/vOTU_abundance_table_json.{sampling}.biom", sampling=SAMPLING_TYPE, confidence=CONFIDENCE_TYPES))
	inputs.extend(expand(dirs_dict["MAPPING_DIR"]+ "/vOTU_summary.{sampling}.txt",sampling=SAMPLING_TYPE))
	if NANOPORE=="True":
		inputs.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_preQC.html",sample_nanopore=NANOPORE_SAMPLES))
		inputs.extend(expand(dirs_dict["QC_DIR"] + "/{sample_nanopore}_nanopore_report_postQC.html", sample_nanopore=NANOPORE_SAMPLES))
	return inputs

rule all:
	input:
		#THA RULES!
		inputAll,
		#TESTING:


# rule qc:ç
# 	input:
# 		dirs_dict["QC_DIR"]+ "/pre_processing_multiqc_report.html"
# 		expand(dirs_dict["QC_DIR"] + "/{sample}_nanopore_report.html",sample=NANOPORE_SAMPLES)
# 		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_norm.{sampling}.fastq", sample=SAMPLES, sampling=SAMPLING_TYPE)
# 		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_nanopore_clean.{sampling}.fastq", sample=SAMPLES, sampling=SAMPLING_TYPE)
# 		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample_nanopore}_nanopore_clean.sub.fastq",sample=NANOPORE_SAMPLES)

# rule assembly:
# 	input:






include: os.path.join(RULES_DIR, 'qualityControl.smk')	
include: os.path.join(RULES_DIR, 'assembly.smk')
include: os.path.join(RULES_DIR, 'vOTUclustering.smk')
include: os.path.join(RULES_DIR, 'viralFiltering.smk')
include: os.path.join(RULES_DIR, 'abundance.smk')
include: os.path.join(RULES_DIR, 'taxonomyAssignment.smk')
include: os.path.join(RULES_DIR, 'resultsParsing.smk')





