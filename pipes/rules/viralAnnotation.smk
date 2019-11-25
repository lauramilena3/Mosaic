# cd TOOLS
# git clone https://github.com/EGTortuero/viga
# chmod 744 viga/create_dbs.sh viga/VIGA.py
# ./viga/create_dbs.sh
#
# cd pilercr1.06
# make
# PILER=/home/lmf/apps/Mosaic/pipes/tools/pilercr1.06
# PATH=$PILER:$PATH
# TRF=/home/lmf/apps/Mosaic/pipes/tools/TRF
# PATH=$TRF:$PATH
#
# cd VIGA
#
# /home/lmf/apps/Mosaic/pipes/tools/viga/VIGA.py --input /home/lmf/03_COLIPHAGES/FASTA/similar.fasta --diamonddb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/RefSeq_Viral_DIAMOND/ --blastdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/RefSeq_Viral_BLAST/ --hmmerdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/pvogs/pvogs.hmm --rfamdb /home/lmf/apps/Mosaic/pipes/tools/viga/databases/rfam/Rfam.cm --modifiers modifiers.txt --threads 16

rule fasta_to_db:
	input:
		fasta= "{directory}" + "/{fasta}.fasta"
	output:
		index= dirs_dict["ASSEMBLY_DIR"] + "/{fasta}.index"
	message:
		"Creating reference and assembly mmseqs db"
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	threads: 2
	shell:
		"""
		mmseqs createdb {input.fasta} {wildcards.fasta}
		"""
