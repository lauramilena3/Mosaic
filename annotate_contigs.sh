#/bin/bash
source /home/lmf/.local/opt/anaconda3/etc/profile.d/conda.sh
cwd=$(pwd)
fasta=$1
fasta_name=$(basename "$fasta" .fasta)
cp $fasta ${fasta_name}.tot.fasta
conda activate Mosaic
cd /home/lmf/apps/Mosaic/workflow/
snakemake --use-conda -p annotate_VIGA --config representative_contigs=${cwd}/${fasta_name}.tot.fasta -k -j 16
cd $cwd
mv 07_ANNOTATION annotation_VIGA

#Get proteins from VIGA annotation
python2 /home/lmf/scripts/gbk_to_faa.py annotation_VIGA/${fasta_name}.tot_annotated.gbk ${fasta_name}.tot.faa

#Annotation with HHpred
sed -i 's/\s.*$//' ${fasta_name}.tot.faa
mkdir annotation_hhpred
ln -s ${cwd}/${fasta_name}.tot.faa annotation_hhpred
cd annotation_hhpred
sed -i 's/ .*//' ${fasta_name}.tot.faa

cat ${fasta_name}.tot.faa | awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename}'

hh_annotation () {
	local aa=$1
	name="${aa%.fa}"
	echo $name
	/home/lmf/scripts/reformat.pl fas a2m ${name}.fa ${name}.a2m
	hhblits -i ${name}.a2m -d /home/lmf/db/uniclust/UniRef30_2020_02 -oa3m MSA_uniref_${name}.a3m -norealign -n 3 -e 1e-3 -qid 0 -cov 20 -cpu 2 -o temp_results${name}.hhr
	hhsearch -i MSA_uniref_${name}.a3m -d /opt/hh-suite/data/pdb70 -d /opt/hh-suite/data/pfam -d /opt/hh-suite/data/scop70_1.75  -d /home/lmf/db/NCBI_CD/NCBI_CD  -o results_${name}.hhr -oa3m results_${name}.a3m -p 20 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -contxt /opt/hh-suite/data/context_data.crf -cpu 2
}
export -f hh_annotation
ls *fa | parallel --lb --jobs 32 hh_annotation {}

for hhr in ./results_*hhr; do name="${hhr%.hhr}"; echo $name; grep "^  [0-9] " $hhr > ${name}.txt; done

#Annotation with BLAST 2.9.0+ viral_dir
cd ..
mkdir annotation_blast
ln -s ${cwd}/${fasta_name}.tot.faa annotation_blast
cd annotation_blast
blastp -num_threads 32 -db /home/lmf/db/NCBI_viral_proteins/NCBI_virus_proteins_Jun_20_2020.faa -query ${fasta_name}.tot.faa -outfmt "6 qseqid sseqid stitle qstart qend qlen slen qcovs evalue length" > blast.out
python /home/lmf/scripts/annotate_contigs.py ${fasta_name}
