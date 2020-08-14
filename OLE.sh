https://schaechter.asmblog.org/schaechter/2015/09/phage-responses-to-the-sos-response.html


hmmsearch --tblout DUF99.tbl DUF99.hmm /home/lmf/db/IMG_VR/IMG_VR.aa > DUF99.out
grep -h -v "^#" DUF99.tbl | awk '{print $1}' | sort | uniq | sort > positive_domain_GGDEF.csv
seqtk subseq /home/lmf/db/IMG_VR/IMG_VR.aa positive_domain_GGDEF.csv > positive_domain_GGDEF.fasta
mmseqs easy-cluster positive_domain_GGDEF.fasta clusterRes_GGDEF tmp --min-seq-id 0.80 -c 0.85 --cov-mode 0
#positive_domain_GGDEF.fasta:		71
#clusterRes_GGDEF_rep_seq.fasta:	34

hmmsearch --tblout EAL.tbl EAL.hmm /home/lmf/db/IMG_VR/IMG_VR.aa > EAL.out
grep -h -v "^#" EAL.tbl | awk '{print $1}' | sort | uniq | sort > positive_domain_EAL.csv
seqtk subseq /home/lmf/db/IMG_VR/IMG_VR.aa positive_domain_EAL.csv > positive_domain_EAL.fasta
mmseqs easy-cluster positive_domain_EAL.fasta clusterRes_EAL tmp --min-seq-id 0.80 -c 0.85 --cov-mode 0
#positive_domain_EAL.fasta:			295
#clusterRes_EAL_rep_seq.fasta:		210

grep -h -v "^#" *tbl | awk '{print $1}' | sort | uniq -c | sort

mkdir hhpred
ln -s /home/lmf/OLE/clusterRes_rep_seq.fasta .
sed -i 's/\s.*$//' clusterRes_rep_seq.fasta
cat clusterRes_rep_seq.fasta  | awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename}'
