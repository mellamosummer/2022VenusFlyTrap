#to create text file compare DEG lists to BLAST hits
paste number_of_genes_DEG_lists.txt number_of_annotated_BLAST_hits_DEG_lists.txt > total_genes_per_list_vs_BLAST_hits.txt


**BLAST against Tair database for GO annotations**

#download TAIR sequences
wget https://www.arabidopsis.org/download_files/Proteins/Araport11_protein_lists/Araport11_pep_20220103.gz

#only keep ID and sequence
sed -e 's/^\(>[^[:space:]]*\).*/\1/' Araport11_pep_20220103 > Araport11_pep_20220103.fasta

#make BLASTdb
makeblastdb.sh

#run gridrunner script with -db tairdb
gridrunnerBLAST.sh

#use same format to concatenate out files
find blast_results_title_oneHSP_tairdb -name '*.fa.OUT' | xargs cat > combined_blast_results_tairdb.out

**GO annotations**

#pull AT gene names from BLAST results
awk '{print $1, $2}' results/combined_blast_results_tairdb.out > Dm_AT_BLASThits.txt

#paste entire list into https://www.arabidopsis.org/tools/go_term_enrichment.jsp


**Traps Overlap BLAST**

#swissprot BLAST
grep -f data/trapspreynopreyoverlap.txt results/Dm_proteins_BLAST_combined/combined_blast_results_oneHSP.out > trapsoverlap_swissprot.txt

#tair BLAST
grep -f data/trapspreynopreyoverlap.txt results/Dm_proteins_BLAST_combined/combined_blast_results_tairdb.out > trapsoverlap_tair.txt

**BLAST results headers**

qseqid 
sseqid 
stitle 
pident 
length 
mismatch 
gapopen 
qstart 
qend 
sstart 
send 
evalue 
bitscore
