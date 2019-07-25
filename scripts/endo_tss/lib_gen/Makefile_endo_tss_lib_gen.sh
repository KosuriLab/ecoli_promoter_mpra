tss_only_wanner.txt: parse_wanner_annot.py
	@ echo "Parsing Wanner annotation..."
	@ python parse_annot.py GSE52059_wanner_txn_annot.txt $@


final_tss_list.txt: tss_only_wanner.txt
	@ echo "Creating TSS list..."
	@ Rscript parse_tss.R


endo_lib_2016.txt: final_tss_list.txt
	@ echo "Generating endo TSS library..."
	@ python lib_gen_tss.py ../../../ref/Escherichia_coli_K-12_MG1655.fasta \
	final_tss_list.txt fwd_primers.fasta rev_primers.fasta $@ tab


# prepare control sequences
control_sequences.txt: endo_lib_2016.txt
	@ echo "Generating control sequences..."
	@ python generate_control_promoters.py synthetic_promoter.csv features.tsv \
	../../../ref/Escherichia_coli_K-12_MG1655.fasta $@


# concatenate control sequences to library, remove duplicates, convert to FASTA
endo_lib_2016_controls.fasta: control_sequences.txt endo_lib_2016.txt
	@ cat $^ > tmp.txt
	@ Rscript remove_duplicates.R tmp.txt endo_lib_2016_controls.txt
	@ python tab_to_fasta.py endo_lib_2016_controls.txt $@
	@ rm tmp.txt


endo_lib_2016_controls_clean.txt: endo_lib_2016_controls.fasta
	@ echo "Cleaning up RE sites..."
	@ python cleanup_REs.py $< RE.fasta $@ tab

