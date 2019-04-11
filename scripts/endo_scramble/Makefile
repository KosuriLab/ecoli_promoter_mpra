MERGE_DEST=../../processed_data/endo_scramble/merged_reads
RAW_DATA=../../raw_data/endo_scramble
DATA=../../processed_data/endo_scramble


all: map map_controls format

merge: $(MERGE_DEST)/endo_scramble_combined_merged.fastq

map: $(DATA)/endo_scramble_combined_bc_map.txt

map_controls: $(DATA)/endo_scramble_combined_controls_bc_map.txt

# use one file to trigger counting
barcode_count: $(DATA)/counts_rLP5_EndoScramble_LB_DNA1_1_S13_R1_001.txt 

expression: $(DATA)/endo_scramble_expression.txt

format: $(DATA)/endo_scramble_expression_formatted.txt


# Data from two separate mapping runs 
#############ls p################### Read merging ################################
$(MERGE_DEST)/pEndoScramble_80_combined.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_80_combined_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pEndoScramble_80_combined_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_80_combined $(MERGE_DEST)


$(MERGE_DEST)/pEndoScramble_81_combined.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndoScramble_81_combined_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pEndoScramble_81_combined_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndoScramble_81_combined $(MERGE_DEST)


$(MERGE_DEST)/endo_scramble_combined_merged.fastq: $(MERGE_DEST)/pEndoScramble_80_combined.paired.assembled.fastq $(MERGE_DEST)/pEndoScramble_81_combined.paired.assembled.fastq
	@ cat $^ > $@


################################ Mapping ####################################
$(DATA)/endo_scramble_combined_bc_map.txt: $(MERGE_DEST)/endo_scramble_combined_merged.fastq
	@ python ../bc_map_consensus.py $< ../../ref/20180507_active_tss_scrambled10_stride5_trimmed.fasta \
	--bc-start -20 --bc-length 20 --proc 30 --verbose --bbmap-procs 30 \
	--var_start 1 --var_length 150 --min-reads 3 --bad-bcs $(DATA)/bc_map_consensus_bad_bcs.txt > $@


$(DATA)/endo_scramble_combined_controls_bc_map.txt: $(MERGE_DEST)/endo_scramble_combined_merged.fastq
	@ python ../barcode_mapping_ecoli_pos_controls.py $< fastq \
	../../ref/synthetic_promoter_pos_controls_formatted.fasta 0 end 20 $(DATA)/endo_scramble_combined_controls
	@ rm -f $(DATA)/endo_scramble_controls_bbmap.sam
	@ rm -fr $(DATA)/ref


############################ Barcode counting ##################################

$(DATA)/counts_rLP5_EndoScramble_LB_DNA1_1_S13_R1_001.txt: $(RAW_DATA)/DNA_RNA/rLP5_EndoScramble_LB_DNA1_1_S13_R1_001.fastq.gz
	@ ../count_barcodes.sh $(RAW_DATA)/DNA_RNA
	@ mv counts* $(DATA)


################################ Expression ####################################
$(DATA)/endo_scramble_expression.txt: $(DATA)/endo_scramble_combined_bc_map.txt \
	$(DATA)/endo_scramble_combined_controls_bc_map.txt
	@ Rscript --vanilla ../barcode_counts_to_expression.R $(DATA) $^ ../../ref/20180507_active_tss_scrambled10_stride5.txt $@

$(DATA)/endo_scramble_expression_formatted.txt: $(DATA)/endo_scramble_expression.txt
	@ Rscript --vanilla endo_scramble_format.R $< $@
