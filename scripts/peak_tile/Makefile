MERGE_DEST=../../processed_data/peak_tile/merged_reads
RAW_DATA=../../raw_data/peak_tile
DATA=../../processed_data/peak_tile


all: map map_controls format

merge: $(MERGE_DEST)/peak_tile_combined_merged.fastq

map: $(DATA)/peak_tile_combined_bc_map.txt

map_controls: $(DATA)/peak_tile_combined_controls_bc_map.txt

barcode_counts: $(DATA)/counts_rLP5_PeakTiling_LB_DNA1_1_S5_R1_001.txt

expression: $(DATA)/peak_tile_expression.txt

format: $(DATA)/peak_tile_expression_formatted.txt

split: $(DATA)/peak_tile_expression_90train.txt


.IGNORE: ../../ref/20180508_lb_peak_tile_lib.txt $(DATA)/peak_tile_expression_formatted_std.txt

# map combined reads
$(MERGE_DEST)/pPeakTiling_82_combined.paired.assembled.fastq: $(RAW_DATA)/bc_map/pPeakTiling_82_combined_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pPeakTiling_82_combined_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/clearpPeakTiling_82_combined $(MERGE_DEST)


$(MERGE_DEST)/pPeakTiling_83_combined.paired.assembled.fastq: $(RAW_DATA)/bc_map/pPeakTiling_83_combined_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pPeakTiling_83_combined_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pPeakTiling_83_combined $(MERGE_DEST)


$(MERGE_DEST)/peak_tile_combined_merged.fastq: $(MERGE_DEST)/pPeakTiling_82_combined.paired.assembled.fastq $(MERGE_DEST)/pPeakTiling_83_combined.paired.assembled.fastq
	@ cat $^ > $@


$(DATA)/peak_tile_combined_bc_map.txt: $(MERGE_DEST)/peak_tile_combined_merged.fastq
	@ python ../bc_map_consensus.py $< ../../ref/20180508_lb_peak_tile_lib_trimmed.fasta \
	--bc-start -20 --bc-length 20 --proc 30 --verbose --bbmap-procs 30 \
	--var_start 1 --var_length 150 --min-reads 3 --bad-bcs $(DATA)/bc_map_consensus_bad_bcs.txt > $@


# peaks < 150bp were synthesized with flanking stuffer sequence and must be mapped
# with positive controls, create combined fasta file
../../ref/peak_tile_short_and_pos_controls.fasta: ../../ref/20180508_lb_peak_tile_lib.txt
	@ Rscript --vanilla format_short_peaks.R
	@ python ../format_converter.py ../../ref/20180508_lb_peak_tile_lib_short_peaks_trimmed.txt \
	tab fasta ../../ref/20180508_lb_peak_tile_lib_short_peaks_trimmed.fasta
	@ cat ../../ref/20180508_lb_peak_tile_lib_short_peaks_trimmed.fasta \
	../../ref/synthetic_promoter_pos_controls_formatted.fasta > $@


$(DATA)/peak_tile_combined_controls_bc_map.txt: $(MERGE_DEST)/peak_tile_combined_merged.fastq
	@ python ../barcode_mapping_ecoli_pos_controls.py $< fastq \
	../../ref/peak_tile_short_and_pos_controls.fasta 0 end 20 $(DATA)/peak_tile_combined_controls
	@ rm -f $(DATA)/peak_tile_combined_controls_bbmap.sam
	

# use one FASTQ file to trigger barcode counting
$(DATA)/counts_rLP5_PeakTiling_LB_DNA1_1_S5_R1_001.txt: $(RAW_DATA)/DNA_RNA/rLP5_PeakTiling_LB_DNA1_1_S5_R1_001.fastq.gz
	@ ../count_barcodes.sh $(RAW_DATA)/DNA_RNA
	@ mv counts* $(DATA)


$(DATA)/peak_tile_expression.txt: $(DATA)/peak_tile_combined_bc_map.txt $(DATA)/peak_tile_combined_controls_bc_map.txt
	@ Rscript --vanilla ../barcode_counts_to_expression.R $(DATA) $^ ../../ref/20180508_lb_peak_tile_lib.txt $@


$(DATA)/peak_tile_expression_formatted.txt: $(DATA)/peak_tile_expression.txt
	@ Rscript --vanilla peak_tile_format.R $< $@


############ Data for figures, depends on standardized expression #############

$(DATA)/peak_tile_expression_90train.txt: $(DATA)/peak_tile_expression_formatted_std.txt
	@ echo "Splitting into test and train..."
	@ Rscript --vanilla peak_tile_train_test_split.R

