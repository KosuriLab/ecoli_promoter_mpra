MERGE_DEST=../../processed_data/endo_tss/lb/merged_reads
RAW_DATA=../../raw_data/endo_tss/lb
DATA=../../processed_data/endo_tss/lb
LIB=../../processed_data/lib_gen
DATE=$(shell date +%Y%m%d)

all: map map_controls format standardize bed
	
expression: $(DATA)/rLP5_Endo2_lb_expression.txt

merge: $(MERGE_DEST)/endo_lb_merged.fastq

map: $(DATA)/endo_lb_bc_map.txt

map_controls: $(DATA)/endo_lb_controls_bc_map.txt

barcode_counts: $(DATA)/counts_rLP5_Endo2_DNA1.txt

format: $(DATA)/rLP5_Endo2_lb_expression_formatted.txt

standardize: $(DATA)/rLP5_Endo2_lb_expression_formatted_std.txt

bed: $(DATA)/lb/negative_TSS_plus.bed

# mapping
$(MERGE_DEST)/pEndo2_84.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndo2_84_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pEndo2_84_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndo2_84 $(MERGE_DEST)


$(MERGE_DEST)/pEndo2_85.paired.assembled.fastq: $(RAW_DATA)/bc_map/pEndo2_85_R1.fastq.gz \
	$(RAW_DATA)/bc_map/pEndo2_85_R2.fastq.gz ../pearallel.sh
	../pearallel.sh $(RAW_DATA)/bc_map/pEndo2_85 $(MERGE_DEST)


$(MERGE_DEST)/endo_lb_merged.fastq: $(MERGE_DEST)/pEndo2_84.paired.assembled.fastq $(MERGE_DEST)/pEndo2_85.paired.assembled.fastq
	@ cat $^ > $@


# # trim reference to remove primer: awk '!/>/{ print substr($1, 25, 150);next}1' ref.fasta
# $(DATA)/endo_lb_mapping_barcode_statistics.txt: $(MERGE_DEST)/endo_lb_merged.fastq
# 	@ python ../barcode_mapping.py $< fastq ../../ref/endo_lib_2016_controls_clean.txt \
# 	tab 150 24 end 20 $(DATA)/endo_lb_mapping --controls ../../ref/synthetic_promoter_pos_controls_formatted.txt \
# 	--controls_primer_len 0

$(DATA)/endo_lb_bc_map.txt: $(MERGE_DEST)/endo_lb_merged.fastq
	@ python ../bc_map_consensus.py $< ../../ref/endo_lib_2016_controls_clean.fasta \
	--bc-start -20 --bc-length 20 --proc 30 --verbose --bbmap-procs 30 \
	--var_start 1 --var_length 150 --min-reads 3 --bad-bcs $(DATA)/bc_map_consensus_bad_bcs.txt > $@


$(DATA)/endo_lb_controls_bc_map.txt: $(MERGE_DEST)/endo_lb_merged.fastq
	@ python ../barcode_mapping_ecoli_pos_controls.py $<  fastq \
	../../ref/synthetic_promoter_pos_controls_formatted.fasta 0 end 20 $(DATA)/endo_lb_controls
	@ rm -f $(DATA)/endo_lb_controls_bbmap.sam

# RNA barcode counting

# use one FASTQ file to trigger barcode counting for all FASTQ files
$(DATA)/counts_rLP5_Endo2_DNA1.txt: | $(RAW_DATA)/DNA_RNA/rLP5_Endo2_DNA1.fastq.gz
	@ ../count_barcodes.sh $(RAW_DATA)/DNA_RNA
	@ mv counts* $(DATA)

# # calculate expression
# $(DATA)/rLP5_Endo2_lb_expression.txt: $(DATA)/counts_rLP5_Endo2_DNA1.txt \
# 	$(DATA)/endo_lb_mapping_barcode_statistics.txt
# 	@ Rscript bc_counts2expression.R $(DATA) $(DATA)/endo_lb_mapping_barcode_statistics.txt \
# 	$(DATA)/endo_lb_mapping_variant_statistics.txt ../../ref/endo_lib_2016_controls_clean.txt \
# 	$@

$(DATA)/rLP5_Endo2_lb_expression.txt: $(DATA)/endo_lb_bc_map.txt $(DATA)/endo_lb_controls_bc_map.txt
	@ Rscript --vanilla ../barcode_counts_to_expression.R $(DATA) $^ ../../ref/endo_lib_2016_controls_clean.txt $@


# format library
$(DATA)/rLP5_Endo2_lb_expression_formatted.txt: $(DATA)/rLP5_Endo2_lb_expression.txt
	@ Rscript --vanilla endo_tss_format.R $< $@


# standardize, all other libraries must be formatted first
$(DATA)/rLP5_Endo2_lb_expression_formatted_std.txt: $(DATA)/rLP5_Endo2_lb_expression_formatted.txt
	@ echo "Standardizing expression files..."
	@ cd ..
	@ Rscript --vanilla standardize_expression.R


# generate bed files for downstream analysis
$(DATA)/lb/negative_TSS_plus.bed: $(DATA)/rLP5_Endo2_lb_expression_formatted_std.txt
	@ echo "Creating bed files ..."
	@ ./create_tss_bed.sh $<