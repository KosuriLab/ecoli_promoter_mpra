DATA=$1

# $5 and $6 are genomic start and end, $4 is strand, $17 is category, $20 is active
# active TSS plus
cat $DATA | awk -v OFS='\t' '{if (($4 == "+") && ($20 == "active")) {print "U00096.2", $5, $6, $1, ".", $4} }' > \
../../processed_data/endo_tss/lb/active_TSS_plus.bed

# inactive TSS plus
cat $DATA | awk -v OFS='\t' '{if (($4 == "+") && ($20 == "inactive")) {print "U00096.2", $5, $6, $1, ".", $4} }' > \
../../processed_data/endo_tss/lb/inactive_TSS_plus.bed

# active TSS plus
cat $DATA | awk -v OFS='\t' '{if (($4 == "-") && ($20 == "active")) {print "U00096.2", $5, $6, $1, ".", $4} }' > \
../../processed_data/endo_tss/lb/active_TSS_minus.bed

# inactive TSS plus
cat $DATA | awk -v OFS='\t' '{if (($4 == "-") && ($20 == "inactive")) {print "U00096.2", $5, $6, $1, ".", $4} }' > \
../../processed_data/endo_tss/lb/inactive_TSS_minus.bed

# negative TSS plus
cat $DATA | awk -v OFS='\t' '{if ($17 == "neg_control") {print "U00096.2", $5, $6, $1, ".", $4} }' > \
../../processed_data/endo_tss/lb/negative_TSS_plus.bed