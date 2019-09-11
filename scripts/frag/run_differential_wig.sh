echo "plus LB over M9"
python ../wig_calculator.py ../../processed_data/frag/lb/plus_frag_pileup.wig \
../../processed_data/frag/m9/plus_frag_pileup_M9.wig division \
../../processed_data/frag/plus_frag_pileup_lb_over_m9.wig

echo "minus LB over M9"
python ../wig_calculator.py ../../processed_data/frag/lb/minus_frag_pileup.wig \
../../processed_data/frag/m9/minus_frag_pileup_M9.wig division \
../../processed_data/frag/minus_frag_pileup_lb_over_m9.wig

echo "plus M9 over LB"
python ../wig_calculator.py ../../processed_data/frag/m9/plus_frag_pileup_M9.wig \
../../processed_data/frag/lb/plus_frag_pileup.wig division \
../../processed_data/frag/plus_frag_pileup_m9_over_lb.wig

echo "minus M9 over LB"
python ../wig_calculator.py ../../processed_data/frag/m9/minus_frag_pileup_M9.wig \
../../processed_data/frag/lb/minus_frag_pileup.wig division \
../../processed_data/frag/minus_frag_pileup_m9_over_lb.wig

echo "Calling peaks..."

python call_peaks.py ../../processed_data/frag/plus_frag_pileup_lb_over_m9.wig \
2.0 10 60 + ../../processed_data/frag/plus_lb_over_m9_peaks.bed

python call_peaks.py ../../processed_data/frag/minus_frag_pileup_lb_over_m9.wig \
2.0 10 60 - ../../processed_data/frag/minus_lb_over_m9_peaks.bed

python call_peaks.py ../../processed_data/frag/plus_frag_pileup_lb_over_m9.wig \
2.0 10 60 + ../../processed_data/frag/plus_m9_over_lb_peaks.bed

python call_peaks.py ../../processed_data/frag/plus_frag_pileup_lb_over_m9.wig \
2.0 10 60 - ../../processed_data/frag/minus_m9_over_lb_peaks.bed