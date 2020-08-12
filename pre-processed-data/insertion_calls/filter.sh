bedtools intersect -u -F 1 -f 1 -a 0.5.21.unfiltered_calls_delta_guts.gff -b manual_call_results/Delta_calls_guts.csv.bed > 0.5.21.filtered_calls_delta_guts.gff
bedtools intersect -u -F 1 -f 1 -a 0.5.21.unfiltered_calls_delta_heads.gff -b manual_call_results/Delta_calls_heads.csv.bed > 0.5.21.filtered_calls_delta_heads.gff
bedtools intersect -u -F 1 -f 1 -a 0.5.21.unfiltered_calls_guts.gff -b manual_call_results/ProsGFP_calls_guts.csv.bed > 0.5.21.filtered_calls_prosgfp_guts.gff
bedtools intersect -u -F 1 -f 1 -a 0.5.21.unfiltered_calls_heads.gff -b manual_call_results/ProsGFP_calls_heads.csv.bed > 0.5.21.filtered_calls_prosgfp_heads.gff
