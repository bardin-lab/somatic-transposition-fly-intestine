for notebook in "Density and Rainfallplot.ipynb" "ideograms_and_insertions.ipynb" "somatic_overview.ipynb"
do
    papermill -p INSERTIONS 0.5.21.filtered_calls_delta_heads.gff -p LABEL 'Delta calls heads' -p SAMPLE_FILE Delta_head_samples.csv "$notebook" "Delta heads $notebook"
    papermill -p INSERTIONS 0.5.21.filtered_calls_delta_guts.gff -p LABEL 'Delta calls guts' -p SAMPLE_FILE Delta_gut_samples.csv "$notebook" "Delta guts $notebook"
    papermill -p INSERTIONS 0.5.21.filtered_calls_prosgfp_heads.gff -p LABEL 'ProsGFP calls heads' -p SAMPLE_FILE ProsGFP_head_samples.csv "$notebook" "ProsGFP heads $notebook"
    papermill -p INSERTIONS 0.5.21.filtered_calls_prosgfp_guts.gff -p LABEL 'ProsGFP calls guts' -p SAMPLE_FILE ProsGFP_gut_samples.csv "$notebook" "ProsGFP gut $notebook"
done
