import os
import papermill as pm

guts_combined = 'nanopore_guts_combined.bed'
heads_combined = 'nanopore_heads_combined.bed'

if not os.path.exists(guts_combined):
    nanopore_data = '../1. nanopore insertion analysis/Variant fragments 25 days normalized/all_singletons_with_tsd.bed'
    with open(guts_combined, 'w') as guts, open(heads_combined, 'w') as heads, open(nanopore_data) as heads_and_guts:
        for line in heads_and_guts:
            if 'Guts' in line:
                guts.write(line)
            else:
                heads.write()

gff_files = [
    '0.5.21.filtered_calls_delta_guts.gff',
    '0.5.21.filtered_calls_delta_heads.gff',
    '0.5.21.filtered_calls_prosgfp_guts.gff',
    '0.5.21.filtered_calls_prosgfp_heads.gff',
    guts_combined,
    heads_combined,
]
INSERTION_SET = ['Delta guts', 'Delta heads', 'ProsGFP guts', 'ProsGFP heads', 'ProsGFP guts nanopore', 'ProsGFP heads nanopore']

for gff_file, insertion_set in zip(gff_files, INSERTION_SET):
    print("running notebook for " + gff_file)
    pm.execute_notebook('overlap_modencode_data.ipynb', "%s_notebook.ipynb" % insertion_set.lower().replace(' ', '_'), parameters=dict(all_insertions=gff_file, INSERTION_SET=insertion_set))
