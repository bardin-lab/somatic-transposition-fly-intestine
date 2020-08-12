import papermill as pm


gff_files = [
    '0.5.21.filtered_calls_delta_guts.gff',
    '0.5.21.filtered_calls_delta_heads.gff',
    '0.5.21.filtered_calls_prosgfp_guts.gff',
    '0.5.21.filtered_calls_prosgfp_heads.gff',
    '../nanopore_guts_combined.bed',
    '../nanopore_heads_combined.bed',
]
INSERTION_SET = ['Delta guts', 'Delta heads', 'ProsGFP guts', 'ProsGFP heads', 'ProsGFP guts nanopore', 'ProsGFP heads nanopore']

for gff_file, insertion_set in zip(gff_files, INSERTION_SET):
    print("running notebook for " + gff_file)
    pm.execute_notebook('gene_overlap_data.ipynb', "%s_notebook.ipynb" % insertion_set.lower().replace(' ', '_'), parameters=dict(all_insertions=gff_file, INSERTION_SET=insertion_set))
