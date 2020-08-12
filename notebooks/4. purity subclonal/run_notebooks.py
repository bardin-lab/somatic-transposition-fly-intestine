import papermill as pm


gff_files = [
    '0.5.21.filtered_calls_delta_guts.gff',
    '0.5.21.filtered_calls_prosgfp_guts.gff',
]
DELETION_READ_RATIOS = [
    'delta_purity_estimates_deletion.tab',
    'pros_purity_estimates_deletion.tab',
]
INSERTION_SET = ['Delta guts', 'ProsGFP guts']


for gff_file, insertion_set, deletion_read_ratio in zip(gff_files, INSERTION_SET, DELETION_READ_RATIOS):
    print("running notebook for " + gff_file)
    pm.execute_notebook('subclonal analysis.ipynb', "%s_notebook.ipynb" % insertion_set.lower().replace(' ', '_'), parameters=dict(ALL_INSERTIONS=gff_file, LABEL=insertion_set, DELETION_READ_RATIO=deletion_read_ratio))
