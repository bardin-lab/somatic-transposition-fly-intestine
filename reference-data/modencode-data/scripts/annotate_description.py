import os
import pandas as pd


DESC = '../source_files/modencode/modencode_chromatin_track_descriptions.tsv'

description = pd.read_csv(DESC, sep='\t', usecols=list(range(8)))
# Create smaller categories
description['Tissue'] = None
description.loc[description.Dataset.str.contains('S2'), 'Tissue'] = 'S2'
description.loc[description.Dataset.str.contains('Kc167'), 'Tissue'] = 'Kc'
description.loc[description.Dataset.str.contains('ML-DmBG3'), 'Tissue'] = 'ML-DmBG3'
description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Embryo'), 'Tissue'] = 'Embryo'
description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Larvae'), 'Tissue'] = 'Larvae'
description.loc[description.Tissue.isnull() & description.Dataset.str.contains('L3'), 'Tissue'] = 'Larvae'
description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Adult'), 'Tissue'] = 'Adult'
description.loc[description.Tissue.isnull() & description.Dataset.str.contains('upae'), 'Tissue'] = 'Pupae'

original_filenames = [os.path.basename(p) for p in os.listdir('../build/modencode/dm6/gff') if p.endswith('gff.gz')]

