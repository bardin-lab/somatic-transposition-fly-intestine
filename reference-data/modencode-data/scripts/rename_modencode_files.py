import click
import os
import pandas as pd


def annotate_description(description, interval_file_path, annotated_output):
    description = pd.read_csv(description, sep='\t', usecols=list(range(8)))
    # Create smaller categories
    description['Tissue'] = None
    description.loc[description.Dataset.str.contains('S2'), 'Tissue'] = 'S2'
    description.loc[description.Dataset.str.contains('Kc167'), 'Tissue'] = 'Kc'
    description.loc[description.Dataset.str.contains('ML-DmBG3'), 'Tissue'] = 'ML-DmBG3'
    description.loc[description.Dataset.str.contains('CME-W1-Cl.8+'), 'Tissue'] = 'CME-W1-Cl.8+'
    description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Embryo'), 'Tissue'] = 'Embryo'
    description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Larvae'), 'Tissue'] = 'Larvae'
    description.loc[description.Tissue.isnull() & description.Dataset.str.contains('L3'), 'Tissue'] = 'Larvae'
    description.loc[description.Tissue.isnull() & description.Dataset.str.contains('Adult'), 'Tissue'] = 'Adult'
    description.loc[description.Tissue.isnull() & description.Dataset.str.contains('upae'), 'Tissue'] = 'Pupae'
    id_file = {}
    for p in os.listdir(interval_file_path):
        if p.endswith('gff.gz') or p.endswith('gff') or p.endswith('gff3.gz') or p.endswith('gff3') or p.endswith('.bed.gz') or p.endswith('.bed'):
            try:
                ID = int(p.split('_')[0])
            except ValueError:
                print("Couldn't read file %s" % p)
                continue
            id_file[p] = ID
    original_filenames = pd.DataFrame.from_dict(id_file, orient='index')
    original_filenames.columns = ['ID']
    original_filenames['original_filename'] = original_filenames.index
    description = original_filenames.merge(description, on='ID')
    description['region_type'] = 'enriched'
    description.loc[description.original_filename.str.contains('depleted'), 'region_type'] = 'depleted'
    description['basename'] = description['original_filename'].str.replace('.enriched.*', '')
    description['basename'] = description['basename'].str.replace('.depleted.*', '')
    description['genotype'] = 'wt'
    description.loc[description['Target Element'] == 'Histone ModificationandNon TF Chromatin binding factor', 'Target Element'] = 'Non TF Chromatin binding factor'
    description.loc[description['Target Element'] == 'Transcriptional FactorandNon TF Chromatin binding factor', 'Target Element'] = 'Non TF Chromatin binding factor'  # TRL
    description.loc[description['Target Element'] == 'Non TF Chromatin binding factorandTranscriptional Factor', 'Target Element'] = 'Non TF Chromatin binding factor'  # PolII
    description.loc[description['Target Element'] == 'DNA ReplicationandNon TF Chromatin binding factor', 'Target Element'] = 'Non TF Chromatin binding factor'  # PolII
    description.loc[description.Dataset.str.contains('utant'), 'genotype'] = description.Dataset.str.split(';').str.get(1)

    # Figure out the replicate number
    grouped = description.groupby('ID')
    fname_rep = {}
    for name, group in grouped:
        group_fnames = group['basename'].unique()
        for i, fname in enumerate(group_fnames):
            fname_rep[fname] = i
    rep_df = pd.DataFrame.from_dict(fname_rep, orient='index')
    rep_df.columns = ['replicate']
    description = description.merge(rep_df, left_on='basename', right_index=True)

    description['new_filename'] = description['Assay Factor'] + '_' + description['genotype'] + '_' + description['Tissue'] + '_' + description['Target Element'] + '_' + description['replicate'].astype(str) + '_'  + description['region_type'] + '_' + description['ID'].astype(str)
    description['new_filename'] = description['new_filename'].str.replace('/', '_').str.replace(' ', '-')
    description.to_csv(annotated_output, sep='\t')
    return description


def rename_files(description, directory):
    files = os.listdir(directory)
    basename_new_filename_lookup = dict(zip(description['original_filename'], description['new_filename']))
    for f in files:
        if f.endswith('.gff'):
            suffix = '.gff'
        elif f.endswith('.gff.gz'):
            suffix = '.gff.gz'
        elif f.endswith('.gff3'):
            suffix = '.gff'
        elif f.endswith('.gff3.gz'):
            suffix = '.gff.gz'
        elif f.endswith('.bed'):
            suffix = '.bed'
        elif f.endswith('.bed.gz'):
            suffix = '.bed.gz'
        else:
            print("%s is not supported, should be a bed/gff/bed.gz/gff.gz file." % f)
            continue
        for original_filename, new_filename in basename_new_filename_lookup.items():
            if f == original_filename:
                new_path = os.path.join(directory, "%s%s" % (new_filename, suffix))
                old_path = os.path.join(directory, f)
                os.rename(src=old_path, dst=new_path)
                break
        else:
            print("Could not find '%s'in annotation file, skipping rename" % f)


@click.command("Annotate and move downloaded modencode files to human-readable names")
@click.option('--description', help="Path to modencode_chromatin_track_descriptions.tsv", required=True, type=click.Path(exists=True))
@click.option('--annotated_output', help="Path to write annotated table file", required=True)
@click.option('--interval_directory', help="Path to directory containing modencode gff/bed files", required=True, type=click.Path(exists=True))
def main(description, interval_directory, annotated_output):
    description = annotate_description(description=description, interval_file_path=interval_directory, annotated_output=annotated_output)
    rename_files(description, directory=interval_directory)

if __name__ == '__main__':
    main()
