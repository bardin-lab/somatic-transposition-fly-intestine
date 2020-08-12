import gzip
import logging
import os
import shlex
import subprocess
from tempfile import TemporaryDirectory
import time

import click
import pysam

logger = logging.getLogger("GFF liftover")

GFF = {'GFF', 'gff', 'GFF3', 'gff3', 'GFF.gz', 'gff.gz', 'GFF3.gz', 'gff3.gz'}
BED = {'BED', 'bed', 'BED.gz', 'bed.gz'}


def get_track_type(infile):
    if any((True for ext in GFF if infile.endswith(ext))):
        return 'gff'
    if any((True for ext in BED if infile.endswith(ext))):
        return 'bed'


def prefix_chromosome(infile, outfile, drop_contig, track_type):
    """
    Prefix first column of infile with 'chr'

    The liftover files from UCSC require chromosome names to start with 'chr'.
    Additionally whitespaces are not tolerated, so we exchange white spaces with underscores.
    """
    with gzip.open(infile, 'rt') as fh_in, open(outfile, 'w') as fh_out:
        i = 0
        for i, line in enumerate(fh_in):
            if line.startswith(('#', '"', 'track', 'browser', 'chr\t', 'CHROM', 'chrom')):
                # This is the header
                continue
            else:
                fields = line.strip().split('\t')
                if len(fields) == 1:
                    fields = fields[0].split(' ')
                if not fields[0] in drop_contig:
                    for n, field in enumerate(fields):
                        # Spaces inside fields make liftOver go boom
                        if not field:
                            field = '.'
                        fields[n] = field.replace(' ', '_')
                    if track_type == 'bed' or len(fields) > 6:  # GFF should have 9 fields, but some gff files are incompatible
                        if not fields[0].startswith('chr'):
                            fields[0] = "chr%s" % fields[0]
                        if track_type == 'gff':
                            while len(fields) < 9:
                                # Insert '.' as missing field values
                                fields.insert(-1, '.')
                        elif track_type == 'bed':
                            fields = fields[:4]
                        fh_out.write("%s\n" % "\t".join(fields))
        logger.info("'%s' has %d records", infile, i)
    return i


def drop_prefix(infile, outfile, drop_contig, contig_whitelist):
    """Remove 'chr' from first column."""
    with open(infile) as fh_in, open(outfile, 'w') as fh_out:
        i = 0
        for i, line in enumerate(fh_in):
            if line.startswith('chr'):
                line = line[3:]
                fields = line.split('\t')
                if fields[0] in drop_contig:
                    continue
                if contig_whitelist and fields[0] not in contig_whitelist:
                    continue
            fh_out.write(line)
        logger.info("After lifting over to new coordinates '%s' has %d records", infile, i)
    return i


def sort_file(infile, outfile):
    """Sort files using bedtools sort."""
    tmp_out = "%s.tmp" % outfile
    with open(tmp_out, 'w') as sorted_out:
        subprocess.check_call(['bedtools', 'sort', '-i', infile], stdout=sorted_out)
    if os.path.exists(outfile):
        os.remove(outfile)
    os.rename(tmp_out, outfile)


def split_depleted_and_enriched(infile, track_type):
    enriched = False
    depleted = False
    ext = track_type
    enriched_out = "%s.enriched.%s" % (os.path.abspath(infile), ext)
    depleted_out = "%s.depleted.%s" % (os.path.abspath(infile), ext)
    out_files = [infile]
    with open(infile) as fh_in, open(enriched_out, 'w') as enriched_fh, open(depleted_out, 'w') as depleted_fh:
        for line in fh_in:
            if not line.startswith('#'):
                fields = line.split('\t')
                if 'enrich' in fields[1]:
                    enriched = True
                    enriched_fh.write("\t".join(fields))
                elif 'deplet' in fields[1]:
                    depleted = True
                    depleted_fh.write("\t".join(fields))
    if enriched:
        out_files.append(enriched_out)
    else:
        os.remove(enriched_out)
    if depleted:
        out_files.append(depleted_out)
    else:
        os.remove(depleted_out)
    if enriched or depleted:
        out_files = out_files[1:]
        os.remove(infile)
    return out_files


def index_file(infile, track_type):
    """Compress and Index GFF files."""
    pysam.tabix_index(infile, preset=track_type)
    return "{infile}.gz".format(infile=infile)


def liftover(chainfile, infile, outfile, tmp_dir, track_type):
    """Lift over input file."""
    name = os.path.basename(infile)
    unmapped = os.path.join(tmp_dir, "%s.unmapped.gff" % name)
    track_switch = '-gff' if track_type == 'gff' else ''
    cmd = shlex.split("liftOver {track_switch} {infile} {chainfile} {outfile} {unmapped}".format(infile=infile,
                                                                                                 track_switch=track_switch,
                                                                                                 chainfile=chainfile,
                                                                                                 outfile=outfile,
                                                                                                 unmapped=unmapped))
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        logger.error(e.output)


def mass_liftover(infiles, output_directory, chainfile, drop_contig, contig_whitelist, write_unlifted=False):
    """Lift over a number of GFF files at infiles and place them into `output_directory`."""
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    with TemporaryDirectory() as tmp_dir:
        for infile in infiles:
            track_type = get_track_type(infile)
            base_name = os.path.basename(infile)
            if base_name.endswith('.gz'):
                base_name = base_name[:-3]
            prefixed_file = os.path.join(tmp_dir, "%s.prefixed.%s" % (base_name, track_type))
            unprefixed_file = os.path.join(tmp_dir, "%s.unprefixed.%s" % (base_name, track_type))
            liftover_file = os.path.join(tmp_dir, "%s.liftover.%s" % (base_name, track_type))
            no_liftover_file = os.path.join(tmp_dir, "%s.no_liftover.%s" % (base_name, track_type))
            output_no_liftover_file = os.path.join(output_directory, "no_liftover.%s" % base_name)
            output_file = os.path.join(output_directory, base_name)
            input_record_number = prefix_chromosome(infile=infile, outfile=prefixed_file, drop_contig=drop_contig, track_type=track_type)
            liftover(chainfile=chainfile, infile=prefixed_file, outfile=liftover_file, tmp_dir=tmp_dir, track_type=track_type)
            output_record_number = drop_prefix(liftover_file, outfile=unprefixed_file, drop_contig=drop_contig, contig_whitelist=contig_whitelist)
            logger.info("Lost %d lines during liftover of %s", input_record_number - output_record_number, infile)
            if write_unlifted:
                output_record_number = drop_prefix(prefixed_file, outfile=no_liftover_file, drop_contig=drop_contig, contig_whitelist=contig_whitelist)
            try:
                sort_file(unprefixed_file, output_file)
                if write_unlifted:
                    sort_file(no_liftover_file, output_no_liftover_file)
            except Exception as e:
                logger.error(e)
                print(unprefixed_file, output_file)
                continue
            output_files = split_depleted_and_enriched(output_file, track_type=track_type)
            if write_unlifted:
                no_liftover_output_files = split_depleted_and_enriched(output_no_liftover_file, track_type=track_type)
            for output_file in output_files:
                index_file(output_file, track_type)
            if write_unlifted:
                for output_file in no_liftover_output_files:
                    index_file(output_file, track_type)


@click.command("Lift over and sanitize many GFF files")
@click.argument('infiles', nargs=-1, type=click.Path(exists=True), required=True)
@click.argument('output_directory', nargs=1, required=True)
@click.option('--chainfile', type=click.Path(exists=True), required=True, help="Chain file for conversion.")
@click.option('--drop_contig', multiple=True)
@click.option('--contig_whitelist', multiple=True)
@click.option('--write_unlifted/--no_write_unlifted', help="Output sanitized but unlifted files", default=False)
@click.option('--logfile', default=None, help="Log messages to this file")
@click.option('--loglevel', default='INFO')
def cli(infiles, output_directory, chainfile, loglevel, drop_contig, contig_whitelist, logfile=None, write_unlifted=False):
    logging.basicConfig(format='%(asctime)s %(name)s %(levelname)s - %(message)s',
                        filename=logfile,
                        level=loglevel)
    mass_liftover(infiles=infiles, output_directory=output_directory, chainfile=chainfile, drop_contig=drop_contig, contig_whitelist=contig_whitelist, write_unlifted=write_unlifted)


if __name__ == '__main__':
    cli()
