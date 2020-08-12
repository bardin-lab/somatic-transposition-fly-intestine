#!/usr/bin/env bash
# Download all GFF tracks
CURRENT_DIR=`dirname "$0"`
grep gff "$CURRENT_DIR"/../source_files/modencode/modencode_chromatin_download_tracks.tsv|cut -f 2| wget --directory-prefix "$CURRENT_DIR"/../build/modencode  -i -
