#!/usr/bin/env bash
set -e
set -o pipefail

GFF_DIR="${GFF_DIR:-../build/modencode/dm6/gff/}"
BED_DIR="${BED_DIR:-../build/modencode/dm6/bed/}"

mkdir "$BED_DIR"
for i in $(ls "$GFF_DIR"/*.gz)
do 
    zcat < "$i" | sed -E 's/ID=(region_.[0-9]*).*/ID=\1/g' |bedtools sort -i -|convert2bed --input=gff - |cut -f1,2,3 > "$i".bed && mv "$i".bed "$BED_DIR"
done
