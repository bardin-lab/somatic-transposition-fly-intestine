modencode_download: 
	scripts/download_modenencode_data.sh

modencode_liftover: modencode_download
	python scripts/sanitize_and_liftover.py --chainfile scripts/dm3ToDm6.over.chain.gz build/modencode/*.gz build/modencode/dm6 --contig_whitelist X --contig_whitelist 2L --contig_whitelist 2R --contig_whitelist 3L --contig_whitelist 3R --contig_whitelist Y --contig_whitelist 4

modencode: modencode_download modencode_liftover
