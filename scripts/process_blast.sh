#!/usr/bin/bash

# Output available S-plates
SPLATES=`cat data/sanger/blast_results.tsv | cut -f 2 | sort | uniq | egrep -v 's_plate|NA' | awk 'NR == 1 { print "CONTROL" } { print }'`

# Run them!
for s in ${SPLATES}; do
    echo "RUNNING ${s}"
    Rscript -e "rmarkdown::render('02-blast.Rmd', output_file = '${s}.html', output_dir ='data/sanger/reports/')" ${s}
done;