# 20191201_Hawaii
Hawaiian nematode sampling conducted by Clay Dilks, Gaotian Zhang, Erik Andersen, and Chris Andersen.

# Using automated blast scripts
1. Move all sequencing files into the Sanger/raw folder in repo. They can be in different subfolders. Can include .ab1, .seq, whatever, best to include all.
2. Run sanger.py script to trim sequences and save in a single .fasta file
3. Load the single fast file in the blastn web interface and run
4. On the output page just to the right of the the RID line (line below “Job Title”) click on “Download All” drop-down.
    * Select Single file JSON
    * Move the downloaded file to the data/sanger folder of the collection repo and rename it as blast_results(unique number).json
    * Use gzip to compress command from the command line to compress the file command: gzip your/JSON/file.json.gz
    * Run the process_blast.py script to create the blast_results(unique number).tsv file.
    * Note, you can run multiple times as seq files accumulate in the Sanger folder. Once all are present you can run one final time to get a full blast summary.
