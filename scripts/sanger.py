#!/usr/bin/env python
"""
    Trims Sanger sequences with low quality.
"""
import glob
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import json
from os.path import basename, dirname, join
from subprocess import Popen, PIPE
from Bio.SeqIO import FastaIO

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_git_dir():
    path = Popen(['git', 'rev-parse', '--git-dir'], stdout=PIPE).communicate()[0]
    return str(dirname(path), encoding='UTF-8')

def get_wd():
    cwd = os.getcwd()
    return cwd

def trim_sanger(seq, minimum = 10):
    """
        Trims the start and end of a sanger sequence
        until quality surpasses the minimum.
    """
    qualities = seq.letter_annotations['phred_quality']
    qual_set = [x <= minimum for x in qualities]
    if False in qual_set:
        start = qual_set.index(False)
        seq_length = len(qualities)
        end = seq_length - [x <= minimum for x in qualities[::-1]].index(False)
        if len(seq[start:end]) > 100:
            return seq[start:end], start, (seq_length - end)
    return None, "-", "-"

#wd = get_git_dir()
wd = get_wd()

for i in [0, 10, 20, 30]:

    # Open trimming file
    sanger_trimming = open(wd + '/data/sanger/sanger_trimming_{}.tsv'.format(i), 'w')
    sanger_trimming.write('seq\ttrim_left\ttrim_right\tprimer_type\n')

    # Open fasta sequences file
    fasta_seqs = open(wd + '/data/sanger/sanger_{}.fasta'.format(i), 'w')
    fasta_out = FastaIO.FastaWriter(fasta_seqs, wrap=None)
    seq_set = []
    for s in glob.glob(wd + '/data/sanger/raw/*/*ab1')[11:]:
        # Only look in confirmationseq
        primer_type = dirname(s).split("/")[-1]
        folder = basename(dirname(s))
        # Trim sanger sequence
        seq = SeqIO.read(s, 'abi')
        sname = basename(s).replace('.ab1', '') + "_" + primer_type
        s_plate, primer, well = sname.split("_")[0:3]
        seq.name = sname
        seq.id = sname
        sanger, trim_left, trim_right = trim_sanger(seq, i)

        sanger_trimming.write('{}\t{}\t{}\t{}\n'.format(sname, trim_left, trim_right, primer_type))
        if sanger:
            seq_set.append(sanger)
    fasta_out.write_file(seq_set)
