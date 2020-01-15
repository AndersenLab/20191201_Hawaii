import os
import json
import gzip
from subprocess import Popen, PIPE
from os.path import dirname
from collections import OrderedDict


def get_git_dir():
    path = Popen(['git',
                  'rev-parse',
                  '--git-dir'], stdout=PIPE).communicate()[0]
    return str(dirname(path), encoding='UTF-8')

def get_wd():
    cwd = os.getcwd()
    return cwd

#wd = get_git_dir()
wd = get_wd()
blast_json = wd + '/data/sanger/blast_results1.json.gz'
blast_results = json.load(gzip.open(blast_json))['BlastOutput2']

header_out = False

with open(wd + '/data/sanger/blast_results1.tsv', 'w') as f:
    for record in blast_results:
        output = OrderedDict()
        record = record['report']['results']['search']
        output.update({'seq': record['query_title']})
        print(record['query_title'])
        if 'control' in record['query_title']:
            record['query_title'] = record['query_title'].replace("____control_","control-")
        elif record['query_title'].count("_") == 4:
            s_plate, primer, tmp1, tmp2, tmp3 = record['query_title'].split("_")
            well, plate, sanger_id, rec = "NA", "NA", "NA", "NA"
        elif record['query_title'].count("_") == 5:
            s_plate, primer, well, tmp1, tmp2, tmp3 = record['query_title'].split("_")
            plate, sanger_id, rec = "NA", "NA", "NA"
        else:
            try:
                s_plate, primer, well, sanger_id = record['query_title'].split("_")
            except:
                s_plate, primer = record['query_title'].split("_")
                well = "NA"
                sanger_id = "NA"
            plate = "NA"
            rec = "NA"
        folder = record['query_title'].split("_")[-1]
        output.update({'s_plate': s_plate,
                    'primer': primer,
                    'well': well,
                    'plate': plate,
                    'sanger_id': sanger_id,
                    'rec': rec,
                    'folder': folder})
        for hit in record['hits'][:3]:
            output.update(hit['description'][0])
            alignment = hit['hsps'][0]
            output.update(alignment)

            if header_out is False:
                f.write('\t'.join(output.keys()) + "\n")
                header_out = True
            f.write('\t'.join(map(str, output.values())) + "\n")