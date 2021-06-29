#!/usr/bin/env python
import sys
import pickle
import urllib.request, urllib.parse, urllib.error
import gzip
from cyvcf2 import VCF, Writer
import os.path


if __name__ == "__main__":
    # Download the current set of gene IDs
    if len(sys.argv) == 1 and not os.path.isfile("gene.pkl"):
        urllib.request.urlretrieve('ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.current.geneIDs.txt.gz', 'gene_ids.txt.gz')
        with gzip.open('gene_ids.txt.gz', 'rb') as f:
            gene_set = [x.decode('utf8').split(',') for x in f.read().splitlines()]
            gene_set = {x[1]:(x[2] or x[3]) for x in gene_set if x[4] == 'Live'}
            pickle.dump(gene_set, open("gene.pkl", 'wb'))
    else:
        v = VCF("-")
        print((v.raw_header.strip()))
        gene_set = pickle.load(open("gene.pkl", 'rb'))
        for line in v:
            ANN = line.INFO.get("ANN")
            BCSQ = line.INFO.get("BCSQ")
            if ANN:
                ann = [x.split("|") for x in ANN.split(",")]
                for x in ann:
                    if x[3] in list(gene_set.keys()):
                        x[3] = gene_set[x[3]]
                ann = ','.join(["|".join(x) for x in ann])
                line.INFO['ANN'] = ann
            if BCSQ:
                bcsq = [x.split("|") for x in BCSQ.split(",")]
                for x in bcsq:
                    if len(x) > 1:
                        if x[1] in list(gene_set.keys()):
                            x[1] = gene_set[x[1]]
                bcsq = ','.join(["|".join(x) for x in bcsq])
                line.INFO['BCSQ'] = bcsq
            print((str(line).strip()))
                    
