#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author: Daniel E. Cook

This script summarizes variation for samples.

"""

import sys
import numpy as np
import json
from cyvcf2 import VCF
from collections import defaultdict


ANN_fields = ["allele",
              "effect",
              "impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "exon_intron_rank",
              "nt_change",
              "aa_change",
              "cDNA_position/cDNA_len",
              "protein_position",
              "distance_to_feature",
              "error"]

BCSQ_fields = ['consequence',
               'gene_name',
               'ensembl_transcript_id',
               'coding_strand',
               'amino_acid_position',
               'variants']

GENOTYPES = {0: 'homozygous_ref',
             1: 'heterozygous',
             2: 'homozygous_alt',
             3: 'missing'}


args = sys.argv
if len(args) > 1:
    args[1]

if not sys.stdin.isatty():
    vcf_in = "-"
else:
    args = sys.argv
    if len(args) > 1:
        vcf_in = args[1]

def int_to_bool_list(num):
    return [bool(num & (1<<n)) for n in range(8)]

vcf = VCF(vcf_in, gts012=False)
tree = lambda: defaultdict(tree)
results = tree()


for line in vcf:
    genotype_counts = dict(zip(*np.unique(line.gt_types, return_counts = True)))
    # For ANN annotations
    ANN_annotations = []
    allele_set = [line.REF] + line.ALT
    if 'ANN' in dict(line.INFO).keys():
        ANN_annotations = [dict(zip(ANN_fields, x.split("|"))) for x in line.INFO['ANN'].split(",")]
    for sample, gt, gt_num in zip(vcf.samples, line.gt_types, line.genotypes):
        gt_name = GENOTYPES[gt]
        if type(results[sample]["gt_count"][gt_name][str(genotype_counts[gt])]) == int:
            results[sample]["gt_count"][gt_name][str(genotype_counts[gt])] += 1
        else:
            results[sample]["gt_count"][gt_name][str(genotype_counts[gt])] = 1
    
        for anno in ANN_annotations:
            try:
                if allele_set.index(anno.get('allele')) in gt_num:
                    ANN = results[sample]['ANN']
                    gt_alleles = [allele_set[x] for x in gt_num[:-1]]
                    # Initialize counters
                    if 'impact' not in ANN.keys():
                        ANN['impact'], ANN['effect'], ANN['transcript_biotype'] = defaultdict(int), defaultdict(int), defaultdict(int)
                        ANN['HIGH_impact_genes'] = []
                    if anno['allele'] in gt_alleles:
                        ANN['impact'][anno['impact']] += 1
                        for eff in anno['effect'].split("&"):
                            ANN['effect'][eff] += 1
                        ANN['transcript_biotype'][anno['transcript_biotype']] += 1
                        if anno['impact'] == 'HIGH':
                            anno['CHROM'] = line.CHROM
                            anno['POS'] = line.POS
                            anno['REF'] = line.REF
                            anno['ALT'] = ','.join(line.ALT)
                            ANN['HIGH_impact_genes'].append(anno)
            except ValueError:
                pass

print(json.dumps(results))


# For BCSQ annotations
#if 'BCSQ' in dict(line.INFO).keys():
#    annotations = [dict(zip(BCSQ_fields, x.split("|"))) for x in line.INFO['BCSQ'].split(",")]
#    haplotypes = line.format('BCSQ')
#    print(np.unique(haplotypes))
#    if len(np.unique(haplotypes)) > 2 and any([a['consequence'].startswith("mis") for a in annotations]):
#        break







