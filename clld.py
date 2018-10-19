#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 00:16:34 2018

@author: chenrui
"""

import vcf
from collections import defaultdict
import pandas as pd
import sys


prefix = sys.argv[1]
vcf_file = sys.argv[2]
bed_file = sys.argv[3]
parents = sys.argv[4].split(',')
offspring = sys.argv[5].split(',')

if len(parents) == 2:
    parent1 = parents[0]
    parent2 = parents[1]
else:
    print('Error: the number of parents is wrong, it should be 2. ')

if len(offspring) < 1:
    print('Error: the number of offsprings is less than 1')

vcf_reader = vcf.Reader(open(vcf_file,'r'))
region = pd.read_table(bed_file,header=None,names=['chr','start','end'])
result = {}

def is_large_deletion(record,sample):
    call = record.genotype(sample)
    flag = False
    if 'DEL' in record.ID and call['FT'] == 'PASS':
        gt = call['GT'].split('/')
        for i in gt:
            if i == '1':
                flag = True
    return flag

def de_novo_deletion(record,proband,parent1,parent2):
    is_de_novo = True
    call_proband = record.genotype(proband)
    call_parent1 = record.genotype(parent1)
    call_parent2 = record.genotype(parent2)
    if is_large_deletion(record,proband) == True and call_parent1['FT'] == 'PASS' and call_parent2['FT'] == 'PASS':    
        gt_proband = call_proband['GT']
        gt_parent1 = call_parent1['GT'].split('/')
        gt_parent2 = call_parent2['GT'].split('/')
        for i in gt_parent1:
            for j in gt_parent2:
                temp = [i, j]
                temp.sort()
                GT = temp[0]+'/'+temp[1]
                if GT == gt_proband:
                    is_de_novo = False
                    break
    else:
        is_de_novo = False
    return is_de_novo

for record in vcf_reader:
    select = region[(record.CHROM == region['chr']) & (record.start <= region['end']) & (record.end >= region['start'])]
    if len(select) > 0:
        for proband in offspring:
            if de_novo_deletion(record, proband, parent1,parent2):
                result['ID'] = record.ID
                result['chromosome'] = record.CHROM
                result['start'] = record.start
                result['end'] = record.end
                result['length'] = record.end - record.POS
                result['sample'] = proband
                result['target'] = select
                result['genotype'] = record.genotype(proband)['GT']               
result = pd.DataFrame(result)
result.to_csv('{0}.csv'.format(prefix),index=False)
