#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import math
import sys
import os
import glob
import pandas_extra
import argparse
from shutil import copyfile

def getOpts():
    parser = argparse.ArgumentParser(description='Small genotyping batch has many monomorphic alleles. Check illumina alleles, account for strand switches and just add a dummy sample so these variants are not discarded.')
    parser.add_argument('--ped', required=True, help="""plink ped file""")
    parser.add_argument('--map', required=True, help="""plink map file""")
    parser.add_argument('--annot', required=True, help="""Illumina annot file""")
    parser.add_argument('--prefix', required=True, help="""Output prefix. Will generate prefix.ped and prefix.map """)
    args = parser.parse_args()
    return args

complement = {'A': 'T',
             'T': 'A',
             'C': 'G',
             'G': 'C'}

def switch(l):
    out = [complement[i] for i in l]
    return out

def compare_alleles(x):
    overlap = len([a for a in x['obs_alleles'] if a in x['exp_alleles']])
    if x['exp_alleles'] == ['9', '9']:
        #not in annotation file
        out = x['obs_alleles']
    elif x['obs_alleles'] == ['0', '0']:
        # homozygous indel - just take the expected
        out = x['exp_alleles']
    elif overlap == 0:
        # alleles switched
        try:
            out = switch(x['exp_alleles'])
        except:
            print(f"check this: {x['obs_alleles']}, {x['exp_alleles']}")
    elif overlap == 2:
        # alleles ok
        out = x['exp_alleles']
    else:
        # some weird issue
        print(f"weird case, check: {x['obs_alleles']}, {x['exp_alleles']}")
        out = x['obs_alleles']
    return out


if __name__ == '__main__':

    args = getOpts()


    """Fix the directory where all figures will be saved"""

    fdir = "."
    pe = pandas_extra.ExtraFunctions(fdir)

    ### get sample count and allele info from ped file
    # ped = "data/6719-NM/Prjt_408_Parker_20220825_6719-NM/Reports/PLINK_290822_0816/Reports.ped"
    # map = "data/6719-NM/Prjt_408_Parker_20220825_6719-NM/Reports/PLINK_290822_0816/Reports.map"
    # annot = "resources/illumina_files/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt"
    # prefix = "results/plink/raw"
    ped = args.ped
    map = args.map
    annot = args.annot
    prefix = args.prefix

    nsamples = len(pandas.read_csv(ped, sep='\t', header=None, usecols=[0,1]).index)
    #nsamples = len(pandas.read_csv(args.ped, sep='\t', header=None, usecols=[0,1]).index)
    m = pandas.read_csv(map, sep='\t', header=None, dtype={0: str})
    #m = pandas.read_csv(args.map, sep='\t', header=None, dtype={0: str})
    print(f"map file had N = {len(m.index)} variants")

    # read first individual from the ped file
    #with open(args.ped, "r") as f:
    with open(ped, "r") as f:
        l1 = f.readline().split("\t")

    # if there are strand switches, adjust alleles accordingly
    ind = 0
    ilist = []
    for i in list(range(6, len(l1))):
        if i % 2 == 0 :
            ilist.append([l1[i], l1[i+1].rstrip()])
        else:
            ind += 1

    m['obs_alleles'] = ilist

    # illumina allele annots
    an = pandas.read_csv(annot, sep='\t', usecols = ['Name','Chr','Alleles'], dtype={'Chr': str})
    #an = pandas.read_csv(args.annot, sep='\t', usecols = ['Name','Chr','Alleles'], dtype={'Chr': str})

    m1 = pandas.merge(m, an, how="left", left_on=1, right_on="Name")[['Name','Chr','Alleles','obs_alleles']]
    m1['Alleles'] = m1['Alleles'].str.replace("[", "").str.replace("]", "")
    #print(f"after merging with annot file on variant name, N = {len(m1.index)} variants")

    # def unique(list1):
    #     unique_list = []
    #     for x in list1:
    #         if x not in unique_list:
    #             unique_list.append(x)
    #     for x in unique_list:
    #         print(x)

    m1['Alleles'] = m1['Alleles'].fillna('9/9')
    m1["exp_alleles"] = m1['Alleles'].str.split("/")
    m1['out_alleles'] = m1.apply(lambda x: compare_alleles(x), axis=1)

    m1['out_alleles_str'] = m1['out_alleles'].map(lambda x : "/".join(x))
    #m2 = pe.explode_df(m1, 'out_alleles_str', "/")
    m2 = m1.explode('out_alleles')

    assert len(m2.index) == 2 * len(m1.index)

    out = "\t".join([f'{nsamples + 1}', "dummy", "0", "0", "0", "0"] + m2['out_alleles'].tolist())

    newped = f"{prefix}.ped"
    newmap = f"{prefix}.map"
    copyfile(ped, newped)
    copyfile(map, newmap)

    with open(newped, "a") as f:
        f.write(out)
