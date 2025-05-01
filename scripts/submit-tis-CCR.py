#!/usr/bin/env python
# coding: utf-8

import requests
import json
import glob
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Submit vcfs for prephasing and imputation')
    parser.add_argument('--jobname', required=True, nargs='+', help="""VCF files""")
    parser.add_argument('--vcf', required=True, nargs='+', help="""VCF files""")
    parser.add_argument('--url', default="https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2", help="""imputation server url""")
    parser.add_argument('--post', default="/jobs/submit/imputationserver", help="""imputation server POST""")
    parser.add_argument('--refpanel', default="apps@1000g-phase-3-v5", help="""reference panel. Options: hrc-r1.1, 1000g-phase-3-v5, genome-asia-panel, 1000g-phase-1, cappa, hapmap-2. [Default = 1000g-phase-3-v5]""")
    parser.add_argument('--population', default="all", help="""Population for AF computation. If there are <100 samples, MIS does not drop variants due to AF differences""")
    parser.add_argument('--build', default="hg19", help="""Build. Options: hg19, hg38. [Default = hg19]""")
    parser.add_argument('--r2Filter', default="0", help="""r-squared threshold for imputed variants to return""")
    parser.add_argument('--token', required=True, help="""API token""")
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = getOpts()

    # imputation server url
    url = args.url
    post = args.post
    token = args.token
    vcf_files = args.vcf

    # add token to header (see Authentication)
    headers = {'X-Auth-Token' : token }
    data = {
        'job-name': args.jobname,
        'refpanel': args.refpanel,
        'population': args.population,
        'build': args.build,
        'r2Filter': args.r2Filter,
        'meta': 'yes',
    }

    # submit new job
    files = [('files', open(vcf, 'rb')) for vcf in vcf_files]
    r = requests.post(url + post, files=files, data=data, headers=headers)

    if r.status_code != 200:
        print(r.json()['message'])
        raise Exception('POST {} {}'.format(post, r.status_code))
    else:
        print(r.json()['message'])
        print(r.json()['id'])
