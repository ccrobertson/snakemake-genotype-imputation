#!/usr/bin/env python3


configfile: "config.yaml"
configfile: "token.yaml"
CHROMOSOMES = list(range(1,23))
CHROMOSOMES.append("X")

def get_1000g_file(chr):
    return config["KGref_hg19"]["chr"+chr]

def get_refpanel_info(refpanel):
    if refpanel == 'topmed':
        out = {
            "url": "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2",
            "post": "/jobs/submit/imputationserver",
            "refpanel": "apps@topmed-r3",
            "token": config["tis_token"]
            }
    elif refpanel == '1000g':
        out = {
            "url": "https://imputationserver.sph.umich.edu/api/v2",
            "post": "/jobs/submit/minimac4",
            "refpanel": "1000g-phase-3-v5",
            "token": config["mis_token"]
            }
    elif refpanel == "hla_4d":
        out = {
            "url": "https://imputationserver.sph.umich.edu/api/v2",
            "post": "/jobs/submit/imputationserver-hla",
            "refpanel": "multiethnic-hla-panel-4digit",
            "token": config["mis_token"]
            }
    elif refpanel == "hla_g":
        out = {
            "url": "https://imputationserver.sph.umich.edu/api/v2",
            "post": "/jobs/submit/imputationserver-hla",
            "refpanel": "multiethnic-hla-panel-Ggroup",
            "token": config["mis_token"]
            }
    else:
        raise ValueError("That reference panel is not supported.")
    return(out)


rule all:
    input:
        #expand("results/imputation_prep/clean-updated-add1000g-chr{chr}.vcf.gz", chr=CHROMOSOMES)
        #"results/imputation_prep/clean-updated-add1000g-chrX.vcf.gz"
        #results/imputation_prep/clean-updated-add1000g-chr{chr}.vcf.gz
        "results/imputation-topmed-r2Filter_0.3/tis.submitted",
        #"results/imputation-1000g-r2Filter/mis.submitted",
        #"results/imputation_hla-hla_4d-r2Filter/mis.submitted",
        #"results/imputation_hla-hla_g-r2Filter/mis.submitted",



rule freq:
    input:
        bim = config["input_plink_prefix"]+".bim",
    output:
        frq = "results/imputation_prep/clean.frq",
        bim = "results/imputation_prep/clean.bim",
    params:
        input_prefix = config["input_plink_prefix"],
        output_prefix = "results/imputation_prep/clean",
    shell:
        """        
        plink1.9 --bfile {params.input_prefix} --freq --make-bed --out {params.output_prefix}
        """

rule align:
    input:
        bim = "results/imputation_prep/clean.bim",
        frq = "results/imputation_prep/clean.frq",
        reference_1000g = "resources/1000GP_Phase3_combined.legend"
    output:
        "results/imputation_prep/clean-updated.bim",
        expand("results/imputation_prep/clean-updated-chr{chr}.vcf", chr=range(1,24))
    params:
        outdir = "results/imputation_prep",
    shell:
        """
        scripts/raynor_imputation_prep/HRC-1000G-check-bim-NoReadKey.pl --bim {input.bim} --frequency {input.frq} --ref {input.reference_1000g} --1000g --output {params.outdir} --verbose
        grep -v '^rm' {params.outdir}/Run-plink.sh | sed 's/--recode vcf/--recode vcf-iid/g'> {params.outdir}/Run-plink_fix.sh
        bash {params.outdir}/Run-plink_fix.sh
        """

rule update_sex_chr:
    input:
        bim = "results/imputation_prep/clean-updated.bim"
    output:
        vcf = "results/imputation_prep/clean-updated-chrX.vcf",
    params:
        input_prefix = "results/imputation_prep/clean-updated",
        output_prefix = "results/imputation_prep/clean-updated-chrX",
    shell:
        """
        plink --bfile {params.input_prefix} --real-ref-alleles --make-bed --chr 23 --output-chr M --out {params.output_prefix}
        plink --bfile {params.output_prefix} --real-ref-alleles --recode vcf-iid --output-chr M --out {params.output_prefix}
        """

rule zip_and_index:
    input:
        vcfs = expand("results/imputation_prep/clean-updated-chr{chr}.vcf", chr=CHROMOSOMES),
    output:
        gz = expand("results/imputation_prep/clean-updated-chr{chr}.vcf.gz", chr=CHROMOSOMES),
        tbi = expand("results/imputation_prep/clean-updated-chr{chr}.vcf.gz.tbi", chr=CHROMOSOMES),
    conda:
        "general"
    shell:
        """
        for i in {input.vcfs}; do bgzip -k $i; tabix ${{i}}.gz ; done;
        """

rule filter_1000g:
    input:
        vcf_array = "results/imputation_prep/clean-updated-chr{chr}.vcf.gz",
        vcf_1000g = lambda wildcards: get_1000g_file(wildcards.chr),
    output:
        vcf = "results/imputation_prep/1000g_bctfools_isec_chr{chr}/0001.vcf.gz",
    params:
        outdir = "results/imputation_prep/1000g_bctfools_isec_chr{chr}",
    conda:
        "general"
    shell:
        """
        #extract records from 1000g vcf that are in array vcf using exact allele match
        bcftools isec -p {params.outdir} -n=2 -c none -w2 {input.vcf_array} {input.vcf_1000g} -Oz
        """


rule merge_with_1000g:
    input:
        vcf_array = "results/imputation_prep/clean-updated-chr{chr}.vcf.gz",
        vcf_1000g = "results/imputation_prep/1000g_bctfools_isec_chr{chr}/0001.vcf.gz",
    output:
        vcf = "results/imputation_prep/clean-updated-add1000g-chr{chr}.vcf.gz",
    conda:
        "general"
    shell:
        """
        bcftools merge --merge none {input.vcf_array} {input.vcf_1000g} -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

### NOTE: it looks like the API doesn't actually support the r2 filter
### Assumes python3
rule impute:
    input:
        vcfs = expand("results/imputation_prep/clean-updated-add1000g-chr{chr}.vcf.gz", chr=CHROMOSOMES)
    output:
        confirmation = "results/imputation-{refpanel}-r2Filter_{r2Filter}/tis.submitted",
    params:
        jobname = config["jobname"],
        url = lambda wildcards: get_refpanel_info(wildcards.refpanel)['url'],
        post = lambda wildcards: get_refpanel_info(wildcards.refpanel)['post'],
        refpanel_keyword = lambda wildcards: get_refpanel_info(wildcards.refpanel)['refpanel'],
        token = lambda wildcards: get_refpanel_info(wildcards.refpanel)['token'],
        vcfstring = ' '.join(expand("results/imputation_prep/clean-updated-add1000g-chr{chr}.vcf.gz", chr=CHROMOSOMES)),
        r2Filter = "{r2Filter}",
    conda:
        "base"
    shell:
        """
        python scripts/submit-tis-CCR.py \
            --jobname {params.jobname} \
            --vcf {params.vcfstring} \
            --url {params.url} \
            --post {params.post} \
            --r2Filter {params.r2Filter} \
            --refpanel {params.refpanel_keyword} \
            --population all \
            --build hg19 \
            --token {params.token} &> {output.confirmation}.log
        echo "FLAG" > {output.confirmation}
        """

# rule impute_hla:
#     input:
#         vcf = "results/imputation_prep/clean-updated-add1000g-chr6.vcf.gz",
#     output:
#         confirmation = "results/imputation_hla-{refpanel}-r2Filter{r2Filter}/mis.submitted",
#     params:
#         url = lambda wildcards: get_refpanel_info(wildcards.refpanel)['url'],
#         post = lambda wildcards: get_refpanel_info(wildcards.refpanel)['post'],
#         refpanel_keyword = lambda wildcards: get_refpanel_info(wildcards.refpanel)['refpanel'],
#         token = lambda wildcards: get_refpanel_info(wildcards.refpanel)['token'],
#         r2Filter = "{r2Filter}",
#     conda:
#         "base"
#     shell:
#         """
#         python scripts/submit-mis-CCR.py \
#             --vcf {input.vcf} \
#             --url {params.url} \
#             --post {params.post} \
#             --r2Filter {params.r2Filter} \
#             --refpanel {params.refpanel_keyword} \
#             --population all \
#             --build hg19 \
#             --token {params.token} &> {output.confirmation}.log
#         echo "FLAG" > {output.confirmation}
#         """
