#!/usr/bin/env python3


configfile = "config.yaml"


def get_refpanel_info(refpanel):
    if refpanel == 'topmed':
        out = {
            "url": "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2",
            "post": "/jobs/submit/imputationserver",
            "refpanel": "apps@topmed-r2@1.0.0",
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
        "results/imputation-topmed-r2filter_0.3/mis.submitted",
        "results/imputation-1000g-r2filter_0/mis.submitted",
        "results/imputation_hla-hla_4d-r2filter_0/mis.submitted",
        "results/imputation_hla-hla_g-r2filter_0/mis.submitted",



rule imputation_prep:
    input:
        bim = "results/plink/formatted.bim",
        frq = "results/plink/formatted.frq",
        reference_1000g = "resources/1000GP_Phase3_combined.legend"
    output:
        "results/imputation_prep/Run-plink.sh",
        expand("results/imputation_prep/formatted-updated-chr{chr}.vcf", chr=range(1,23))
    params:
        outdir = "results/imputation_prep",
    shell:
        """
        scripts/raynor_imputation_prep/HRC-1000G-check-bim-NoReadKey.pl --bim {input.bim} --frequency {input.frq} --ref {input.reference_1000g} --1000g --output {params.outdir} --verbose --noexclude
        grep -v '^rm' {params.outdir}/Run-plink.sh | sed 's/--recode vcf/--recode vcf-iid/g'> {params.outdir}/Run-plink_fix.sh
        bash {params.outdir}/Run-plink_fix.sh
        """


rule zip_and_index:
    input:
        vcfs =expand("results/imputation_prep/formatted-updated-chr{chr}.vcf", chr=range(1,23)),
    output:
        gz = expand("results/imputation_prep/formatted-updated-chr{chr}.vcf.gz", chr=range(1,23)),
        tbi = expand("results/imputation_prep/formatted-updated-chr{chr}.vcf.gz.tbi", chr=range(1,23)),
    conda:
        "general"
    shell:
        """
        for i in {input.vcfs}; do bgzip $i; tabix ${{i}}.gz ; done;
        """

rule filter_1000g:
    input:
        vcf_array = "results/imputation_prep/formatted-updated-chr{chr}.vcf.gz",
        vcf_1000g = lambda wildcards: config["KGref_hg19"][wildcards.chr],
    output:
        vcf = "results/imputation_prep/1000g_bctfools_isec_chr{chr}/0001.vcf.gz",
        tbi = "results/imputation_prep/1000g_bctfools_isec_chr{chr}/0001.vcf.gz.tbi",
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
        vcf_array = "results/imputation_prep/formatted-updated-chr{chr}.vcf.gz",
        vcf_1000g = "results/imputation_prep/1000g_bctfools_isec_chr{chr}/0001.vcf.gz",
    output:
        vcf = "results/imputation_prep/formatted-updated-add1000g-chr{chr}.vcf.gz",
        tbi = "results/imputation_prep/formatted-updated-add1000g-chr{chr}.vcf.gz.tbi",
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
        vcfs = expand("results/imputation_prep/formatted-updated-add1000g-chr{chr}.vcf.gz", chr=range(1,23))
    output:
        confirmation = "results/imputation-{refpanel}-r2filter_{r2filter}/mis.submitted",
    params:
        url = lambda wildcards: get_refpanel_info(wildcards.refpanel)['url'],
        post = lambda wildcards: get_refpanel_info(wildcards.refpanel)['post'],
        refpanel_keyword = lambda wildcards: get_refpanel_info(wildcards.refpanel)['refpanel'],
        token = lambda wildcards: get_refpanel_info(wildcards.refpanel)['token'],
        vcfstring = ' '.join(expand("results/imputation_prep/formatted-updated-add1000g-chr{chr}.vcf.gz", chr=range(1,23))),
        r2filter = "{r2filter}",
    conda:
        "general"
    shell:
        """
        python scripts/submit-mis-CCR.py \
            --vcf {params.vcfstring} \
            --url {params.url} \
            --post {params.post} \
            --r2filter {params.r2filter} \
            --refpanel {params.refpanel_keyword} \
            --population mixed \
            --build hg19 \
            --mode imputation \
            --token {params.token} &> {output.confirmation}.log
        echo "FLAG" > {output.confirmation}
        """

rule impute_hla:
    input:
        vcf = "results/imputation_prep/formatted-updated-add1000g-chr6.vcf.gz",
    output:
        confirmation = "results/imputation_hla-{refpanel}-r2filter_{r2filter}/mis.submitted",
    params:
        url = lambda wildcards: get_refpanel_info(wildcards.refpanel)['url'],
        post = lambda wildcards: get_refpanel_info(wildcards.refpanel)['post'],
        refpanel_keyword = lambda wildcards: get_refpanel_info(wildcards.refpanel)['refpanel'],
        token = lambda wildcards: get_refpanel_info(wildcards.refpanel)['token'],
        r2filter = "{r2filter}",
    conda:
        "general"
    shell:
        """
        python scripts/submit-mis-CCR.py \
            --vcf {input.vcf} \
            --url {params.url} \
            --post {params.post} \
            --r2filter {params.r2filter} \
            --refpanel {params.refpanel_keyword} \
            --population mixed \
            --build hg19 \
            --mode imputation \
            --token {params.token} &> {output.confirmation}.log
        echo "FLAG" > {output.confirmation}
        """
