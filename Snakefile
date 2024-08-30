import pandas as pd


configfile: "config/config.yaml"


(IDS,) = glob_wildcards("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf")
ORIGIN = ["luebeck", "wurzburg", "portugal"]
outdir = config["output-dir"]

IDS_gwas = pd.read_csv('vcf_Haemophilus/annotation_files/HLR-metadata-298.tsv', sep='\t')['libID'].to_list()
(A,B,C,) = glob_wildcards('vcf_Haemophilus/{libid}_{machine}_{bla}_151bp.gatk.vcf')
IDS_gwas = ['_'.join(i) for i in zip(A,B,C) if i[0] in IDS_gwas]

IDS_all = pd.read_csv('vcf_Haemophilus/annotation_files/HLR-metadata-322.tsv', sep='\t')['libID'].to_list()
(A,B,C,) = glob_wildcards('vcf_Haemophilus/{libid}_{machine}_{bla}_151bp.gatk.vcf')
IDS_all = ['_'.join(i) for i in zip(A,B,C) if i[0] in IDS_all]



rule all:
    input:
        outdir + "/samples_AMP_nonNAN.csv",
        outdir + "/stats.tsv",
        outdir + "/Hinf_norm_genes.vcf.gz",
        outdir + "/Hinf_norm_mic.vcf.gz",
        outdir + "/Hinf_norm_genes_mic.vcf.gz",
        outdir + "/Hinf_norm_genes.vcf.gz.csi",
        outdir + "/feather/Hinf_norm.feather",
        outdir + "/feather/Hinf_norm_mic.feather",
        outdir + "/feather/Hinf_norm_bin.feather",
        outdir + "/feather/Hinf_norm_genes.feather",
        outdir + "/feather/Hinf_norm_genes_mic.feather",
        outdir + "/feather/Hinf_norm_genes_bin.feather",
        outdir + "/linreg_logscaled.csv",
        outdir + "/logreg.csv",

rule extract_AMP_nonNAN:
    input:
        "vcf_Haemophilus/annotation_files/HLR-metadata-298.tsv",
    output:
        outdir + "/samples_AMP_MIC_nonNAN.csv",
        outdir + "/samples_AMP_nonNAN.csv",
    params:
        maximum_mic=config["sample-filtering"]["maximum-mic"],
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/extract_sample_ids.py"


rule calculate_statistics:
    input:
        expand("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf", id=IDS),
    output:
        outdir + "/stats.tsv",
    run:
        import re
        import os
        from contextlib import redirect_stdout

        with open(output[0], "w") as f_out:
            for filen in input:
                with open(filen, "r") as f_in:
                    lines = f_in.readlines()
                    datalines = [line for line in lines if re.search("^(?!#).*$", line)]
                    matched_lines_01 = [
                        line for line in lines if re.search("0/1", line)
                    ]
                    matched_lines_12 = [
                        line for line in lines if re.search("1/2", line)
                    ]
                    # print(os.path.basename(f_in.name), len(datalines), len(matched_lines_01), len(matched_lines_12), sep='\t', file=f_out)
                    f_out.write(
                        "\t".join(
                                (
                                    os.path.basename(f_in.name),
                                    str(len(datalines)),
                                    str(len(matched_lines_01)),
                                    str(len(matched_lines_12)),
                                )
                            )
                        + "\n"
                    )

use rule calculate_statistics as calculate_statistics_gwas_cohort with:
    input:
        expand("vcf_Haemophilus/{id}_151bp.gatk.vcf", id=IDS_gwas)
    output:
        outdir + "/stats_gwas.tsv"

rule vcf_to_vcfgz:
    input:
        "vcf_Haemophilus/{id}.gatk.vcf",
    output:
        temp(outdir + "/vcf/{id}.gatk.vcf.gz"),
        temp(outdir + "/vcf/{id}.gatk.vcf.gz.tbi"),
    conda: "envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}; tabix -p vcf {output[0]}"


rule merge_vcfs:
    input:
        vcf=expand("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf", id=IDS),
    output:
        outdir + "/Hinf.vcf.gz",
    conda: "envs/bcftools.yaml"
    shell:
        "bcftools merge --no-index -0 {input.vcf} -Oz -o {output}; tabix -p vcf {output[0]}"


rule normalize:
    input:
        outdir + "/Hinf.vcf.gz",
    output:
        outdir + "/Hinf_norm.vcf.gz",
    conda: "envs/bcftools.yaml"
    shell:
        "bcftools norm  -m- --multi-overlaps 0 {input} -Oz -o {output}; tabix -p vcf {output[0]}"


# Only include variants that are located within genes. I uncommented the
# first line in the gff file, because it spans the whole genome
rule filter_intergenic_variants:
    input:
        vcf=outdir + "/Hinf_norm.vcf.gz",
        gff="vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15.gff3",
    output:
        outdir + "/Hinf_norm_genes.vcf.gz",
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.vcf} -b {input.gff} -header -u | bgzip > {output}"


rule create_index_for_genes:
    input:
        outdir + "/Hinf_norm_genes.vcf.gz",
    output:
        outdir + "/Hinf_norm_genes.vcf.gz.csi",
    conda: "envs/bcftools.yaml"
    shell:
        "tabix -C {input}"



rule filter_samples_whole_genome:
    input:
        outdir + "/Hinf_norm.vcf.gz",
        outdir + "/samples_AMP_MIC_nonNAN.csv",
        outdir + "/samples_AMP_nonNAN.csv",
    output:
        outdir + "/Hinf_norm_mic.vcf.gz",
        outdir + "/Hinf_norm_bin.vcf.gz",
    conda: "envs/bcftools.yaml"
    shell:
        "bcftools view -S {input[1]} --force-samples {input[0]} -Oz -o {output[0]}; \
            bcftools view -S {input[2]} --force-samples {input[0]} -Oz -o {output[1]}; \
            tabix -p vcf {output[0]}; tabix -p vcf {output[1]}"


rule filter_samples_genes:
    input:
        outdir + "/Hinf_norm_genes.vcf.gz",
        outdir + "/samples_AMP_MIC_nonNAN.csv",
        outdir + "/samples_AMP_nonNAN.csv",
    output:
        outdir + "/Hinf_norm_genes_mic.vcf.gz",
        outdir + "/Hinf_norm_genes_bin.vcf.gz",
    conda: "envs/bcftools.yaml"
    shell:
        "bcftools view -S {input[1]} --force-samples {input[0]} -Oz -o {output[0]} \
            && bcftools view -S {input[2]} --force-samples {input[0]} -Oz -o {output[1]}; \
            tabix -p vcf {output[0]}; tabix -p vcf {output[1]}"


rule create_zarrs:
    input:
        outdir + "/{name}.vcf.gz",
    output:
        directory(outdir + "/zarrs/{name}.zarr"),
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/create_zarrs.py"


rule filter_heterozygous_calls_and_map:
    input:
        outdir + "/zarrs/{name}.zarr",
        "vcf_Haemophilus/annotation_files/HLR-metadata-298.tsv",
        "vcf_Haemophilus/annotation_files/amended_positions.txt"
    output:
        directory(outdir + "/zarrs/{name}_temp.zarr"),
    params:
        minimum_allele_count=config["variant-filtering"]["minimum-allele-count"],
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/filter_vcf.py"


rule output_for_r:
    input:
        outdir + "/zarrs/{name}_temp.zarr",
        "vcf_Haemophilus/annotation_files/HLR-metadata-298.tsv",
    output:
        outdir + "/feather/{name}.feather",
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/convert_to_feather.py"


rule zarr_to_vcf:
    input:
        outdir + "/zarrs/Hinf_norm_temp.zarr",
    output:
        outdir + "/Hinf_norm_temp.vcf",  # set temp
    run:
        import sgkit as sg
        from sgkit.io.vcf import zarr_to_vcf

        zarr_to_vcf(input[0], output[0])



rule linear_regression_logscaled:
    input:
        outdir + "/feather/Hinf_norm_mic.feather",
        "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt",
    output:
        outdir + "/regression/linear_regression_logscaled.csv",
    conda:
        "envs/regression.yaml"
    script:
        "scripts/linear_regression_logscaled.R"



rule logistic_regression:
    input:
        outdir + "/feather/Hinf_norm_bin.feather",
        "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt",
    output:
        outdir + "/regression/logistic_regression.csv",
    conda:
        "envs/regression.yaml"
    script:
        "scripts/logistic_regression.R"


# Additionally save the original file as feather
rule create_position_to_mutation_mapping:
    input:
        "vcf_Haemophilus/annotation_files/HLR_final_logreg_cf4_cr4_fr75_ph8_l0_x1_322_combined_amended.csv",
    output:
        outdir + "/mapping.csv",
        outdir + "/HLR_final_logreg_cf4_cr4_fr75_ph8_l0_x1_322_combined_amended.feather",
    conda:
        "envs/pandas.yaml"
    script:
        "scripts/create_position_to_mutation_mapping.py"

# Use rule inheritance maybe
rule create_final_linear_regression_result_object:
    input:
        mapping=outdir + "/mapping.csv",
        gwas_results=outdir + "/regression/linear_regression{scaling}.csv",
        input_dataset=outdir + "/zarrs/Hinf_norm_mic_temp.zarr/",
    output:
        directory(outdir + "/zarrs/Hinf_norm_mic_linreg{scaling}_results.zarr"),
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/annotate_with_mutations.py"


rule create_final_logistic_regression_result_object:
    input:
        mapping=outdir + "/mapping.csv",
        gwas_results=outdir + "/regression/logistic_regression.csv",
        input_dataset=outdir + "/zarrs/Hinf_norm_bin_temp.zarr/",
    output:
        directory(outdir + "/zarrs/Hinf_norm_bin_logreg_results.zarr"),
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/annotate_with_mutations.py"

# Rule inheritance also
rule create_results_table_linreg:
    input:
        outdir + "/zarrs/Hinf_norm_mic_linreg{scaling}_results.zarr",
    output:
        outdir + "/linreg{scaling}.csv",
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/create_results_table_linreg.py"


rule create_results_table_logreg:
    input:
        outdir + "/zarrs/Hinf_norm_bin_logreg_results.zarr",
    output:
        outdir + "/logreg.csv",
    conda:
        "envs/sgkit.yaml"
    script:
        "scripts/create_results_table_logreg.py"


rule extract_ftsI_vcf:
    input: outdir + "/zarrs/Hinf_norm_mic_linreg_logscaled_results.zarr/"
    output: outdir + "/vcf/linreg_nonsynonymous_ftsI.vcf"
    conda: "envs/sgkit.yaml"
    script: "scripts/extract_ftsI_vcf.py"

rule vcf2plink:
    input: outdir + "/vcf/linreg_nonsynonymous_ftsI.vcf"
    output: temp(multiext(outdir + "/tmp/linreg_nonsynonymous_ftsI.", "ped", "log", "map"))
    params: prefix=lambda wildcards, output: output[0][:-4]
    conda: "envs/vcftools.yaml"
    shell: "vcftools --vcf {input[0]} --plink --out {params.prefix}"

rule ped2bed:
    input: multiext(outdir + "/tmp/linreg_nonsynonymous_ftsI.", "ped", "log", "map")
    output: multiext(outdir + "/bed/linreg_nonsynonymous_ftsI.", "bed", "bim", "fam")
    params:
        in_prefix=lambda wildcards, input: input[0][:-4],
        out_prefix=lambda wildcards, output: output[0][:-4]
    conda: "envs/plink.yaml"
    shell: "plink --file {params.in_prefix} --make-bed --out {params.out_prefix}"

rule calc_ld_plink:
    input: outdir + "/bed/linreg_nonsynonymous"
    output: outdir + "/ld/ld_results"
    conda: "envs/plink.yaml"
    shell: "plink --bfile {input[0]} --ld-window-kb 10000 --ld-window-r2 0 --r2 --out {output[0]}"
