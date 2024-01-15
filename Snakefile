import pandas as pd


configfile: "config/config.yaml"

IDS, = glob_wildcards("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf")
ORIGIN = ['luebeck', 'wurzburg', 'portugal']

# import pandas as pd
# df_negative = (pd.read_excel('metadata_HLR_extern.xlsx', sheet_name='blac_negative')
#                     .assign(AMP_MIC= lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
#                     # .query('AMP_MIC < 200')
# )

# print(df_negative.SampleID.shape)

# container: "docker://continuumio/miniconda3"

def get_all_blac_negative_ids():
    df_negative = (pd.read_excel("vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx", sheet_name='blac_negative')
                        .dropna(subset='AMP_MIC')
                        .assign(AMP_MIC = lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
                        .query('AMP_MIC < 8')
                        .query('SampleID != "HLR-103"')
    )
    return (df_negative['FullID']).tolist()

# print(get_all_blac_negative_ids())

rule all:
    input:  "results/samples_AMP_nonNAN.csv",
            "results/stats.tsv",
            # "results/mash.tsv",
            # "results/resistance_bin.pheno",
            # "results/resistance_mic.pheno",
            # "results/ampicillin_SNPs.txt",    
            "results/Hinf_norm_genes.vcf.gz",
            "results/Hinf_norm_mic.vcf.gz",
            "results/Hinf_norm_genes_mic.vcf.gz",
            "results/Hinf_norm_genes.vcf.gz.csi",
            "results/feather/Hinf_norm.feather",
            "results/feather/Hinf_norm_mic.feather",
            # expand("results/feather/Hinf_norm_mic_{origin}.feather",origin=ORIGIN),
            "results/feather/Hinf_norm_bin.feather",
            "results/feather/Hinf_norm_genes.feather",
            "results/feather/Hinf_norm_genes_mic.feather",
            "results/feather/Hinf_norm_genes_bin.feather",
            # "results/regression/linear_regression.csv",
            # "results/regression/logistic_regression.csv",
            # "results/zarrs/Hinf_norm_mic_linreg_results.zarr",
            # "results/zarrs/Hinf_norm_bin_logreg_results.zarr",
            "results/linreg_unscaled.csv",
            "results/linreg_logscaled.csv",
            "results/linreg_rankscaled.csv",
            "results/logreg.csv",
            # "results/sigsnps_linreg_unscaled.txt",
            # "results/sigsnps_linreg_logscaled.txt",
            # "results/sigsnps_linreg_rankscaled.txt",
            # "results/sigsnps_logreg.txt"


rule extract_AMP_nonNAN:
    input: "vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx"
    output: "results/samples_AMP_MIC_nonNAN.csv",
            "results/samples_AMP_nonNAN.csv"
    params: maximum_mic=config["sample-filtering"]["maximum-mic"]
    script: "scripts/extract_sample_ids.py"



rule calculate_statistics:
    input:  expand("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf", id=IDS)
    output: "results/stats.tsv"
    run:
        import re
        import os
        from contextlib import redirect_stdout
        with open(output[0], 'w') as f_out:
            for filen in input:
                with open(filen, 'r') as f_in:
                    lines = f_in.readlines()
                    datalines = [line for line in lines if re.search('^(?!#).*$',line)]
                    matched_lines_01 = [line for line in lines if re.search('0/1',line)]
                    matched_lines_12 = [line for line in lines if re.search('1/2',line)]
                    # print(os.path.basename(f_in.name), len(datalines), len(matched_lines_01), len(matched_lines_12), sep='\t', file=f_out)
                    f_out.write('\t'.join((os.path.basename(f_in.name), str(len(datalines)), str(len(matched_lines_01)), str(len(matched_lines_12)))) + '\n')


rule vcf_to_vcfgz:
    input: "vcf_Haemophilus/{id}.gatk.vcf"
    output: temp("results/vcf/{id}.gatk.vcf.gz"),
            temp("results/vcf/{id}.gatk.vcf.gz.tbi")
    shell: "bgzip -c {input} > {output}; tabix -p vcf {output[0]}"

rule create_consensus:
    input: vcf="results/vcf/{id}.gatk.vcf.gz",
            index="results/vcf/{id}.gatk.vcf.gz.tbi",
            fasta="vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15.fasta"
    output: "results/pyseer/consensus/{id}.fa"
    shell: "bcftools consensus -H R -f {input.fasta} {input.vcf} > {output}"

rule create_sketch:
    input: expand("results/pyseer/consensus/{id}.fa", id=get_all_blac_negative_ids())
    output: "results/pyseer/samples.msh"
    # conda: "envs/pyseer.yaml"
    shell: "mash sketch -s 1000 -o {output} {input}"

rule create_distance_matrix:
    input: "results/pyseer/samples.msh"
    output: "results/pyseer/mash.tsv"
    conda: "envs/pyseer.yaml"
    shell: "mash dist {input} {input} | square_mash > {output}"

rule rename_samples_in_distance_matrix:
    input: "results/pyseer/mash.tsv"
    output: "results/pyseer/mash_normed.tsv"
    conda: "envs/pandas.yaml"
    script: "scripts/rename_samples_in_distance_matrix.py"


rule create_phenotypefile_pyseer:
    input: "vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx"
    output: "results/pyseer/resistance_bin.pheno",
            "results/pyseer/resistance_mic.pheno",
            "results/pyseer/resistance_mic_log2.pheno"
    run:
        import pandas as pd
        import numpy as np
        ll = ['HLR-453', 'HLR-503', 'HLR-513']
        df_negative = pd.read_excel(input[0], sheet_name='blac_negative')
        (df_negative
            .dropna(subset=['AMP'])
            .rename(columns={'SampleID':'samples'})
            .replace(['S','R'],[0,1])
            .query('samples != "HLR-103"')
            .query('samples not in @ll')
            [['samples', 'AMP']]
            .to_csv(output[0], sep='\t', index=False)
        )

        (df_negative
            .dropna(subset=['AMP_MIC'])
            .rename(columns={'SampleID':'samples'})
            .replace('>8','8')
            .astype({'AMP_MIC': 'float'})
            .query('AMP_MIC <= 8')
            .query('samples != "HLR-103"')
            .query('samples not in @ll')
            [['samples', 'AMP_MIC']]
            .to_csv(output[1], sep='\t', index=False)
        )

        (df_negative
            .dropna(subset=['AMP_MIC'])
            .rename(columns={'SampleID':'samples'})
            .replace('>8','8')
            .astype({'AMP_MIC': 'float'})
            .query('AMP_MIC <= 8')
            .query('samples != "HLR-103"')
            .query('samples not in @ll')
            .assign(AMP_MIC= lambda x: np.log2(x['AMP_MIC']))
            [['samples', 'AMP_MIC']]
            .to_csv(output[2], sep='\t', index=False)
        )

rule pyseer_create_covariate_files:
    input: "vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx"
    output: "results/pyseer/covariates_bin.txt",
            "results/pyseer/covariates_mic.txt",
    run:
        import pandas as pd
        import numpy as np

        df_negative = pd.read_excel(input[0], sheet_name='blac_negative')
        ll = ['HLR-453', 'HLR-503', 'HLR-513']

        (df_negative
            .dropna(subset=['AMP'])
            .rename(columns={'SampleID':'samples'})
            .query('samples not in @ll')
            [['samples', 'origin']]
            .to_csv(output[0], sep='\t', index=False)
        )

        (df_negative
            .dropna(subset=['AMP_MIC'])
            .rename(columns={'SampleID':'samples'})
            .replace('>8','8')
            .astype({'AMP_MIC': 'float'})
            .query('AMP_MIC <= 8')
            .query('samples != "HLR-103"')
            .query('samples not in @ll')
            [['samples', 'origin']]
            .to_csv(output[1], sep='\t', index=False)
        )



# pyseer only uses samples that are present in the phenotypefile, it can filter according to maf, and
# 1/1 0/1 1/0 is interpreted as 1 everthing else as 0
rule pyseer_mic:
    input:  vcf="results/Hinf_norm_temp.vcf",
            pheno="results/pyseer/resistance_mic.pheno",
            dist="results/pyseer/mash_normed.tsv"
    output: "results/pyseer/ampicillin_SNPs.txt"
    conda: "envs/pyseer.yaml"
    shell: "pyseer --phenotypes {input.pheno} --vcf {input.vcf} --distances {input.dist} --max-dimensions 7 --min-af 0.05 --max-af 0.95 > {output}"

# rule create_pyseer_plot:
#     input: "results/ampicillin_SNPs.txt"
#     output: "results/ampicillin_snps.plot"
#     shell: "cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' ampicillin_SNPs.txt | cut -d "_" -f 3) <(sed '1d' ampicillin_SNPs.txt | cut -f 4) | awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) | tr ' ' '\t' > ampicillin_snps.plot"

rule pyseer_samples_for_kinship:
    input: mic="results/pyseer/resistance_mic_log2.pheno",
            bin="results/pyseer/resistance_bin.pheno"
    output: mic="results/pyseer/samples_mic.txt",
            bin="results/pyseer/samples_bin.txt"
    shell: 
            """
            cat {input.mic} | cut -f 1 | tail -n +2 > {output.mic}
            cat {input.bin} | cut -f 1 | tail -n +2 > {output.bin}    
            """


rule pyseer_kinship_matrix:
    input: vcf="results/Hinf_norm_temp.vcf",
            samples_mic= "results/pyseer/samples_mic.txt",
            samples_bin= "results/pyseer/samples_bin.txt"
    output: kinship_mic="results/pyseer/genotype_kinship_mic.tsv",
            kinship_bin="results/pyseer/genotype_kinship_bin.tsv"
    conda: "envs/pyseer.yaml"
    shell: 
            """
            similarity_pyseer --vcf {input.vcf} {input.samples_mic} > {output.kinship_mic}
            similarity_pyseer --vcf {input.vcf} {input.samples_bin} > {output.kinship_bin}
            """


rule pyseer_lmm_mic_log2:
    input: pheno="results/pyseer/resistance_mic_log2.pheno",
            vcf="results/Hinf_norm_temp.vcf",
            kinship="results/pyseer/genotype_kinship_mic.tsv",
            covariates="results/pyseer/covariates_mic.txt"
    output:
            "results/pyseer/ampicillin_snps_mic_log2.txt"
    conda: "envs/pyseer.yaml"
    shell: "pyseer --lmm --phenotypes {input.pheno} --vcf {input.vcf} --similarity {input.kinship} --covariates {input.covariates} --use-covariates 2 --min-af 0.03775 --max-af 0.96225 > {output}"  # based on sample size 267 

rule pyseer_lmm_mic_bin:
    input: pheno="results/pyseer/resistance_bin.pheno",
            vcf="results/Hinf_norm_temp.vcf",
            kinship="results/pyseer/genotype_kinship_bin.tsv",
            covariates="results/pyseer/covariates_bin.txt"
    output:
            "results/pyseer/ampicillin_snps_bin.txt"
    conda: "envs/pyseer.yaml"
    shell: "pyseer --lmm --phenotypes {input.pheno} --vcf {input.vcf} --similarity {input.kinship} --covariates {input.covariates} --use-covariates 2 --min-af 0.03106 --max-af 0.96894 > {output}" # based on sample size 326


rule pyseer_annotate_results_with_geneinformation:
    input: gwas="results/pyseer/ampicillin_snps_mic_log2.txt",
            genes="vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/pyseer/ampicillin_snps_mic_log2_genes.txt"
    conda: "envs/regression.yaml"
    script: "scripts/annotate_with_genes.R"

rule pyseer_annotate_results_with_geneinformation2:
    input: gwas="results/pyseer/ampicillin_snps_bin.txt",
            genes="vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/pyseer/ampicillin_snps_bin_genes.txt"
    conda: "envs/regression.yaml"
    script: "scripts/annotate_with_genes.R"

rule merge_vcfs:
    input:
            vcf=expand("vcf_Haemophilus/HLR-{id}_151bp.gatk.vcf", id=IDS)
    output:
            "results/Hinf.vcf.gz"
    # conda: "envs/bcftools.yaml"
    shell: "bcftools merge --no-index -0 {input.vcf} -Oz -o {output}; tabix -p vcf {output[0]}"


rule normalize:
    input: "results/Hinf.vcf.gz"
    output: "results/Hinf_norm.vcf.gz"
    # conda: "envs/bcftools.yaml"
    shell: "bcftools norm  -m- --multi-overlaps 0 {input} -Oz -o {output}; tabix -p vcf {output[0]}"


# Only include variants that are located within genes. I uncommented the
# first line in the gff file, because it spans the whole genome
rule filter_intergenic_variants:
    input: vcf="results/Hinf_norm.vcf.gz",
            gff="vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15.gff3"
    output: "results/Hinf_norm_genes.vcf.gz"
    # conda: "envs/bedtools.yaml"
    shell: "bedtools intersect -a {input.vcf} -b {input.gff} -header -u | bgzip > {output}"

rule create_index_for_genes:
    input: "results/Hinf_norm_genes.vcf.gz"
    output: "results/Hinf_norm_genes.vcf.gz.csi"
    # conda: "envs/bcftools.yaml"
    shell: "tabix -C {input}"



# HLR-103 not provided as VCF, so use --force-samples
rule filter_samples_whole_genome:
    input: "results/Hinf_norm.vcf.gz",
            "results/samples_AMP_MIC_nonNAN.csv",
            "results/samples_AMP_nonNAN.csv"
    output: "results/Hinf_norm_mic.vcf.gz",
            "results/Hinf_norm_bin.vcf.gz"
    # conda: "envs/bcftools.yaml"
    shell: "bcftools view -S {input[1]} --force-samples {input[0]} -Oz -o {output[0]}; \
            bcftools view -S {input[2]} --force-samples {input[0]} -Oz -o {output[1]}; \
            tabix -p vcf {output[0]}; tabix -p vcf {output[1]}"

# HLR-103 not provided as VCF, so use --force-samples
rule filter_samples_genes:
    input: "results/Hinf_norm_genes.vcf.gz",
            "results/samples_AMP_MIC_nonNAN.csv",
            "results/samples_AMP_nonNAN.csv"
    output: "results/Hinf_norm_genes_mic.vcf.gz",
            "results/Hinf_norm_genes_bin.vcf.gz"
    shell: "bcftools view -S {input[1]} --force-samples {input[0]} -Oz -o {output[0]} \
            && bcftools view -S {input[2]} --force-samples {input[0]} -Oz -o {output[1]}; \
            tabix -p vcf {output[0]}; tabix -p vcf {output[1]}"

rule create_zarrs:
    input: "results/{name}.vcf.gz"
    output: directory("results/zarrs/{name}.zarr")
    run: 
        import sgkit
        from sgkit.io.vcf import vcf_to_zarr
        
        vcf_to_zarr(input[0], output[0])


rule filter_heterozygous_calls_and_map:
    input: "results/zarrs/{name}.zarr",
            "vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx"
    output: directory("results/zarrs/{name}_temp.zarr"),
        #     directory("results/zarrs/{name}_luebeck_temp.zarr"),
        #     directory("results/zarrs/{name}_wurzburg_temp.zarr"),
        #     directory("results/zarrs/{name}_portugal_temp.zarr"),
    params: minimum_allele_count=config["variant-filtering"]["minimum-allele-count"]
    conda: "envs/sgkit.yaml"
    script: "scripts/filter_vcf.py"

rule output_for_r:
    input: "results/zarrs/{name}_temp.zarr",
            # "results/zarrs/{name}_luebeck_temp.zarr",
            # "results/zarrs/{name}_wurzburg_temp.zarr",
            # "results/zarrs/{name}_portugal_temp.zarr",
            "vcf_Haemophilus/annotation_files/metadata_HLR_extern.xlsx"
    output: "results/feather/{name}.feather",
            # "results/feather/{name}_luebeck.feather",
            # "results/feather/{name}_wurzburg.feather",
            # "results/feather/{name}_portugal.feather"
    script: "scripts/convert_to_feather.py"

rule zarr_to_vcf:
    input: "results/zarrs/Hinf_norm_temp.zarr"
    output: "results/Hinf_norm_temp.vcf" # set temp
    run:
        import sgkit as sg
        from sgkit.io.vcf import zarr_to_vcf

        zarr_to_vcf(input[0], output[0])


rule linear_regression:
    input: "results/feather/Hinf_norm_mic.feather",
            "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/regression/linear_regression_unscaled.csv"
    conda: "envs/regression.yaml"
    script: "scripts/linear_regression.R"

rule linear_regression_logscaled:
    input: "results/feather/Hinf_norm_mic.feather",
            "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/regression/linear_regression_logscaled.csv"
    conda: "envs/regression.yaml"
    script: "scripts/linear_regression_logscaled.R"


rule linear_regression_rankscaled:
    input: "results/feather/Hinf_norm_mic.feather",
            "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/regression/linear_regression_rankscaled.csv"
    conda: "envs/regression.yaml"
    script: "scripts/linear_regression_rankscaled.R"

rule logistic_regression:
    input: "results/feather/Hinf_norm_bin.feather",
            "vcf_Haemophilus/annotation_files/Hinf_Rd-KW20v3_DSM11121_2023-06-15_genes_adjst.txt"
    output: "results/regression/logistic_regression.csv"
    conda: "envs/regression.yaml"
    script: "scripts/logistic_regression.R"

# Additionally save the original file as feather
rule create_position_to_mutation_mapping:
    # input: "vcf_Haemophilus/annotation_files/Hflu_association_data_cf4_cr4_fr75_ph8_l0_x0_271_combined_amended_u95_phylo.csv"
    input: "vcf_Haemophilus/annotation_files/HLR_final_logreg_cf4_cr4_fr75_ph8_l0_x1_322_combined_amended.csv"
    output: "results/mapping.csv",
            "results/HLR_final_logreg_cf4_cr4_fr75_ph8_l0_x1_322_combined_amended.feather"
            # "results/Hflu_association_data_cf4_cr4_fr75_ph8_l0_x0_271_combined_amended_u95_phylo.feather"
    conda: "envs/pandas.yaml"
    script: "scripts/create_position_to_mutation_mapping.py" 


rule create_final_linear_regression_result_object:
    input:  mapping="results/mapping.csv",
            gwas_results="results/regression/linear_regression{scaling}.csv",
            input_dataset="results/zarrs/Hinf_norm_mic_temp.zarr/"
    output: directory("results/zarrs/Hinf_norm_mic_linreg{scaling}_results.zarr")
    conda: "envs/sgkit.yaml"
    script: "scripts/annotate_with_mutations.py" 


rule create_final_logistic_regression_result_object:
    input:  mapping="results/mapping.csv",
            gwas_results="results/regression/logistic_regression.csv",
            input_dataset="results/zarrs/Hinf_norm_bin_temp.zarr/"
    output: directory("results/zarrs/Hinf_norm_bin_logreg_results.zarr")
    conda: "envs/sgkit.yaml"
    script: "scripts/annotate_with_mutations.py" 


rule create_results_table_linreg:
    input: "results/zarrs/Hinf_norm_mic_linreg{scaling}_results.zarr"
    output: "results/linreg{scaling}.csv"
    conda: "envs/sgkit.yaml"
    script: "scripts/create_results_table_linreg.py"

rule create_results_table_logreg:
    input: "results/zarrs/Hinf_norm_bin_logreg_results.zarr"
    output: "results/logreg.csv"
    conda: "envs/sgkit.yaml"
    script: "scripts/create_results_table_logreg.py"

rule create_genotype_boxplot:
    input: "results/feather/Hinf_norm_mic.feather"
    output: "results/genotype_boxplots/boxplot_genotype_{variant_id}.pdf"
    log:
        # optional path to the processed notebook
        notebook="logs/notebooks/processed_notebook_{variant_id}.ipynb"
    notebook: "notebooks/genotype_boxplot.r.ipynb"

# https://www.biostars.org/p/250212/
# rule dip_to_hap:
#     input: "results/Hinf_norm_temp.vcf"
#     output: "results/Hinf_norm_hap.vcf"
#     run: 
#         dip_to_hap = {'0/0':'0', '1/1':'1'}

#         with open(input[0],"r") as file, open(output[0],"w") as s_file:
#                 for line in file:
#                     if line[0] == '#':
#                         s_file.write(line)
#                     else:
#                         rowstr = line.split('\t')[9:]
#                         for i in range(len(rowstr)):
#                             rowstr[i] = dip_to_hap[rowstr[i][:3]] + rowstr[i][3:]
#                         row = '\t'.join(str(line).split('\t')[:9]) + '\t' + '\t'.join(rowstr)
#                         s_file.write(row)


rule extract_most_significant_snps_linreg:
    input: "results/linreg{scaling}.csv"
    output: "results/sigsnps_linreg{scaling}.txt"
    shell: "cat {input[0]} | grep  nonsyn | cut -f 1,2,4,10,11,13  | grep CDS | cut -d _ -f 3 | head -n 53 | cut -f 4,5,9 | sort | uniq  | cut -f 1,2 > {output[0]}"





rule extract_most_significant_snps_logreg:
    input: "results/logreg.csv"
    output: "results/sigsnps_logreg.txt"
    shell: "cat {input[0]} | grep  nonsyn | cut -f 1,2,4,9,10,11 | grep CDS | cut -d _ -f 3 | head -n 53 | cut -f 4,5,9 | sort | uniq  | cut -f 1,2 > {output[0]}"