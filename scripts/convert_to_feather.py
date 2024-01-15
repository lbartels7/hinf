import pandas as pd
import sgkit as sg

# genotype_data = snakemake.input[0]
# phenotype_data = snakemake.input[4]
# output_file = snakemake.output[0]


for i in range(1):

    genotype_data = snakemake.input[i] # type: ignore
    phenotype_data = snakemake.input[1] # type: ignore
    output_file = snakemake.output[i] # type: ignore
    #  Load genotype data as xarray, set samples as index for the join
    ds_genotype = sg.load_dataset(genotype_data)
    # ds_genotype = ds_genotype.set_index({"samples": "sample_id", 'variants': 'variant_id'})
    ds_genotype = ds_genotype.set_index({'variants': 'variant_id', 'samples': 'sample_id'})


    df_negative = (pd.read_excel(phenotype_data, sheet_name='blac_negative', index_col="SampleID")
                   .assign(AMP_MIC = lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
    )
    # df_positive = pd.read_excel(phenotype_data, sheet_name='excluded_blac_positive', index_col="SampleID")



    # Reduce diploid calls to haploid, all calls should be pseudo homozygous 
    ds_genotype['genotype'] = ds_genotype.call_genotype.max(dim='ploidy')
    df_genotype = ds_genotype.genotype.transpose().to_pandas()

    # join on samples/SampleID (implicit)
    df_final = df_genotype.join(df_negative, how='left')

    # Write df to feather format, first take care of feathers requirements:
    # no index, str columnnames etc.
    df_final.reset_index(inplace=True)
    df_final.columns = df_final.columns.astype('str')
    (
        df_final.astype({'origin': str, 'AMP': str, 'AMP_MIC': str, 'serotype': str, 'beta_lactamase': str})
                .to_feather(output_file)
    )
    # df_final.to_feather(output_file)