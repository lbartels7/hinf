import sgkit as sg
import pandas as pd


input = snakemake.input[0] # type: ignore
output = snakemake.output[0] # type: ignore
ds_logreg = sg.load_dataset(input)


df_variants = (ds_logreg[['variants', 'variant_position', 'variant_quality', 'pvalue', 'p_values.adj', 'mutation', 'odds.ratio', 'effect',  'gene.name', 'gene.product', 'gene.type', 'gene.id']]
            .to_pandas()
            .assign(mutation_type= lambda df: ((df['mutation'] == '') | (df['mutation'].str.split().str[0].str[0] == df['mutation'].str.split().str[0].str[-1]))
                    .apply(lambda x: 'synonymous' if x else 'nonsynonymous'))
)

df_genotypes = ds_logreg.set_index({'samples': 'sample_id'}).call_genotype.squeeze().to_pandas()


df_allele_name = (ds_logreg.variant_allele
                  .to_pandas()
                  .drop([2,3], axis=1)
                  .rename(columns={0: 'ref', 1: 'alt'})
)

df_allele_count = (ds_logreg.pipe(sg.variant_stats).variant_allele_count
                   .to_pandas()
                   .drop([2,3], axis=1)
                   .rename(columns={0: 'ref_count', 1: 'alt_count'})
)

df_allele = pd.merge(df_allele_name,df_allele_count, on='variants')

bigtable = pd.merge(df_variants, df_allele, on='variants')
bigtable = pd.merge(df_variants, df_genotypes, on='variants').sort_values('pvalue')
bigtable.to_csv(output, sep='\t')