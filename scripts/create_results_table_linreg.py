import sgkit as sg
import pandas as pd


input = snakemake.input[0]
output = snakemake.output[0]
ds_linreg = sg.load_dataset(input)


df_variants = (ds_linreg[['variants', 'variant_position', 'variant_quality', 'pvalue', 'p_values.adj', 'mutation', 'r.squared', 'adj.r.squared', 'effect', 'gene.name', 'gene.product', 'gene.type']]
            .to_pandas()
            .assign(mutation_type= lambda df: ((df['mutation'] == '') | (df['mutation'].str.split().str[0].str[0] == df['mutation'].str.split().str[0].str[-1]))
                    .apply(lambda x: 'synonymous' if x else 'nonsynonymous'))
)
df_genotypes = ds_linreg.set_index({'samples': 'sample_id'}).call_genotype.squeeze().to_pandas()


df_allele_name = ds_linreg.variant_allele.to_pandas().drop([2,3], axis=1).rename(columns={0: 'ref', 1: 'alt'})
df_allele_count = ds_linreg.pipe(sg.variant_stats).variant_allele_count.to_pandas().drop([2,3], axis=1).rename(columns={0: 'ref_count', 1: 'alt_count'})
df_allele = pd.merge(df_allele_name,df_allele_count, on='variants')

bigtable = pd.merge(df_variants, df_allele, on='variants')
bigtable = pd.merge(df_variants, df_genotypes, on='variants').sort_values('pvalue')
bigtable.to_csv(output, sep='\t')