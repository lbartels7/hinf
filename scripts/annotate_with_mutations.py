import pandas as pd
import sgkit as sg

df_mapping = pd.read_csv(snakemake.input[0]) # type: ignore
df_gwas = pd.read_csv(snakemake.input[1]) # type: ignore
ds_input = sg.load_dataset(snakemake.input[2]) # genotype information etc. # type: ignore



# Merging the gwas resulst with the mapping table on the id (id is {position}_{alternative allele}), so we have
# the corresponding mutations for each SNP
df_gwas_mutation = (pd.merge(df_gwas, df_mapping, how='left', on='id', validate='one_to_one')
        .drop(columns=['position']) # 'pvalue', 'p_values.adj', 'effect', 'r.squared', 'adj.r.squared', 'mutation'
)
# could be a mistake, maybee I wanted to remove the whole column instead
# df_gwas_mutation = df_gwas_mutation.dropna(subset='allele')





# Merging the regression results with the starting xarray which holds all the
# other variant- and sample-related data, like MIC and variant calls
ds_input = (ds_input.reset_index('variants').reset_coords('variants')
            .set_index({'variants': 'variant_id'})
            .set_coords({'variants': 'variant_id'})
)




ds_gwas = (df_gwas_mutation.to_xarray()
           .reset_index('index').reset_coords('index').rename_dims({'index' : 'variants'})
           .rename_vars({'id': 'variant_id'}).set_coords('variant_id').set_index({'variants' : 'variant_id'})
           .drop_vars('index')
)


ds_final = (ds_gwas.merge(ds_input, join='left')
            .drop(['allele', 'Unnamed: 0_x', 'Unnamed: 0_y', 'pos'])
)

sg.save_dataset(ds_final, snakemake.output[0], auto_rechunk=True) # type: ignore

