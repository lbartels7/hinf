import pandas as pd

csv = pd.read_csv(snakemake.input[0], sep='\t', skiprows=[0], na_values=['-'], nrows=10) # type: ignore

variant_cols = ['#Pos.', 'Ref', 'Gene', 'GeneName', 'Annotation', 'Main_type', 'Main_allel']
strain_cols = csv.columns.to_series().filter(regex='^(HLR)').tolist()

# csv = pd.read_csv(snakemake.input[0], sep='\t', skiprows=[0], na_values=['-'], usecols=(variant_cols + strain_cols), low_memory=False)




csv = (pd.read_csv(snakemake.input[0], sep='\t', skiprows=[0], na_values=['-'], # type: ignore
                    usecols=(variant_cols + strain_cols), low_memory=False)
            .rename(columns={'#Pos.': 'position'})
            .filter(regex='(position)|^(HLR)', axis=1)
)

csv.to_feather(snakemake.output[1]) # type: ignore

allele_cols = csv.columns.to_series().filter(regex='^(HLR).*(bp)$')
alleles_df = (csv
        .melt(id_vars=['position'], value_vars=allele_cols) 
        .set_index(['position', 'variable'])
)

mutation_cols = csv.columns.to_series().filter(regex='^(HLR).*(1)$')
mutations_df = (csv
         .melt(id_vars=['position'], value_vars=mutation_cols)
         .assign(variable= lambda x: x['variable'].str.slice(stop=-2))
         .set_index(['position', 'variable'])
)


mapping = (pd.merge(alleles_df, mutations_df, left_index=True, right_index=True)
           .reset_index()
           .drop(columns='variable')
           .rename(columns={'value_x': 'allele', 'value_y' : 'mutation'})
           .assign(allele= lambda df: df['allele'].str.upper())
           .drop_duplicates()
           .dropna(subset='mutation')
           .assign(id = lambda df: 'ID_' + df['position'].astype(str) + '_' + df['allele'])
)

mapping.to_csv(snakemake.output[0]) # type: ignore
