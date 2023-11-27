import pandas as pd

input = snakemake.input[0]
output = snakemake.output[0]

mash = pd.read_csv(input, sep='\t', index_col=0)
mash.index = mash.index.str.split('_').str[0]
mash.columns = mash.columns.str.split('_').str[0]
mash.to_csv(output, sep='\t')