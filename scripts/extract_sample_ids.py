import pandas as pd



metadata_file = snakemake.input[0] # type: ignore
maximum_mic = snakemake.params["maximum_mic"] # type: ignore

output_csv_mic = snakemake.output[0] # type: ignore
output_csv_amp = snakemake.output[1] # type: ignore

# samples with beta_lactamase status negative
df_negative = pd.read_csv(metadata_file, sep='\t')

# Keep all samples with MIC values above a defined threshold
df_amp_mic = (df_negative
              .dropna(subset='AMP_MIC')
              .assign(AMP_MIC = lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
              .query('AMP_MIC <= @maximum_mic')
            #   .query("Institute != 'Lisbon'")
)

IDs_mic = df_amp_mic.ID
IDs_mic.to_csv(output_csv_mic, index=False, header=False) # type: ignore

# Keep all samples with resistance status
df_amp_binary = (df_negative
                 .dropna(subset='AMP')
                #  .query("Institute != 'Lisbon'")
)
IDs_binary = df_amp_binary.ID
IDs_binary.to_csv(output_csv_amp, index=False, header=False) # type: ignore

