import pandas as pd



excel_file = snakemake.input[0] # type: ignore
maximum_mic = snakemake.params["maximum_mic"] # type: ignore

output_csv_mic = snakemake.output[0] # type: ignore
output_csv_amp = snakemake.output[1] # type: ignore

# samples with beta_lactamase status negative and positive
df_negative = pd.read_excel(excel_file, sheet_name='blac_negative')
df_positive = pd.read_excel(excel_file, sheet_name='excluded_blac_positive')

# Isolates that harbor bla genes, identified in hindsight
bla_pos = ['HLR-453', 'HLR-503', 'HLR-513']
# Keep all samples with MIC values above a defined threshold
df_amp_mic = (df_negative
              .dropna(subset='AMP_MIC')
              .assign(AMP_MIC = lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
            #   .query('AMP_MIC < ' + str(MAXIMUM_MIC))
              .query('AMP_MIC <= @maximum_mic')
              .query("origin == 'Portugal'")
              .query('SampleID not in @bla_pos')
)

IDs_mic = df_amp_mic.SampleID
IDs_mic.to_csv(output_csv_mic, index=False, header=False) # type: ignore

# Keep all samples with resistance status
df_amp_binary = (df_negative
                 .dropna(subset='AMP')
                 .query("origin == 'Portugal'")
                 .query('SampleID not in @bla_pos')
)
IDs_binary = df_amp_binary.SampleID
IDs_binary.to_csv(output_csv_amp, index=False, header=False) # type: ignore

