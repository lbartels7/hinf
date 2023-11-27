import pandas as pd
# import snakemake

excel_file = snakemake.input[0]

# samples with beta_lactamase status negative and positive
df_negative = pd.read_excel(excel_file, sheet_name='blac_negative')
df_positive = pd.read_excel(excel_file, sheet_name='excluded_blac_positive')


# Keep all samples with MIC values above a defined threshold
MAXIMUM_MIC = 8
df_amp_mic = (df_negative
              .dropna(subset='AMP_MIC')
              .assign(AMP_MIC = lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
            #   .query('AMP_MIC < ' + str(MAXIMUM_MIC))
              .query('AMP_MIC < @MAXIMUM_MIC')
)

IDs_mic = df_amp_mic.SampleID
IDs_mic.to_csv(snakemake.output[0], index=False, header=False)

# Keep all samples with resistance status
df_amp_binary = df_negative.dropna(subset='AMP')
IDs_binary = df_amp_binary.SampleID
IDs_binary.to_csv(snakemake.output[1], index=False, header=False)

