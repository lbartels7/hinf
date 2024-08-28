import sgkit as sg
from sgkit.io.vcf import write_vcf, zarr_to_vcf


input_zarr = snakemake.input[0] # type: ignore
output_vcf = snakemake.output[0] # type: ignore

ds_results = sg.load_dataset(input_zarr)

# Extract just the ftsI positions
# ftsI 1688288 to 1690120, transpeptidase 1688533 to 1689074
ds_ftsI = ds_results.sel(variants=((1688288 < ds_results.variant_position) & (ds_results.variant_position < 1690120)))

# We just want to have the nonsynonymous variants
ref = ds_ftsI['mutation'].to_series().str.split().str[0].str[0]
alt = ds_ftsI['mutation'].to_series().str.split().str[0].str[-1]
mask_mut = ~(ref == alt)

ds_ftsI_nonsym = ds_ftsI.sel(variants=mask_mut.values)

# artificially fix the zarr, so that a vcf can be generated
# not perfect, but should have no consequences on the subsequent steps
ds_ftsI.attrs = {'contig_lengths': [1830701],
                 'contigs': ['Hinf_Rd-KW20v3'],
                 'filters': ['PASS', 'LowQual'],
                 'max_alt_alleles_seen': 1,
                 'source': 'sgkit-0.8.0',
                 'vcf_zarr_version': '0.2'}

write_vcf(ds_ftsI, output_vcf)