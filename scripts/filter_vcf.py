import sgkit as sg
from sgkit.io.vcf import vcf_to_zarr

import numpy as np
import pandas as pd
import xarray as xr

input_zarr = snakemake.input[0] # type: ignore
excel_file = snakemake.input[1] # type: ignore
minimum_allele_count = snakemake.params["minimum_allele_count"] # type: ignore

output_zarr = snakemake.output[0] # type: ignore

def remove_heterozgous_variants(ds : xr.Dataset) -> xr.Dataset:
    
    # variables_to_drop = ['sample_n_called', 'sample_call_rate', 'sample_n_het', 'sample_n_hom_ref',
    #                       'sample_n_hom_alt', 'sample_n_non_ref', 'variant_n_called', 'variant_call_rate',
    #                       'variant_n_het', 'variant_n_hom_ref', 'variant_n_hom_alt', 'variant_n_non_ref',
    #                       'variant_allele_count', 'variant_allele_total', 'variant_allele_frequency']

    return(
        ds
        # .pipe(sg.variant_stats)
        # .pipe(sg.sample_stats)
        .pipe(lambda ds: ds.sel(variants=(ds.call_genotype[:,:,0] == ds.call_genotype[:,:,1]).all(dim='samples').compute()))
        # .drop(variables_to_drop)
        )

def filter_minimum_AC(ds: xr.Dataset, minimum_AC: int = 20) -> xr.Dataset:

    variables_to_drop = ['variant_n_called', 'variant_call_rate',
                          'variant_n_het', 'variant_n_hom_ref', 'variant_n_hom_alt', 'variant_n_non_ref',
                          'variant_allele_count', 'variant_allele_total', 'variant_allele_frequency']

    return(
        ds
        .pipe(sg.variant_stats)
        # .pipe(sg.sample_stats)
        .pipe(lambda ds: ds.sel(variants=(ds.variant_allele_count[:,:2].min(dim='alleles') >= minimum_AC).compute()))
        .drop(variables_to_drop)
        )

def filter_maximum_hetero(ds: xr.Dataset, maximum_hetero: int = 10_000) -> xr.Dataset:

    variables_to_drop = ['sample_n_called', 'sample_call_rate', 'sample_n_het', 'sample_n_hom_ref', 'sample_n_hom_alt',
                          'sample_n_non_ref']

    return(
        ds
        .pipe(sg.sample_stats)
        .pipe(lambda ds: ds.sel(samples=(ds.sample_n_het < 10_000).compute()))
        .drop(variables_to_drop)
        )

# MAximum MIC is 256, so 300 means no filtering
# def filter_maximum_MIC(ds: xr.Dataset, maximum_mic: int = 300) -> xr.Dataset:
#     return(
#         ds
#         .pipe(lambda ds: ds.sel(samples=(ds.AMP_MIC < maximum_mic).compute()))
#         )

# Drop all variables with the ploidy dimension + variant_id_mask and reduce call_genotype to haploid
ds = (sg.load_dataset(input_zarr) 
      .drop_vars(['call_genotype_mask', 'call_genotype_phased', 'variant_id_mask'])
      .assign(call_genotype = lambda ds: ds.call_genotype.max(dim='ploidy')))

# Recreate the ploidy dimension, because downstream method calls expect it to be present
# Same for the genotype mask
ds['call_genotype'] = ds.call_genotype.expand_dims({'ploidy':1},2)
ds['call_genotype_mask'] = xr.zeros_like(ds.call_genotype).astype(bool)


# samples with beta_lactamase status negative and positive
df_negative = pd.read_excel(excel_file, sheet_name='blac_negative').assign(AMP_MIC= lambda df: df['AMP_MIC'].replace('>8','8').astype('float'))
# df_positive = pd.read_excel(excel_file, sheet_name='excluded_blac_positive').assign(AMP_MIC= lambda df: df['AMP_MIC'].replace(['>8', '>256'],['8', '256']).astype('float'))

ds_negative = df_negative.to_xarray().rename({"SampleID":"samples"}).swap_dims({'index': 'samples'}).drop(labels='index')
ds = ds.set_index({"samples": "sample_id"})
ds = ds.merge(ds_negative, join="left")

# ds = ds.reset_index("samples").reset_coords(drop=True)
# ds = (ds.pipe(sg.sample_stats)
#       .pipe(sg.variant_stats))

# mask_n_het = (ds.sample_n_het < 10000).values
# ds = ds.sel(samples=mask_n_het)

# variant and sample stats recalculate
# mask_allele_freq = (ds.variant_allele_frequency[:,:2].min(dim='alleles') > 0.1)
# mask_allele_count = (ds.variant_allele_count[:,:2].min(dim='alleles') > 20)
# mask_hom = (ds.call_genotype[:,:,0] == ds.call_genotype[:,:,1]).all(dim='samples')
# ds =  ds.sel(variants=((mask_hom & mask_allele_count).compute()))

ds = ds.reset_index('samples').reset_coords('samples').rename_vars({'samples': 'sample_id'})
# ds_luebeck = ds.sel(samples =(ds.origin == 'Lübeck'))
# ds_wurzburg = ds.sel(samples =(ds.origin == 'Wurzburg'))
# ds_portugal = ds.sel(samples =(ds.origin == 'Portugal'))


ds = (
    ds
    # .pipe(filter_maximum_hetero,10_000)
    # .pipe(remove_heterozgous_variants)
    .pipe(filter_minimum_AC, minimum_allele_count)
    .assign(variant_id=lambda ds: ( 'ID_' + ds.variant_position.to_series().astype(str) + '_' + ds.variant_allele[:,1].to_series()))
)

# ds_luebeck = (
#     ds_luebeck
#     .pipe(filter_maximum_hetero,10_000)
#     .pipe(remove_heterozgous_variants)
#     .pipe(filter_minimum_AC,20)
#     .assign(variant_id=lambda ds: ( ds.variant_position.to_series().astype(str) + '_' + ds.variant_allele[:,1].to_series()))
# )

# ds_wurzburg = (
#     ds_wurzburg
#     .pipe(filter_maximum_hetero,10_000)
#     .pipe(remove_heterozgous_variants)
#     .pipe(filter_minimum_AC,20)
#     .assign(variant_id=lambda ds: ( ds.variant_position.to_series().astype(str) + '_' + ds.variant_allele[:,1].to_series()))
# )

# ds_portugal = (
#     ds_portugal
#     .pipe(filter_maximum_hetero,10_000)
#     .pipe(remove_heterozgous_variants)
#     .pipe(filter_minimum_AC,20)
#     .assign(variant_id=lambda ds: ( ds.variant_position.to_series().astype(str) + '_' + ds.variant_allele[:,1].to_series()))
# )

# ds['variant_id'].values = np.char.add(np.char.add(ds.variant_position.astype(str), '_'), ds.variant_allele[:,1].astype(str))


# Save to zarr
sg.save_dataset(ds, output_zarr, auto_rechunk=True) # type: ignore
# sg.save_dataset(ds_luebeck, snakemake.output[1], auto_rechunk=True)
# sg.save_dataset(ds_wurzburg, snakemake.output[2], auto_rechunk=True)
# sg.save_dataset(ds_portugal, snakemake.output[3], auto_rechunk=True)



