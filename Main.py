import os
import pandas as pd

from Variant import Variant
from helpers.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output
from helpers.eqtlgen import eqtlgen_cis_LDblock_query_formatted_output
from helpers.tokyo import tokyo_eqtl_LDblock_query_formatted_output, tokyo_eqtl_LDblock_query
from helpers.ldlink import ldtrait

"""
Main script for generating a summary report for a given variant.
"""

#var = Variant("rs149143617", 1, 777870, "C", "G")
var = var = Variant.from_rsid('rs9272363')


def generate_summary_report(variant_object: Variant):
    eqtl_cat_df = eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output(variant_object)
    eqtlgen_cis_df = eqtlgen_cis_LDblock_query_formatted_output(variant_object)
    tokyo_df = tokyo_eqtl_LDblock_query_formatted_output(variant_object)
    ldtrait_df = ldtrait(variant_object.get_rsid())

    all_eqtl_df = pd.concat([eqtl_cat_df, eqtlgen_cis_df, tokyo_df], axis=0)

    return all_eqtl_df, ldtrait_df


eqtl_df, pheno_df = generate_summary_report(var)
eqtl_df.to_csv("/home/antton/Desktop/eqtl_df.csv", index=False)
pheno_df.to_csv("/home/antton/Desktop/pheno_df.csv", index=False)
