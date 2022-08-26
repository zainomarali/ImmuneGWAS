import pandas as pd

from Variant import Variant
from resources.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output
from resources.eqtlgen import eqtlgen_cis_LDblock_query_formatted_output
from resources.tokyo import tokyo_eqtl_LDblock_query_formatted_output
from helpers.ldlink import ldtrait

"""
Main script for generating a summary report for a given variant.
"""


def generate_summary_report(variant_object: Variant):
    eqtl_cat_df = eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output(variant_object)
    eqtlgen_cis_df = eqtlgen_cis_LDblock_query_formatted_output(variant_object)
    tokyo_df = tokyo_eqtl_LDblock_query_formatted_output(variant_object)
    ldtrait_df = ldtrait(variant_object.get_rsid())

    all_eqtl_df = pd.concat([eqtl_cat_df, eqtlgen_cis_df, tokyo_df], axis=0)

    return all_eqtl_df, ldtrait_df


if __name__ == '__main__':

    # var = Variant("rs149143617", 1, 777870, "C", "G")
    var = Variant.from_rsid('rs9272363')
    var.get_LDblock().to_csv('/home/antton/Desktop/00-LDblock.csv')
    eqtl_df, pheno_df = generate_summary_report(var)
    eqtl_df.to_csv("/home/antton/Desktop/01-eqtl_df.csv", index=False)
    pheno_df.to_csv("/home/antton/Desktop/02-pheno_df.csv", index=False)
