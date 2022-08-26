from pandas import ExcelWriter
import pandas as pd

from Variant import Variant
from resources.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multitype
from resources.eqtlgen import eqtlgen_cis_LDblock_query
from resources.tokyo import tokyo_eqtl_LDblock_query
from helpers.ldlink import ldtrait

"""
Main script for generating a summary report for a given variant.
"""


def generate_full_excel_file(variant, output_file):
    with ExcelWriter(f'{output_file}.xlsx') as writer:
        eqtl_cat_df = eqtl_catalogue_LDblock_query_type_restricted_multitype(variant, ['ge', 'microarray'])
        eqtl_cat_df.to_excel(writer, sheet_name='eqtl_cat', index=False)

        eqtlgen_cis_df = eqtlgen_cis_LDblock_query(variant)
        eqtlgen_cis_df.to_excel(writer, sheet_name='eqtlgen_cis', index=False)

        tokyo_df = tokyo_eqtl_LDblock_query(variant)
        tokyo_df.to_excel(writer, sheet_name='tokyo', index=False)

        ldtrait_df = ldtrait(variant.get_rsid())
        ldtrait_df.to_excel(writer, sheet_name='ldtrait', index=False)


if __name__ == '__main__':

    # var = Variant("rs149143617", 1, 777870, "C", "G")
    var = Variant.from_rsid('rs9272363')
    var.get_LDblock().to_csv('/home/antton/Desktop/00-LDblock.csv')
    generate_full_excel_file(var, '/home/antton/Desktop/01-FullReport')
