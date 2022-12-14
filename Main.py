from pandas import ExcelWriter
import pandas as pd
import logging

from ImmuneGWAS.variant import Variant
from ImmuneGWAS.resources.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multitype
from ImmuneGWAS.resources.eqtlgen import eqtlgen_cis_LDblock_query
from ImmuneGWAS.resources.tokyo_eqtl import tokyo_eqtl_LDblock_query
from ImmuneGWAS.resources.tokyo_ge import tokyo_ge_query
from ImmuneGWAS.helpers.ldlink import ldtrait
from ImmuneGWAS.config import output_folder
from ImmuneGWAS.plotting import plot_tokyo_ge

"""
Main script for generating a summary report for a given variant.
"""


def generate_full_excel_file(variant, output_file=output_folder+'full_report.xlsx'):
    """
    Generate an Excel file with all the information for a given variant.

    :param variant: Variant object to generate the report for.
    :param output_file: Name of the output file. By default, it will be called 'full_report.xlsx' and be saved in the
    output folder specified in the config directory.
    """
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

    logging.info('Started new run.')

    variant_obj = Variant("rs1354034", 3, 56815721, "T", "C")

    #variant_obj = Variant.from_rsid("rs10982456")
    #ldtrait(variant_obj)
    #variant_obj.results.ldtrait()
    eqtl_catalogue_LDblock_query_type_restricted_multitype(variant_obj, ['ge', 'microarray'])
    df = variant_obj.results.eqtl_cat()
    print(df)
    #print(df.drop_duplicates())
    df.to_csv("/home/antton/Desktop/eqtlcat.csv")
    #var = Variant("rs149143617", 1, 777870, "C", "G")
    #var = Variant("rs301816", 1, 8444998, "A", "G")
    #var = Variant.from_rsid('rs9272363')
    #print(var.get_gwas_phenotypes())
    #var.get_LDblock().to_csv(output_folder+'00-LDblock.csv')
    #generate_full_excel_file(var, output_folder+'01-FullReport')
    #print(plot_tokyo_ge(["ENSG00000000419", "ENSG00000004776", "ENSG00000284725", "ENSG00000284747"]))
    logging.info('Finished run.')
