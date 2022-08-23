import os
import pandas as pd

from Variant import Variant
from helpers.dbsnp import dbsnp_single_position_query
from helpers.ldlink import ldproxy, ldtrait
from helpers.eqtl_cat import get_eqtl_cat_file_list, eqtl_catalogue_LDblock_query_type_restricted,\
    eqtl_catalogue_LDblock_query_type_restricted_all_types
import tabix

"""
This script exists in place of a main function for now to test new functionality.
"""

#var = Variant.from_rsid("rs476373")
#var = Variant.from_rsid("rs482811")
var = Variant.from_rsid("rs9272363")

print("Position is: ", var.get_fullpos())


with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print("eqtl_catalogue_LDblock_query_type_restricted_all_types(var): ",
          eqtl_catalogue_LDblock_query_type_restricted_all_types(var).head())
