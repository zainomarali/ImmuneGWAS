import os
import pandas as pd

from Variant import Variant
from helpers.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multiple_types
from helpers.eqtlgen import eqtlgen_cis_LDblock_query

"""
This script exists in place of a main function for now to test new functionality.
"""

#var = Variant.from_rsid("rs476373")
#var = Variant.from_rsid("rs482811")
var = Variant.from_rsid("rs9272363")

print("Position is: ", var.get_fullpos())

eqtlgen_cis_LDblock_query(var).to_csv("/home/antton/Desktop/test.csv")
