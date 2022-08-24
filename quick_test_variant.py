import os
import pandas as pd

from Variant import Variant
from helpers.eqtl_cat import eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output
from helpers.eqtlgen import eqtlgen_cis_LDblock_query, eqtlgen_cis_LDblock_query_formatted_output
from helpers.tokyo import tokyo_eqtl_LDblock_query_formatted_output

"""
This script exists in place of a main function for now to test new functionality.
"""

#var = Variant.from_rsid("rs476373")
#var = Variant.from_rsid("rs482811")
var = Variant.from_rsid("rs9272363")
#var = Variant("rs149143617", 1, 777870, "C", "G")

print("Position is: ", var.get_fullpos())
df1 = eqtlgen_cis_LDblock_query_formatted_output(var)
df2 = tokyo_eqtl_LDblock_query_formatted_output(var)
df3 = eqtl_catalogue_LDblock_query_type_restricted_multiple_types_formatted_output(var, ["ge", "microarray"])

df = pd.concat([df1, df2, df3])
df.to_csv("/home/antton/Desktop/test.csv", index=False)
#eqtl_catalogue_LDblock_query_type_restricted_multiple_types(var, ["ge", "microarray"])
#eqtlgen_cis_LDblock_query(var)
#tokyo_eqtl_LDblock_query_formatted_output(var)

#.to_csv("/home/antton/Desktop/test.csv", index=False)
