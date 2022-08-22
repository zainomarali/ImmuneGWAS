import os

from Variant import Variant
from helpers.dbsnp import dbsnp_single_position_query
from helpers.ldlink import ldproxy, ldtrait
from helpers.eqtl_cat import get_eqtl_cat_file_list, eqtl_catalogue_gene_expression_LDblock_query
import tabix

"""
This script exists in place of a main function for now to test new functionality.
"""

#var = Variant.from_rsid("rs476373")
#var = Variant.from_rsid("rs482811")
var = Variant.from_rsid("rs9272363")

print("Position is: ", var.get_fullpos())
print("LDblock columns: ", var.get_LDblock().columns)

print("eqtl_catalogue_gene_expression_LDblock_query(variant): ", eqtl_catalogue_gene_expression_LDblock_query(var))
