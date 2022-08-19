from Variant import Variant
from helpers.dbsnp import dbsnp_single_position_query
from helpers.ldlink import ldproxy, ldtrait
from helpers.eqtl_cat import get_eqtl_cat_file_list

"""
This script exists in place of a main function for now to test new functionality.
"""
var = Variant.from_rsid("rs7292711")
print("Position is: ", var.get_pos())

df = ldproxy(var.get_rsid())
print("LDproxy query output: ",df)

print("dbsnp_single_position_query(6, 569665): ", dbsnp_single_position_query(6, 569665))
# var = Variant.from_rsid("rs943")  # This must fail

print("eQTL Catalogue file list: ", get_eqtl_cat_file_list())
