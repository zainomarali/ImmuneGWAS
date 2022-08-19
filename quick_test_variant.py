from Variant import Variant
from helpers.dbsnp import dbsnp_single_position_query
from helpers.ldlink import ldproxy, ldtrait

"""
This script exists in place of a main function for now to test new functionality.
"""
var = Variant.from_rsid("rs7292711")
print(var.get_pos())

df = ldproxy(var.get_rsid())
print(df)

print("dbsnp_single_position_query(6, 569665): ", dbsnp_single_position_query(6, 569665))
# var = Variant.from_rsid("rs943")  # This must fail



