from Variant import Variant
from helpers.dbsnp import dbsnp_single_position_query

"""
This script exists in place of a main function for now to test new functionality.
"""

var = Variant.from_rsid("rs943")
print(var.get_pos())

print(dbsnp_single_position_query(6, 569665))



