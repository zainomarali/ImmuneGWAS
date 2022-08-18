from helpers.dbsnp import dbsnp_single_position_query


def test_dbsnp_single_position_query():
    assert dbsnp_single_position_query("NW_021160031.1", 17109) is not None  # Should return a row of dbsnp as a list
    assert type(dbsnp_single_position_query("NW_021160031.1", 17109)) == list

    assert dbsnp_single_position_query("NW_021160031.1", 17110) is None  # Should return None

