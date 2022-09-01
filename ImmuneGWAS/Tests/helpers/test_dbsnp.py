from ImmuneGWAS.helpers.dbsnp import dbsnp_single_position_query


def test_dbsnp_single_position_query():
    assert dbsnp_single_position_query(6, 569665) is not None  # Should return a row of dbsnp as a list
    assert type(dbsnp_single_position_query(6, 569665)) == list

    assert dbsnp_single_position_query(6, 569666) is None  # Should return None

