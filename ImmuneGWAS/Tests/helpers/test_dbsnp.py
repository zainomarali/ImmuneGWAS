from ImmuneGWAS.helpers.dbsnp import *
from ImmuneGWAS.variant import Variant


def test_dbsnp_single_position_query():
    assert dbsnp_single_position_query(6, 569665) is not None  # Should return a row of dbsnp as a list
    assert type(dbsnp_single_position_query(6, 569665)) == list

    assert dbsnp_single_position_query(6, 569666) is None  # Should return None


def test_replace_rsid_column_with_dbsnp():
    df = pd.DataFrame({'rsid': ['rsWRONG1', 'rsWRONG2', 'rsWRONG3']})
    assert False
