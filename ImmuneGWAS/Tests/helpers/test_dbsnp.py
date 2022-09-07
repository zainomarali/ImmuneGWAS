from ImmuneGWAS.helpers.dbsnp import *
from ImmuneGWAS.variant import Variant


def test_dbsnp_single_position_query():
    assert dbsnp_single_position_query(6, 569665) is not None  # Should return a row of dbsnp as a list
    assert type(dbsnp_single_position_query(6, 569665)) == list

    assert dbsnp_single_position_query(6, 569666) is None  # Should return None


def test_replace_rsid_column_with_dbsnp():
    df = pd.DataFrame({'rsid': ['rsWRONG1', 'rsWRONG2', 'rsWRONG3'],
                       'chrom': [1, 17, 22],
                       'pos': [10015, 43074471, 42129809]})

    replaced_df = replace_rsid_column_with_dbsnp(df, 'chrom', 'pos')
    correct_rsid = ['rs1570391706', 'rs1800744', 'rs28371704']

    assert replaced_df['rsid'].tolist() == correct_rsid


def test_replace_rsid_column_with_dbsnp_single_coord_column():
    """Same test nut setting 'single_position_column' to True instead.
    TODO: write this test"""
    assert False
