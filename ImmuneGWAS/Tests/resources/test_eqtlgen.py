from ImmuneGWAS.variant import Variant
from ImmuneGWAS.resources.eqtlgen import *


# Tests for CIS-eQTL related functions
def test_single_eqtlgen_cis_query():
    var = Variant.from_rsid('rs9272363')
    df = single_eqtlgen_cis_query(var.get_chrom(), var.get_pos())

    assert df is not None
    assert type(df) == pd.DataFrame


def test_single_eqtlgen_cis_query_EA_check():
    """Make sure that the EA check works correctly. The Zscore value should be 'flipped' (multiplied by -1) if the
    EA corresponds to the OtherAllele instead"""

    # Case where EA and AssessedAllele are the same
    var1 = Variant("rs1354034", 3, 56815721, 'T', 'C')
    df1 = single_eqtlgen_cis_query(var1.get_chrom(), var1.get_pos(), 'T')

    # Case where EA corresponds to OtherAllele instead
    var2 = Variant("rs1354034", 3, 56815721, 'T', 'C')
    df2 = single_eqtlgen_cis_query(var2.get_chrom(), var2.get_pos(), 'C')

    assert df1.Zscore.astype(float).equals(df2.Zscore.astype(float) * -1)


def test_eqtlgen_cis_LDblock_query():
    var = Variant.from_rsid('rs9272363')
    eqtlgen_cis_LDblock_query(var)
    df = var.results.eqtlgen_cis()

    assert df is not None
    assert type(df) == pd.DataFrame


# Tests for TRANS-eQTL related functions
def test_single_eqtlgen_trans_query():
    var = Variant.from_rsid('rs9272363')
    df = single_eqtlgen_trans_query(var.get_chrom(), var.get_pos())

    assert df is not None
    assert type(df) == pd.DataFrame


def test_single_eqtlgen_trans_query_EA_check():
    """Make sure that the EA check works correctly. The Zscore value should be 'flipped' (multiplied by -1) if the
    EA corresponds to the OtherAllele instead"""

    # Case where EA and AssessedAllele are the same
    var1 = Variant("rs1354034", 3, 56815721, 'T', 'C')
    df1 = single_eqtlgen_trans_query(var1.get_chrom(), var1.get_pos(), 'T')

    # Case where EA corresponds to OtherAllele instead
    var2 = Variant("rs1354034", 3, 56815721, 'T', 'C')
    df2 = single_eqtlgen_trans_query(var2.get_chrom(), var2.get_pos(), 'C')

    assert df1.Zscore.astype(float).equals(df2.Zscore.astype(float) * -1)


def test_eqtlgen_trans_LDblock_query():

    var = Variant("rs1354034", 3, 56815721, 'T', 'C')
    eqtlgen_trans_LDblock_query(var)
    df = var.results.eqtlgen_trans()

    assert df is not None
    assert type(df) == pd.DataFrame
