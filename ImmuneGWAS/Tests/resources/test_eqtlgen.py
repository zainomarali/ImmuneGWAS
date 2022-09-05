import ImmuneGWAS.variant
from ImmuneGWAS.resources.eqtlgen import *


def test_single_eqtlgen_cis_query():
    var = Variant.from_rsid('rs9272363')
    df = single_eqtlgen_cis_query(var.get_chrom(), var.get_pos())

    assert df is not None
    assert type(df) == pd.DataFrame


def test_eqtlgen_cis_LDblock_query():
    var = ImmuneGWAS.variant.Variant.from_rsid('rs9272363')
    eqtlgen_cis_LDblock_query(var)
    df = var.results.eqtlgen_cis()

    assert df is not None
    assert type(df) == pd.DataFrame
