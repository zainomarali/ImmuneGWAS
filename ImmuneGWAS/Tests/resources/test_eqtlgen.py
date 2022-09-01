import ImmuneGWAS.Variant
from ImmuneGWAS.resources.eqtlgen import *

def test_single_eqtlgen_cis_query():
    assert False

def test_eqtlgen_cis_LDblock_query():

    var = ImmuneGWAS.Variant.Variant.from_rsid('rs9272363')
    eqtlgen_cis_LDblock_query(var)
    df = var.results.eqtlgen_cis()
    assert df is not None
    assert type(df) == pd.DataFrame
