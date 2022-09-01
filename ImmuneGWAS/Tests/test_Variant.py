import pytest
from ImmuneGWAS.Variant import Variant
import pandas as pd


def test_variant():
    var = Variant("rs7292711", 22, 22716968, "A", "G")
    assert var.get_rsid() == "rs7292711"
    assert var.get_chrom() == 22
    assert var.get_pos() == 22716968
    assert var.get_EA() == "A"
    assert var.get_OA() == "G"


def test_from_rsid():
    """
    Test the alternate 'from_rsid' constructor. This should work as well as specifying the data manually, provided we
    have the variant in our sumstats file.
    """
    input_rsid = "rs7292711"
    var = Variant.from_rsid(input_rsid)

    assert var.get_rsid() == input_rsid
    assert var.get_chrom() == 22
    assert var.get_pos() == 22716968
    assert var.get_EA() == "A"
    assert var.get_OA() == "G"


def test_from_rsid_not_in_sumstats():
    """
    Test the alternate 'from_rsid' constructor when the given rsID is not in the sumstats file.
    Should raise a ValueError.
    """
    with pytest.raises(ValueError):
        Variant.from_rsid("rs943")


@pytest.fixture
def variant_object():
    """
    Create a variant object for testing all the getters.
    """
    return Variant("rs7292711", 22, 22716968, "A", "G")


def test_get_rsid(variant_object):
    assert variant_object.get_rsid() == "rs7292711"


def test_get_fullpos(variant_object):
    assert variant_object.get_fullpos() == (22, 22716968)


def test_get_pos(variant_object):
    assert variant_object.get_pos() == 22716968


def test_get_chrom(variant_object):
    assert variant_object.get_chrom() == 22


def test_get_EA(variant_object):
    assert variant_object.get_EA() == "A"


def test_get_OA(variant_object):
    assert variant_object.get_OA() == "G"


def test_get_ldblock():
    var = Variant("rs7292711", 22, 22716968, "G", "A")
    df = var.get_LDblock()
    assert df is not None
    assert type(df) == pd.DataFrame


def test_set_ldblock():
    var = Variant("rs7292711", 22, 22716968, "G", "A")
    var.set_LDblock()
    assert var.get_LDblock() is not None


def test_set_ldblock_not_found():
    """
    If the rsID is not included in LDlink, an exception should be raised.
    """
    with pytest.raises(ValueError):
        Variant.from_rsid("rs1308164722")
