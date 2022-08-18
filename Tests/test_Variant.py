import pytest
from Variant import Variant


def test_variant():
    var = Variant("rs1308164722", 1, 25231562, "T", "G")
    assert var.get_rsid() == "rs1308164722"
    assert var.get_chrom() == 1
    assert var.get_pos() == 25231562
    assert var.get_EA() == "T"
    assert var.get_OA() == "G"


def test_from_rsid():
    """
    Test the alternate 'from_rsid' constructor.
    """
    input_rsid = "rs9438875"
    var = Variant.from_rsid(input_rsid)

    assert var.get_rsid() == input_rsid
    assert var.get_chrom() == 1
    assert var.get_pos() == 24905071
    assert var.get_EA() == "T"
    assert var.get_OA() == "G"


def test_cross_reference_dbsnp():
    assert False


@pytest.fixture
def variant_object():
    return Variant("rs1308164722", 1, 25231562, "T", "G")


def test_get_rsid(variant_object):
    assert variant_object.get_rsid() == "rs1308164722"


def test_get_fullpos(variant_object):
    assert variant_object.get_fullpos() == (1, 25231562)


def test_get_pos(variant_object):
    assert variant_object.get_pos() == 25231562


def test_get_chrom(variant_object):
    assert variant_object.get_chrom() == 1


def test_get_EA(variant_object):
    assert variant_object.get_EA() == "T"


def test_get_OA(variant_object):
    assert variant_object.get_OA() == "G"
