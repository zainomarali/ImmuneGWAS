import Variant


# TODO: When the code gets moved out of 'prototyping', this test should load the correct class.

def test_from_rsid():
    input_rsid = "rs9438875"
    var = Variant.from_rsid(input_rsid)

    assert var.get_rsid() == input_rsid
    assert var.get_chrom() == "1"
    assert var.get_pos() == "25231562"
    assert var.get_ea() == "T"
    assert var.get_oa() == "G"


def test_get_rsid():
    assert False


def test_get_fullpos():
    assert False


def test_get_pos():
    assert False


def test_get_chrom():
    assert False


def test_get_ea():
    assert False


def test_get_oa():
    assert False
