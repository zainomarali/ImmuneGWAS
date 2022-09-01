from ImmuneGWAS.resources.eqtl_cat import *
from ImmuneGWAS.Variant import Variant
import pytest


def test_get_eqtl_cat_file_list():
    file_list = get_eqtl_cat_file_list()
    assert file_list is not None
    assert type(file_list) == list


def test_get_studies_of_type():
    """
    Call get_studies_of_type() for every possible study type.
    """
    for study_type in ['ge', 'exon', 'tx', 'txrev', 'microarray']:
        output = get_studies_of_type(study_type)
        assert output is not None
        assert type(output) == list


def test_get_studies_of_type_wrong_type():
    """
    Call get_studies_of_type() with a wrong study type.
    Should raise a ValueError.
    """
    with pytest.raises(ValueError):
        get_studies_of_type('wrong_type')


def test_single_eqtl_catalogue_query_type_restricted():
    var = Variant.from_rsid('rs9272363')
    df = single_eqtl_catalogue_query_type_restricted(var.get_chrom(), var.get_pos(), 'ge')
    assert df is not None
    assert type(df) == pd.DataFrame


def test_eqtl_catalogue_LDblock_query_type_restricted_multitype():
    var = Variant.from_rsid('rs9272363')
    eqtl_catalogue_LDblock_query_type_restricted_multitype(var, ["ge"])
    df = var.results.eqtl_cat()
    assert df is not None
    assert type(df) == pd.DataFrame
