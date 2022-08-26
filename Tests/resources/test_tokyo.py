import pytest
from resources.tokyo import *
from Variant import Variant
import pandas as pd


def test_get_tokyo_eqtl_file_list():
    file_list = get_tokyo_eqtl_file_list()
    assert file_list is not None
    assert type(file_list) == list


@pytest.fixture
def variant_object():
    """
    Create a variant object that can be reused in multiple tests.
    """
    return Variant("rs149143617", 1, 777870, "C", "G")  # Known to have an effect in multiple Tokyo eqtl cell types


def test_single_tokyo_eqtl_query(variant_object):
    var = variant_object
    df = single_tokyo_eqtl_query(var.get_chrom(), var.get_pos())
    assert df is not None
    assert type(df) == pd.DataFrame


def test_single_tokyo_eqtl_query_EA_check():
    """
    Test that the EA check in the single_tokyo_eqtl_query() function works. If the EA of the query does not match either
    the EA or OA of the dataset, an error should be raised.
    """
    var = Variant("rs149143617", 1, 777870, "A", "G")  # EA in Tokyo is C, not A
    with pytest.raises(ValueError):
        single_tokyo_eqtl_query(var.get_chrom(), var.get_pos(), EA=var.get_EA())


def test_tokyo_eqtl_ldblock_query(variant_object):
    var = variant_object
    df = tokyo_eqtl_LDblock_query(var)
    assert df is not None
    assert type(df) == pd.DataFrame


def test_tokyo_eqtl_ldblock_query_missing_ldblock(variant_object):
    """
    Test the tokyo_eqtl_LDblock_query() function when the LDblock is missing. It should still work and return a
    DataFrame with the output of the lead variant alone.
    """
    var = variant_object
    var.set_LDblock(False, pd.DataFrame())  # Set the LDblock to be missing
    df = tokyo_eqtl_LDblock_query(var)
    assert df is not None
    assert type(df) == pd.DataFrame


def test_tokyo_eqtl_ldblock_query_corrupted_ldblock(variant_object):
    """
    Test the tokyo_eqtl_LDblock_query() function when the LDblock is corrupted. It should still work and return a
    DataFrame with the output of the lead variant alone.
    """
    var = variant_object
    # Set an incorrect LDblock dataframe that lacks the necessary columns:
    var.set_LDblock(False, pd.DataFrame(columns=['wrong', 'column', 'names'], data=[[1, 2, 3]]))
    with pytest.raises(ValueError):
        tokyo_eqtl_LDblock_query(var)
