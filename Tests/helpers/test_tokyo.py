from helpers.tokyo import *
from Variant import Variant
import pandas as pd


def test_get_tokyo_eqtl_file_list():
    file_list = get_tokyo_eqtl_file_list()
    assert file_list is not None
    assert type(file_list) == list


def test_single_tokyo_eqtl_query():
    var = Variant("rs149143617", 1, 777870, "C", "G")  # Known to have an effect in multiple Tokyo eqtl cell types
    df = single_tokyo_eqtl_query(var.get_chrom(), var.get_pos())
    assert df is not None
    assert type(df) == pd.DataFrame
