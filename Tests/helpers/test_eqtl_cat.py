from helpers.eqtl_cat import *


def test_get_eqtl_cat_file_list():
    file_list = get_eqtl_cat_file_list()
    assert file_list is not None
    assert type(file_list) == list
