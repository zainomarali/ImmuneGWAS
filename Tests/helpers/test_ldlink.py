from helpers.ldlink import *


def test_ldproxy() -> None:
    """
    Test that the dataframe is created and is not empty.
    """
    df = ldproxy('rs333')
    assert not df.empty


def test_ldtrait():
    """
    Test if LDtrait works.
    """
    df = ldproxy('rs333')
    assert not df.empty
