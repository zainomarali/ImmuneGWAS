from helpers.ldlink import *

"""
Testing for the LDLink related functions. Note that these functions are quite slow to run because they need to connect
to the LDLink API.
"""


def test_ldproxy() -> None:
    """
    Test that the dataframe is created and is not empty.
    """
    df = ldproxy('rs16886165')
    assert not df.empty


def test_ldtrait() -> None:
    """
    Test that the dataframe is created and is not empty.
    """
    df = ldtrait('rs16886165')
    assert not df.empty
