import pytest
from ImmuneGWAS.helpers.ldlink import *

"""
Testing for the LDLink related functions. Note that these functions are quite slow to run because they need to connect
to the LDLink API.
We use the SNP rs16886165 because we know it has output for both LDproxy and LDtrait. This could in the future, however,
and a reason for the tests failing could be the the LDlink database itself has changed!
"""


@pytest.fixture(scope="module")
def ldproxy_df():
    return ldproxy('rs16886165')


@pytest.fixture(scope="module")
def ldtrait_df():
    return ldtrait('rs16886165')


def test_ldproxy(ldproxy_df) -> None:
    """
    Test that the dataframe is created and is not empty.
    """
    assert not ldproxy_df.empty


def test_ldproxy_columns(ldproxy_df) -> None:
    """
    Test that the dataframe has the necessary columns
    """
    assert "R2" in ldproxy_df.columns


def test_ldtrait(ldtrait_df) -> None:
    """
    Test that the dataframe is created and is not empty.
    """
    assert not ldtrait_df.empty


def test_ldtrait_columns(ldtrait_df) -> None:
    """
    Test that the dataframe has the necessary columns.
    """
    assert "R2" in ldtrait_df.columns


# TODO: Add tests for when the variant has not LDblock or output is faulty.
