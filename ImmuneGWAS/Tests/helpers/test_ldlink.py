import pytest
from ImmuneGWAS.helpers.ldlink import *
from ImmuneGWAS.variant import Variant

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
    # 'rs624896' has LDtrait output and is in dbSNP
    variant_obj = Variant("rs624896", 5, 114520357, "A", "G")  # Create new Variant object
    ldtrait(variant_obj)  # Update the .results object
    return variant_obj.results.get_ldtrait()  # This is a pandas dataframe


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


def test_ldtrait_missing_variant():
    """
    Test that the ldtrait function returns an empty dataframe when the variant is not found in the LDtrait database.
    """
    variant_obj = Variant("rs149143617", 1, 777870, "C", "G")  # This variant should not have any matches in LDtrait
    ldtrait(variant_obj)
    assert variant_obj.results.get_ldtrait().empty
