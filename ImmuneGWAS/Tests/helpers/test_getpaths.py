from ImmuneGWAS.helpers.getpaths import *


def test_get_paths():
    """
    Test that the file paths are correct. Beware that the filepaths might be changed in cbio3, making this test useless.
    """
    p_dict = get_paths("/path/to/dir/")

    assert p_dict['eqtl_cat'] == "/path/to/dir/cbio3/data/eQTL_DB/"
    assert p_dict['dbsnp'] == "/path/to/dir/cbio3/data/dbSNP/GCF_000001405.39.gz"
    assert p_dict['eqtl_tokyo'] == "/path/to/dir/cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021" \
                                   "/eQTL_summarystats_with_alleles/"
    assert p_dict['ge_tokyo'] == "/path/to/dir/cbio3/projects/Zain_2021/ImmunexUT_GE/E-GEAD-397.processed" \
                                 "/antton_reprocessing/mean_table_TPM.txt"
    assert p_dict['eqtlgen_cis'] == "/path/to/dir/cbio3/projects/Zain_2021/eQTLgen/data/2019-12-11-cis-eQTLsFDR0.05" \
                                    "-ProbeLevel-CohortInfoRemoved-BonferroniAdded_REsorted_hg38.txt.gz"
    assert p_dict['eqtlgen_trans'] == "/path/to/dir/cbio3/projects/Zain_2021/eQTLgen/data/2018-09-04-trans-eQTLsFDR0" \
                                      ".05-CohortInfoRemoved-BonferroniAdded_sorted_hg38.txt.gz"


def test_get_sumstats_path():
    assert get_sumstats_path("/path/to/dir/") == "/path/to/dir/cbio3/projects/antton/Immune_cell_GWAS/data" \
                                                 "/hits_only_table_hg38.txt"
