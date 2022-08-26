
"""
File paths to all the different local resources on cbio3.
As the exact filepath to cbio3 will differ for every user, this functions will return the correct filepaths if the root
is provided.
"""


def get_paths(root):
    """
    Returns a dictionary where the keys correspond to the resource and the values are file paths.
    Some of these resources are folders with multiple files and some are single files.

    Keys:
        eqtl_cat: eQTL Catalogue, multiple eQTL studies
        dbsnp: dbSNP, for rsid harmonization etc
        eqtl_tokyo: significant SNPs from ImmunexUT
        ge_tokyo: gene expression from ImmunexUT
        eqtlgen_cis: cis-eQTLs from eQTLgen
        eqtlgen_trans: trans-eQTLs from eQTLgen
    """
    res_dict = {'eqtl_cat': root + "cbio3/data/eQTL_DB/",
                'dbsnp': root + "cbio3/data/dbSNP/GCF_000001405.39.gz",
                'eqtl_tokyo': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/eQTL_summarystats_with_alleles/",
                'ge_tokyo': root + "cbio3/projects/Zain_2021/ImmunexUT_GE/E-GEAD-397.processed/tpm/",
                'eqtlgen_cis': root + "cbio3/projects/Zain_2021/eQTLgen/data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded_REsorted_hg38.txt.gz",
                'eqtlgen_trans': root + "cbio3/projects/Zain_2021/eQTLgen/data/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded_sorted_hg38.txt.gz",
                'tokyo_alleles': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/metadata/full_allele_defs.txt"}

    return res_dict


def get_sumstats_path(root):
    return root + "cbio3/projects/antton/Immune_cell_GWAS/output/GWAS/combined_output/combined_hits_only_hg38.txt"
