"""
File paths to all the different local resources on cbio3.
As the exact filepath to cbio3 will differ for every user, this function will return the correct filepaths if the root
is provided.
Some of these resources are folders with multiple files and some are single files.

Returns a dictionary where the keys correspond to the resource and the values are file paths.

Keys:
    eqtl_cat: eQTL Catalogue, multiple eQTL studies
    dbsnp: dbSNP, for rsid harmonization etc
    eqtl_tokyo: significant SNPs from ImmunexUT
    ge_tokyo: gene expression from ImmunexUT
    eqtlgen_cis: cis-eQTLs from eQTLgen
    eqtlgen_trans: trans-eQTLs from eQTLgen
"""


def get_paths(root):
    res_dict = {'eqtl_cat': root + "cbio3/data/eQTL_DB/",
                'dbsnp': root + "cbio3/data/dbSNP/",
                'eqtl_tokyo': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/eQTL_summarystats/",
                'ge_tokyo': root + "cbio3/projects/Zain_2021/ImmunexUT_GE/E-GEAD-397.processed/tpm/",
                'eqtlgen_cis': root + "cbio3/projects/Zain_2021/eQTLgen/data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded_sorted.txt ",
                'eqtlgen_trans': root + "cbio3/projects/Zain_2021/eQTLgen/data/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt ",
                'tokyo_alleles': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/metadata/full_allele_defs.txt"}

    return res_dict


def get_sumstats_path(root):
    return root + "cbio3/projects/antton/Immune_cell_GWAS/output/GWAS/combined_output/combined_hits_only.txt"
