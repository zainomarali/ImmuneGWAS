
"""
File paths to all the different local resources on cbio3.
As the exact filepath to cbio3 will differ for every user, this functions will return the correct filepaths if the root
is provided.
"""


def get_paths(root: str) -> dict:
    """
    Returns a dictionary where the keys correspond to the resource and the values are file paths.
    Some of these resources are folders with multiple files and some are single files.

    Keys:
        - eqtl_cat: eQTL Catalogue, multiple eQTL studies
        - dbsnp: dbSNP, for rsid harmonization etc
        - ensembl: Ensembl, list with aliases for gene id to gene symbol conversion
        - eqtl_tokyo: SNPs with significant eQTL signal from ImmunexUT
        - ge_tokyo: gene expression from ImmunexUT
        - eqtlgen_cis: cis-eQTLs from eQTLgen
        - eqtlgen_trans: trans-eQTLs from eQTLgen
        - tokyo_alleles: alleles corresponding to each SNP in eqtl_tokyo. NOTE: eqtl_tokyo was updated to include
          alleles itself

    :param root: path to directory where cbio3 is located in current computer. Should end in '/'
    """

    if root.endswith('cbio3/'):
        root = root[:-6]  # Remove the cbio3/ from the end of the path

    res_dict = {'eqtl_cat': root + "cbio3/data/eQTL_DB/",
                'dbsnp': root + "cbio3/data/dbSNP/GCF_000001405.39.gz",
                'ensembl': root + "cbio3/data/ensembl_biomart/gene_aliases.json",
                'eqtl_tokyo': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/eQTL_summarystats_with_alleles/",
                'ge_tokyo': root + "cbio3/projects/Zain_2021/ImmunexUT_GE/E-GEAD-397.processed/antton_reprocessing/mean_table_TPM.txt",
                'eqtlgen_cis': root + "cbio3/projects/Zain_2021/eQTLgen/data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded_REsorted_hg38.txt.gz",
                'eqtlgen_trans': root + "cbio3/projects/Zain_2021/eQTLgen/data/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded_sorted_hg38.txt.gz",
                'tokyo_alleles': root + "cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/metadata/full_allele_defs.txt",
                'immuneGWAS': root + "cbio3/projects/antton/Immune_cell_GWAS/data/hits_only_table_hg38.txt.gz"}

    return res_dict


def get_sumstats_path(root):
    """
    Returns the path to the combined Immune cell GWAS output summary statistics file.
    TODO: this file should be updated.
    """
    return root + "cbio3/projects/antton/Immune_cell_GWAS/output/GWAS/combined_output/combined_hits_only_hg38.txt"
