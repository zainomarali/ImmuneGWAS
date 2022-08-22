from helpers.getpaths import get_paths
import config
import os
import tabix


eqtl_cat_path = get_paths(config.cbio_root)['eqtl_cat']+"/curated_crediblesets"  # Path to eQTL catalogue directory. Contains several files.


def get_eqtl_cat_file_list() -> list:
    """
    Returns a list of all the files in the eQTL catalogue folder.
    """
    return os.listdir(eqtl_cat_path)


def get_gene_expression_studies() -> list:
    """
    Returns a list of all the gene expression studies in the eQTL catalogue.
    """
    all_studies = get_eqtl_cat_file_list()
    gene_expression_studies_list = []
    for study in all_studies:
        if "_ge_" in study:
            gene_expression_studies_list.append(study)
    return gene_expression_studies_list

# TODO: start lookup with the regular 'gene expression' studies (tagged as '_ge_' in the file names)
