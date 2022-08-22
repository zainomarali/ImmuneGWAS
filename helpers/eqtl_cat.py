from helpers.getpaths import get_paths
from Variant import Variant
import config
import os
import pandas as pd
import tabix


eqtl_cat_path = get_paths(config.cbio_root)['eqtl_cat']+"/curated_crediblesets"  # Path to eQTL catalogue directory. Contains several files.


def get_eqtl_cat_file_list() -> list:
    """
    Returns a list of ALL the files in the eQTL catalogue folder.
    There are 4 files per study.
    """
    return os.listdir(eqtl_cat_path)


def get_gene_expression_studies() -> list:
    """
    Returns a list of all the gene expression studies in the eQTL catalogue. These files should already have a
    corresponding tabix index file in the same folder.
    """
    all_studies = get_eqtl_cat_file_list()
    gene_expression_studies_list = []
    for study in all_studies:
        if "_ge_" in study and study.endswith("txt.gz"):
            gene_expression_studies_list.append(study)
    return gene_expression_studies_list


def single_eqtl_catalogue_gene_expression_query(chromosome, position) -> pd.DataFrame:
    """
    o through the eQTL catalogue gene expression studies and return the eQTLs that are significant for the variant at
    the inputed position.
    This function is called by the function 'eqtl_catalogue_gene_expression_LDblock_query' once for every variant in the LD block.

    :param chromosome: chromosome number
    :param position: position on the chromosome
    :return: DataFrame with all the matches from the eQTL catalogue for the variant at the inputed position
    """

    gene_expression_studies_list = get_gene_expression_studies()
    list_of_study_match_dfs = []  # List of dataframes, one for each study that contains the variant
    for i, study in enumerate(gene_expression_studies_list):
        eqtl_cat_file = eqtl_cat_path + "/" + study
        eqtl_cat_tabix_file = tabix.open(eqtl_cat_file)
        try:
            # Make the query. Returns iterator with all matches from tabix lookup.
            matches = eqtl_cat_tabix_file.querys(f"{chromosome}:{position}-{position}")
            match_list = [x for x in matches]  # Convert the generator to a list of matches for the queried position
            if match_list:
                headers = ["molecular_trait_id", "variant", "chromosome", "position", "ref", "alt", "cs_id", "cs_index",
                           "finemapped_region", "pip", "z", "cs_min_r2", "cs_avg_r2", "cs_size", "posterior_mean",
                           "posterior_sd", "cs_log10bf", "int_chrom"]
                # TODO: header could be a lookup from the bare file instead of hard coding.

                df = pd.DataFrame(match_list, columns=headers)  # DataFrame with all the matches from the study
                # TODO: instead of adding 'study' alone, add several columns (split the file name).
                #  add at least: study_name and cell_type.
                df["study"] = study  # Add the study name to the dataframe
                list_of_study_match_dfs.append(df)

        except tabix.TabixError:
            print("Something went wrong when trying to access the file ", eqtl_cat_file, ". Skipping this study.")
            continue
    if list_of_study_match_dfs:
        return pd.concat(list_of_study_match_dfs)
    else:
        return pd.DataFrame()  # Return an empty dataframe if no matches were found.


def eqtl_catalogue_gene_expression_LDblock_query(variant_object: Variant) -> pd.DataFrame:
    """
    Go through the eQTL catalogue gene expression studies and return the eQTLs that are significant for the variant and
    every other variant in its LD block. This function calls 'single_eqtl_catalogue_gene_expression_query()' once for
    every variant in the LD block, and concatenates the resulting dataframes into a single one.

    :param variant_object: Variant object
    :return: DataFrame with all the matches from the eQTL catalogue for every SNP in the LDblock
    """
    lead_variant_df = single_eqtl_catalogue_gene_expression_query(variant_object.get_chrom(), variant_object.get_pos())
    eqtl_cat_matches_list = [lead_variant_df]  # List of dataframes, used to concat all the dfs together.
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[['chrom', 'hg38_pos']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # This check could be better. If the columns don't exist probably means dataframe is empty.
        print("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all the matches from eQTL cat studies for the variant at the position (list[0], list[1])
            variant_df = single_eqtl_catalogue_gene_expression_query(variant_pos_list[0], variant_pos_list[1])
            eqtl_cat_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(eqtl_cat_matches_list)

    return concat_df
