import tabix
import pandas as pd
import logging

import ImmuneGWAS.config as config
from ImmuneGWAS.helpers.getpaths import get_paths
import ImmuneGWAS.Variant as Variant

"""
This script contains functions used to query our own Immune cell GWAS hits file.
"""

immuneGWAS_path = get_paths(config.cbio_root)["immuneGWAS"]


def single_immuneGWAS_query(chromosome, position) -> pd.DataFrame:
    """
    Returns a dataframe with all the rows in Immune cell GWAS output table related to the SNP in the queried position.
    It uses Tabix to look for the inputted chromosome:position in the Immune cell GWAS hits file, and then stores all
    the matches into a dataframe.

    :param chromosome: Chromosome number
    :param position: Position in chromosome
    :return: DataFrame with all the matches from Immune cell GWAS hits file for variant at (chromosome, position)
    """

    try:
        immuneGWAS_tabix_file = tabix.open(immuneGWAS_path)
        # Make the tabix query. Returns iterator with all matches from tabix lookup.
        matches = immuneGWAS_tabix_file.querys(f"{chromosome}:{position}-{position}")
        match_list = [x for x in matches]  # Convert the generator to a list of matches for the queried position
        if match_list:
            headers = ["rsid", "rsid_short", "chromosome", "position", "A1", "A2", "pval", "beta", "tstat", "n",
                       "Phenotype"]
            # TODO: header could be a lookup from the bare file instead of hard coding, or at least cross check with it
            df = pd.DataFrame(match_list, columns=headers)  # DataFrame with all the matches
        else:  # If match list is empty, we return an empty dataframe
            df = pd.DataFrame()
    except tabix.TabixError:
        logging.exception(f"'tabix.TabixError' raised for tabix lookup for position {chromosome}:{position}-{position}"
                          f"in file {immuneGWAS_path}")
        df = pd.DataFrame()

    return df


def immuneGWAS_LDblock_query(variant_object: Variant) -> pd.DataFrame:
    """
    Query our immune cell GWAS results for every SNP in the LDblock of a Variant object.

    :param variant_object: Variant
    :return: DataFrame with all the matches from Immune cell GWAS hits file for all SNPs in the LDblock of the variant
    """
    logging.info(f"Querying immune-cell GWAS output file for full LD block of Variant {variant_object.rsid}")
    lead_variant_df = single_immuneGWAS_query(variant_object.get_chrom(), variant_object.get_pos())
    immuneGWAS_matches_list = [lead_variant_df]  # List of dataframes, used to concat all the dfs together.
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be better. If the columns doesn't exist probably means dataframe is empty.
        logging.warning("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all the matches from the GWAS file for the variant at the position (list[0], list[1])
            variant_df = single_immuneGWAS_query(variant_pos_list[0], variant_pos_list[1])
            immuneGWAS_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(immuneGWAS_matches_list)
    logging.info(f"Query to immuneGWAS output file complete.")

    return concat_df
