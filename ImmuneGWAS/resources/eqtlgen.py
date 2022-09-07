import pandas as pd
import logging

from ImmuneGWAS.helpers.getpaths import get_paths
from ImmuneGWAS.variant import Variant
import ImmuneGWAS.config as config
import tabix

eqtlgen_cis_path = get_paths(config.cbio_root)['eqtlgen_cis']
eqtlgen_trans_path = get_paths(config.cbio_root)['eqtlgen_trans']


def single_eqtlgen_cis_query(chromosome: int, position: int, EA: str = None) -> pd.DataFrame:
    """
    Single eQTLgen cis-eQTL query. Uses Tabix to look for the inputted chromosome:position in the eqtlgen cis-eQTL file.
    If EA is given, makes sure that the given EA is the same as the AssessedAllele in the eqtlgen cis-eQTL file.
    If it is not, the Zscore value is flipped.
    :param chromosome: Chromosome number
    :param position: Position in chromosome
    :param EA: Effect allele. If given, it checks that the AssessedAllele in the eqtlgen cis-eQTL file is the same as EA
    :return: DataFrame with all the matches from eQTLgen for the variant at the position (chromosome, position)
    """
    try:
        eqtlgen_tabix_file = tabix.open(eqtlgen_cis_path)
        # Make the tabix query. Returns iterator with all matches from tabix lookup.
        matches = eqtlgen_tabix_file.querys(f"{chromosome}:{position}-{position}")
        match_list = [x for x in matches]  # Convert the generator to a list of matches for the queried position
        if match_list:
            headers = ["Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene",
                       "GeneSymbol", "GeneChr", "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP"]
            # TODO: header could be a lookup from the bare file instead of hard coding. Would also be more robust.
            df = pd.DataFrame(match_list, columns=headers)  # DataFrame with all the matches
        else:  # If match list is empty, we return an empty dataframe
            df = pd.DataFrame()
    except tabix.TabixError:
        logging.exception(f"'tabix.TabixError' raised for tabix lookup for position {chromosome}:{position}-{position}"
                          f"in file {eqtlgen_cis_path}")
        df = pd.DataFrame()
    if EA:
        if EA == df.iloc[0]['AssessedAllele']:
            pass
        elif EA == df.iloc[0]['OtherAllele']:  # If it matches the OtherAllele instead, we flip the Zscore
            df['Zscore'] = -(df['Zscore'].astype(float))
        else:
            raise ValueError(f"ERROR: LDblock EA {EA} does not correspond to either eQTLGen AssessedAllele "
                             f"{df.iloc[0]['AssessedAllele']} or OtherAllele {df.iloc[0]['OtherAllele']}")

    return df


def eqtlgen_cis_LDblock_query(variant_object: Variant):
    """
    Query the eQTLgen cis-eQTL file for the entire LD block stored in the variant object.
    """
    logging.info(f"Querying eQTLgen cis-eQTL file for full LD block of Variant {variant_object.rsid}")
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position, EA], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos', 'EA']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be better. If the columns doesn't exist probably means dataframe is empty.
        logging.warning("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")
    eqtlgen_cis_matches_list = []  # List of dataframes, used to concat all the dfs together.
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all the matches from eQTL cat studies for the variant at the position (list[0], list[1])
            variant_df = single_eqtlgen_cis_query(variant_pos_list[0], variant_pos_list[1], variant_pos_list[2])
            eqtlgen_cis_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(eqtlgen_cis_matches_list)
    logging.info(f"Query to eQTLgen cis-eQTL file complete.")
    variant_object.results.set_eqtlgen_cis_df(concat_df)
    return


def single_eqtlgen_trans_query(chromosome: int, position: int, EA: str = None) -> pd.DataFrame:
    """
    Single eQTLgen trans-eQTL query. It uses tabix to look up the given position and extract trans-eQTL information.
    If EA is given, it ensures it is the same as the AssessedAllele in the table.

    :param chromosome: Chromosome
    :param position: Position
    :param EA: Effect Allele to check for.
    :return: DataFrame with all matches for given position
    """
    try:
        eqtlgen_tabix_file = tabix.open(eqtlgen_trans_path)
        # Make the query. Returns iterator with all matches from tabix lookup.
        matches = eqtlgen_tabix_file.querys(f"{chromosome}:{position}-{position}")
        match_list = [x for x in matches]  # Convert the generator to a list of matches for the queried position
        if match_list:
            headers = ["Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene",
                       "GeneSymbol", "GeneChr", "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP"]
            # TODO: header could be a lookup from the bare file instead of hard coding. Would also be more robust.
            df = pd.DataFrame(match_list, columns=headers)  # DataFrame with all the matches
        else:
            df = pd.DataFrame()
    except tabix.TabixError:
        logging.exception(f"'tabix.TabixError' raised for tabix lookup for position {chromosome}:{position}-{position}"
                          f"in file {eqtlgen_cis_path}")
        df = pd.DataFrame()
    if EA:
        if EA == df.iloc[0]['AssessedAllele']:
            pass
        elif EA == df.iloc[0]['OtherAllele']:
            df['Zscore'] = -(df['Zscore'].astype(float))
        else:
            raise ValueError(f"ERROR: LDblock EA {EA} does not correspond to either eQTLGen AssessedAllele "
                             f"{df.iloc[0]['AssessedAllele']} or OtherAllele {df.iloc[0]['OtherAllele']}")

    return df


def eqtlgen_trans_LDblock_query(variant_object: Variant):
    """
    TRANS-eQTL lookup for entire LDblock.
    Query the eQTLgen trans-eQTL file for the entire LD block stored in the variant object.

    :param variant_object: Variant object
    """
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position, EA], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos', 'EA']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be better. If the columns doesn't exist probably means dataframe is empty.
        logging.warning("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")

    eqtlgen_trans_matches_list = []  # List of dataframes, used to concat all the dfs together.
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all matches from eQTLGen (trans) for the variant at the position (list[0], list[1]):
            variant_df = single_eqtlgen_trans_query(variant_pos_list[0], variant_pos_list[1], variant_pos_list[2])
            eqtlgen_trans_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(eqtlgen_trans_matches_list)
    variant_object.results.set_eqtlgen_trans_df(concat_df)
    return


def eqtlgen_cis_to_summary_table(eqtlgen_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the dataframe with all the matches from the eQTLGen lookup for the LD block of the Variant object to a
    summary table.
    The summary table should have the following columns:

    - gene_id: gene id
    - cell_type: cell type
    - p_value: p value
    - beta : effect size. Harmonized so that the sign corresponds to the EA
      NOTE: we call this 'beta' but the data we are using are actually z-scores! TODO : decide if we want to change this
    - resource : database/resource where the data is coming from


    :param eqtlgen_df: DataFrame outputted by the eQTLGen lookup functions
    :return: Summary table with the above columns.
    """
    if eqtlgen_df.empty:
        logging.warning("No matches found for eqtlgen_cis dataset. Returning empty dataframe.")
        return pd.DataFrame()
    elif not {"Gene", "Pvalue", "Zscore"}.issubset(eqtlgen_df.columns):
        raise ValueError("Dataframe does not have required columns.")

    df = eqtlgen_df[["Gene", "Pvalue", "Zscore"]].copy()
    df.columns = ["gene_id", "p_value", "beta"]
    df.insert(1, "cell_type", "whole_blood")  # Manually insert cell type, since all of eqtlgen is whole blood
    df["resource"] = "eqtlgen_cis"

    # TODO: Harmonize the sign of the beta so that the sign corresponds to the EA

    return df


def eqtlgen_cis_LDblock_query_simple(variant_object: Variant) -> pd.DataFrame:
    """
    Call :py:func:`eqtlgen_cis_LDblock_query` and :py:func:`eqtlgen_cis_to_summary_table` to get the summary table
    directly, without the full table.

    :param variant_object: Variant object
    """
    df = eqtlgen_cis_to_summary_table(eqtlgen_cis_LDblock_query(variant_object))

    return df
