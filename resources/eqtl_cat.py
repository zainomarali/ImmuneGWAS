import os
import pandas as pd
import tabix
import logging

from helpers.getpaths import get_paths
from helpers.ensembl import get_gene_symbol
from Variant import Variant
import config


eqtl_cat_path = get_paths(config.cbio_root)[
                    'eqtl_cat'] + "/curated_crediblesets"  # Path to eQTL catalogue directory. Contains several files.


def get_eqtl_cat_file_list() -> list:
    """
    Returns a list of ALL the files in the eQTL catalogue folder.
    There are 4 files per study. Studies are divided into 5 types: 'ge', 'exon', 'tx', 'txrev' and 'microarray'.
    """
    return os.listdir(eqtl_cat_path)


def get_studies_of_type(study_type_key: str) -> list:
    """
    Returns a list of all the studies in eQTL catalogue that are of the specified type.

    :param study_type_key: key for the study type. Options are 'ge', 'exon', 'tx', 'txrev' and 'microarray'
    :return: list of all the studies of the requested type in the eQTL catalogue data
    """
    if study_type_key not in ['ge', 'exon', 'tx', 'txrev', 'microarray']:
        raise ValueError("Invalid study type specified. study_type_key must be one of 'ge', 'exon', 'tx', 'txrev' or "
                         "'microarray'")
    all_studies = get_eqtl_cat_file_list()
    studies_of_type = []
    for study in all_studies:
        if "_" + study_type_key + "_" in study and study.endswith("txt.gz"):
            if os.path.exists(eqtl_cat_path + "/" + study + ".tbi"):
                studies_of_type.append(study)
            else:
                print(f"WARNING: no tabix index file could be found for {study} .")
    return studies_of_type


def single_eqtl_catalogue_query_type_restricted(chromosome: int, position: int, study_type_key: str, EA: str = None) \
        -> pd.DataFrame:
    """
    Single eQTL catalogue query, type restricted.
    Go through the eQTL catalogue studies of the selected type and return the eQTLs that are significant for the variant
    located at chromosome:position.
    This function is called by the function  :py:func:`eqtl_catalogue_LDblock_query_type_restricted` once for every
    variant in the LD block of a Variant object.

    :param chromosome: chromosome number
    :param position: position on the chromosome
    :param study_type_key: key for the study type. Options are 'ge', 'exon', 'tx', 'txrev' and 'microarray'
    :param EA: Effect Allele. If specified, the alleles will be checked against the ref & alt in the eQTL catalogue and
    the sign of the effect size will be adjusted accordingly.
    :return: DataFrame with all the matches from the eQTL catalogue for the variant at the inputted position
    """

    if study_type_key not in ['ge', 'exon', 'tx', 'txrev', 'microarray']:
        raise ValueError("Invalid study type specified. study_type_key must be one of: 'ge', 'exon', 'tx', 'txrev' or "
                         "'microarray'")
    studies_of_requested_type_list = get_studies_of_type(study_type_key=study_type_key)
    list_of_study_match_dfs = []  # List of dataframes, one for each study that contains the variant
    for i, study in enumerate(studies_of_requested_type_list):
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

                # Convert gene ID to gene symbol
                df.insert(1, "gene_symbol", df["molecular_trait_id"].apply(get_gene_symbol))
                # Get study name, type and cell-type from file name
                study_name_split = study.split(f"_{study_type_key}_")
                df["study"] = study_name_split[0]  # Add the study name to the dataframe
                df["study_type"] = 'ge'  # Add the study type to the dataframe
                df["cell_type"] = study_name_split[1].split('.')[0]  # Add the cell type to the dataframe
                list_of_study_match_dfs.append(df)

        except tabix.TabixError:
            # TODO: better message. Pretty sure this exception is thrown when the variant is not in the file,
            #  so it's not really an error, just that the variant is not there
            print("Something went wrong when trying to access the file ", eqtl_cat_file, ". Skipping this study.")
            continue
    if list_of_study_match_dfs:
        concatenated_df = pd.concat(list_of_study_match_dfs)
        # EA check!:
        if EA:  # In eQTLcat ALT should always be the effect allele. (https://www.ebi.ac.uk/eqtl/Data_access/)
            if EA == concatenated_df["alt"].iloc[0]:
                pass
            elif concatenated_df["alt"].iloc[0][0] == concatenated_df["ref"].iloc[0][0]:  # Deletion/addition case
                # LDlink represents -/T while eQTL-cat does G/GT
                if EA == concatenated_df["alt"].iloc[0][1:]:  # Addition:
                    pass  # We compare T(EA) and GT[1:](alt)
                elif EA == '-' and concatenated_df["alt"].iloc[0] == concatenated_df["ref"].iloc[0][0]:  # Deletion
                    pass  # We check EA is '-' and compare G(alt) and GT[0](ref)
                else:
                    raise ValueError("Something went wrong when checking EA harmony with INSERTION/DELETION.\n"
                                     f"Variant located at: {chromosome}:{position} (hg38)\n"
                                     f"LDblock EA = {EA} ,\teQTL catalogue ALT = {concatenated_df['alt'].iloc[0]} ,\t"
                                     f"LDblock ref = {concatenated_df['ref'].iloc[0]}")
            else:
                raise ValueError("The effect allele is not the same as the ALT allele in the eQTL catalogue.\n"
                                 f"Variant located at: {chromosome}:{position} (hg38)\n"
                                 f"LDblock EA = {EA}, \teQTL catalogue ALT = {concatenated_df['alt'].iloc[0]} ,\t"
                                 f"LDblock REF = {concatenated_df['ref'].iloc[0]}")
        return concatenated_df
    else:
        return pd.DataFrame()  # Return an empty dataframe if no matches were found.


def eqtl_catalogue_LDblock_query_type_restricted(variant_object: Variant, study_type_key: str) -> pd.DataFrame:
    """
    Go through the eQTL catalogue studies of the specified type and return the eQTLs that are significant for every
    variant in LD block. This function calls 'single_eqtl_catalogue_query_type_restricted()' once for
    every variant in the LD block, and concatenates the resulting dataframes into a single one.

    :param variant_object: Variant object
    :param study_type_key: key for the study type. Options are 'ge', 'exon', 'tx', 'txrev' and 'microarray'
    :return: DataFrame with all the matches from the eQTL catalogue for every SNP in the LDblock
    """
    lead_variant_df = single_eqtl_catalogue_query_type_restricted(variant_object.get_chrom(), variant_object.get_pos(),
                                                                  study_type_key)
    eqtl_cat_matches_list = [lead_variant_df]  # List of dataframes, used to concat all the dfs together.
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos', 'EA']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be better. If the columns doesn't exist probably means dataframe is empty.
        print("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all the matches from eQTL cat studies for the variant at the position (list[0], list[1])
            try:
                variant_df = single_eqtl_catalogue_query_type_restricted(variant_pos_list[0], variant_pos_list[1],
                                                                         study_type_key, EA=variant_pos_list[2])
            except ValueError as err_msg:  # Avoid breaking if alleles don't match, but report it and don't add to list
                print("WARNING: ", err_msg)
                continue
            eqtl_cat_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(eqtl_cat_matches_list)

    return concat_df


def eqtl_catalogue_LDblock_query_type_restricted_multitype(variant_object: Variant,
                                                           input_study_list: list = None) -> pd.DataFrame:
    """
    Go through the eQTL catalogue studies of all types and return the eQTLs that are significant for every variant.
    :param variant_object: Variant object
    :param input_study_list: list of study types to query. Should be a list containing a slice of the list
    ['ge', 'exon', 'tx', 'txrev' and 'microarray']. If none is specified, all study types are queried.
    """
    logging.info(f"Querying eQTL catalogue files for full LD block of Variant {variant_object.rsid}, for study types: "
                 f"{', '.join(input_study_list)}")

    df_list = []
    if input_study_list is None:
        study_list = ['ge', 'exon', 'tx', 'txrev', 'microarray']
    else:
        study_list = list(set(input_study_list))  # Remove duplicates just in case
    for study_type_key in study_list:
        if study_type_key not in ['ge', 'exon', 'tx', 'txrev', 'microarray']:
            raise ValueError("Invalid study type specified. study_type_key must be one of 'ge', 'exon', 'tx', "
                             "'txrev' or 'microarray'")
        df_list.append(eqtl_catalogue_LDblock_query_type_restricted(variant_object, study_type_key))
    logging.info("Query to eQTL catalogue finished.")
    return pd.concat(df_list)


def eqtl_catalogue_to_summary_table(eqtl_cat_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the dataframe with all the matches from the eQTL catalogue for the LD block of the Variant object to a
    summary table.
    The summary table has the following columns:

    - gene_id: gene id
    - cell_type: cell type
    - p_value: p value NOTE: we call this p-value but we are using the 'pip'  # TODO: decide what to do with this
    - beta : effect size. Harmonized so that the sign corresponds to the EA
    - resource : database/resource where the data is coming from

    :param eqtl_cat_df: DataFrame outputted by the eQTL catalogue lookup functions
    :return: Summary table with the above columns.
    """
    if eqtl_cat_df.empty:
        print("WARNING: No matches found in eQTL catalogue. Returning empty dataframe.")
        return pd.DataFrame()
    elif not {"molecular_trait_id", "cell_type", "pip", "z"}.issubset(eqtl_cat_df.columns):
        raise ValueError("Dataframe does not have the required columns. Required columns are: "
                         "molecular_trait_id, cell_type, pip, z")

    df = eqtl_cat_df[["molecular_trait_id", "cell_type", "pip", "z"]].copy()
    df.columns = ["gene_id", "cell_type", "p_value", "beta"]
    df["resource"] = "eqtl_cat"

    # TODO: Harmonize the sign of the beta so that the sign corresponds to the EA

    return df


def eqtl_catalogue_LDblock_query_type_restricted_multitype_simple(variant_object: Variant,
                                                                  input_study_list: list = None) -> pd.DataFrame:
    """
    Call :py:func:`eqtl_catalogue_LDblock_query_type_restricted_multitype` and
    :py:func:`eqtl_catalogue_to_summary_table` to get the simplified summary table directly, without the full table.

    :param variant_object: Variant object
    :param input_study_list: list of study types to query. Should be a list containing a slice of the list
    ['ge', 'exon', 'tx', 'txrev' and 'microarray']. If none is specified, all study types are queried.
    """
    df = eqtl_catalogue_to_summary_table(eqtl_catalogue_LDblock_query_type_restricted_multitype(variant_object,
                                                                                                input_study_list))
    df.drop_duplicates(inplace=True)  # Remove duplicates

    return df
