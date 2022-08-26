import pandas as pd
import os
from helpers.getpaths import get_paths
from Variant import Variant
import config
import tabix

tokyo_eqtl_path = get_paths(config.cbio_root)['eqtl_tokyo']


def get_tokyo_eqtl_file_list() -> list:
    """
    Returns a list of the files in the Tokyo eQTL folder that have the structure "sorted_<name>_txt.gz" and a
    corresponding tabix index. There should be 28(?) files.
    Each cell type has its own file.

    :return: list of files in the Tokyo eQTL folder that have the structure "sorted_<name>_txt.gz" and .tbi file
    """
    all_files_in_dir = os.listdir(tokyo_eqtl_path)
    tabix_indexed_files = []
    for file in all_files_in_dir:
        if file.startswith('sorted_') and file.endswith(".txt.gz"):
            if os.path.exists(tokyo_eqtl_path + "/" + file + ".tbi"):
                tabix_indexed_files.append(file)
            else:
                print(f"WARNING: {file} is sorted and bgzipped, but no tabix index file exists.")

    return tabix_indexed_files


def single_tokyo_eqtl_query(chromosome, position, EA=None) -> pd.DataFrame:
    """
    Single Tokyo eQTL query.
    Using Tabix, go through the Tokyo eQTL files (there is one per cell type) and get the entries corresponding to the
    variant located at queried chromosome:position. Fill a dataframe with the matches and return it.

    :param chromosome: chromosome number
    :param position: position on the chromosome
    :param EA: Effect allele. Used to flip beta if necessary. If None, no flipping is done.
    :return: DataFrame with all the matches from the eQTL catalogue for the variant at the inputted position
    """
    tabix_indexed_files = get_tokyo_eqtl_file_list()
    list_of_celltype_match_dfs = []  # List of dataframes, one for each cell type in Tokyo that contains the variant
    for i, celltype_file in enumerate(tabix_indexed_files):
        tokyo_eqtl_file = tokyo_eqtl_path + "/" + celltype_file
        tokyo_eqtl_tabix_file = tabix.open(tokyo_eqtl_file)
        try:
            # Make the query. Returns iterator with all matches from tabix lookup.
            matches = tokyo_eqtl_tabix_file.querys(f"{chromosome}:{position}-{position}")
            # Convert iterator to pandas dataframe.
            match_list = [x for x in matches]  # Convert the generator to a list
            if match_list:
                header = ["Gene_id", "Gene_name", "CHR", "TSS_position", "Number_of_variants_cis", "Variant_ID", "OA",
                          "EA", "Variant_CHR", "Variant_position_start", "Variant_position_end", "Rank_of_association",
                          "Forward_nominal_P", "Forward_slope", "Backward_P", "Backward_slope", "int_chrom"]
                # TODO: hard-coding these feels wrong. Should at least double check this with file header.
                match_df = pd.DataFrame(match_list, columns=header)
                match_df["cell_type"] = celltype_file.split('_conditional_')[0].replace("sorted_", "")

                list_of_celltype_match_dfs.append(match_df)
        except tabix.TabixError:
            # TODO: better message. Pretty sure this exception is thrown when the variant is not in the file,
            #  so it's not really an error, just that the variant is not there (?)
            print("No output for file ", tokyo_eqtl_file)
            continue
    if list_of_celltype_match_dfs:
        concatenated_df = pd.concat(list_of_celltype_match_dfs)

        if EA:
            if EA == concatenated_df["EA"].iloc[0]:
                print(f"LDblock EA {EA} corresponds to Tokyo EA {concatenated_df['EA'].iloc[0]}")
            elif EA == concatenated_df["OA"].iloc[0]:
                print(
                    f"LDblock EA {EA} corresponds to Tokyo OA {concatenated_df['OA'].iloc[0]}. Flipping sign of the beta.")
                concatenated_df.Forward_slope = concatenated_df.Forward_slope * -1
            else:
                raise ValueError(f"WARNING: LDblock EA {EA} does not correspond to either Tokyo EA"
                                 f"{concatenated_df['EA'].iloc[0]} or Tokyo OA {concatenated_df['OA'].iloc[0]}")
        return concatenated_df  # Results from every cell type for a single SNP
    else:
        return pd.DataFrame()  # Return an empty dataframe if no matches were found.


def tokyo_eqtl_LDblock_query(variant_object: Variant) -> pd.DataFrame:
    """
    Query the Tokyo eQTL catalogue for each variant in the LD block of the Variant object.
    It calls :py:func:`single_tokyo_eqtl_query` for each variant in the LD block, stores each output dataframe in a
    list, and then concatenates the list into a single dataframe.

    :param variant_object: Variant object
    :return: DataFrame with all the matches from the eQTL catalogue for the LD block of the Variant object
    """
    # Lookup for the user-inputted variant, a.k.a. the lead variant, on its own
    lead_variant_df = single_tokyo_eqtl_query(variant_object.get_chrom(), variant_object.get_pos(),
                                              EA=variant_object.get_EA())
    tokyo_eqtl_matches_list = [lead_variant_df]  # List of dfs, one per variant in LD block, used for concatenation

    LDblock_df = variant_object.get_LDblock()
    if LDblock_df.empty:  # If the LD block is empty, return the lead variant dataframe alone.
        print("WARNING: The Variant object has no LDblock attribute. Returning lead variant only.")
        return lead_variant_df

    variant_positions_list_of_lists = []  # [[chromosome, position, EA], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos', 'EA',
             'RS_Number']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be more thorough.
        raise ValueError("LDblock dataframe has no 'chrom' or 'hg38_pos' columns.")
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:  # Iterate over the LDblock, make a query for each SNP
            # Create a DataFrame with all the matches for the variant at the position (list[0], list[1])
            variant_df = single_tokyo_eqtl_query(variant_pos_list[0], variant_pos_list[1], EA=variant_pos_list[2])
            if not variant_df.empty:
                print(  # Companion message for the print statement in single_tokyo_eqtl_query. Should be stored in log
                    f" for variant {variant_pos_list[3]} at position {variant_pos_list[1]} on chromosome {variant_pos_list[0]}:")

            tokyo_eqtl_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(tokyo_eqtl_matches_list)

    return concat_df


def tokyo_eqtl_to_summary_table(tokyo_eqtl_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the dataframe with all the matches from the eQTL catalogue for the LD block of the Variant object to a
    summary table.
    The summary table has the following columns:

    - gene_id: gene id
    - cell_type: cell type
    - p_value: p value
    - beta : effect size. Harmonized so that the sign corresponds to the EA
    - resource : database/resource where the data is coming from

    :param tokyo_eqtl_df: DataFrame outputted by the Tokyo lookup functions
    :return: Summary table with the above columns.
    """
    if tokyo_eqtl_df.empty:
        print("WARNING: No matches found in Tokyo eQTL dataset. Returning empty dataframe.")
        return pd.DataFrame()
    elif not {"Gene_id", "cell_type", "Forward_nominal_P", "Forward_slope"}.issubset(tokyo_eqtl_df.columns):
        raise ValueError("Input dataframe does not have the required columns. The required columns are: Gene_id, "
                         "cell_type, Forward_nominal_P, Forward_slope")

    df = tokyo_eqtl_df[["Gene_id", "cell_type", "Forward_nominal_P", "Forward_slope"]].copy()
    df.columns = ["gene_id", "cell_type", "p_value", "beta"]
    df["resource"] = "Tokyo"

    return df


def tokyo_eqtl_LDblock_query_simple(variant_object: Variant) -> pd.DataFrame:
    """
    Call :py:func:`tokyo_eqtl_LDblock_query` and :py:func:`tokyo_eqtl_to_summary_table` to get the summary table
    directly, without the full table. It is essentially the same as calling :py:func:`tokyo_eqtl_LDblock_query` but
    the output table only contains a subset of the columns.

    :param variant_object: Variant object
    """
    df = tokyo_eqtl_to_summary_table(tokyo_eqtl_LDblock_query(variant_object))

    return df
