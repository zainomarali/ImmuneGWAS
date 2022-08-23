import pandas as pd
from helpers.getpaths import get_paths
from Variant import Variant
import config
import tabix

eqtlgen_cis_path = get_paths(config.cbio_root)['eqtlgen_cis']


def single_eqtlgen_cis_query(chromosome: int, position: int) -> pd.DataFrame:
    """
    Single eQTLgen cis-eQTL query.
    """
    try:
        eqtlgen_tabix_file = tabix.open(eqtlgen_cis_path)
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
        print(f"No output for tabix lookup for position {chromosome}:{position}-{position} in file {eqtlgen_cis_path}")
        df = pd.DataFrame()

    return df


def eqtlgen_cis_LDblock_query(variant_object: Variant):
    """
    Query the eQTLgen cis-eQTL file for the entire LD block stored in the variant object.
    """
    lead_variant_df = single_eqtlgen_cis_query(variant_object.get_chrom(), variant_object.get_pos())
    eqtlgen_cis_matches_list = [lead_variant_df]  # List of dataframes, used to concat all the dfs together.
    LDblock_df = variant_object.get_LDblock()
    variant_positions_list_of_lists = []  # [[chromosome, position], ...]. We'll iterate over this list later
    if 'chrom' in LDblock_df.columns.to_list() and 'hg38_pos' in LDblock_df.columns.to_list():
        variant_positions_list_of_lists = LDblock_df[
            ['chrom', 'hg38_pos']].values.tolist()  # List of lists with [chr, position] for each variant
    else:  # TODO: This check could be better. If the columns doesn't exist probably means dataframe is empty.
        print("No keys 'chrom', 'hg38_pos' in LDblock_df.columns. Returning lead variant only.")
    if variant_positions_list_of_lists:
        for variant_pos_list in variant_positions_list_of_lists:
            # DataFrame with all the matches from eQTL cat studies for the variant at the position (list[0], list[1])
            variant_df = single_eqtlgen_cis_query(variant_pos_list[0], variant_pos_list[1])
            eqtlgen_cis_matches_list.append(variant_df)  # Add the dataframe to the list of dataframes
    concat_df = pd.concat(eqtlgen_cis_matches_list)

    return concat_df

# TODO: write lookup functions for eQTLgen Trans-eQTLs as well!

