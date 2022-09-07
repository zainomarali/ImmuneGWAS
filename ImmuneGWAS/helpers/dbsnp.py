import pandas as pd
import tabix
import logging

from ImmuneGWAS.helpers.getpaths import get_paths
import ImmuneGWAS.config as config

"""
Access the dbSNP file in cbio3.
This file was downloaded from the dbSNP database, release 155. It is in hg38 format.
The chromosomes are named by accession names instead of chromosome numbers. E.g. 'NC_000001.11' instead of 1.
For more information, check README file in cbio3/data/dbSNP/
"""

dbsnp_path = get_paths(config.cbio_root)['dbsnp']  # Path to dbSNP file

path_to_chr_RefSeq = '/'.join(
    dbsnp_path.split('/')[:-1] + ['chr_to_RefSeq.txt'])  # Path to file with chr to RefSeq mappings
with open(path_to_chr_RefSeq, 'r') as f:  # Load the file with the mappings into a dictionary
    chr_to_RefSeq_dict = {int(line.split('\t')[0]): line.split('\t')[1].rstrip() for line in f}


def dbsnp_single_position_query(SNP_chr: int, SNP_pos: int):
    """
    Query the dbSNP file for a single position.
    It converts the chromosome number to the RefSeq chromosome name before querying.

    :param SNP_chr: chromosome number or name
    :param SNP_pos: position on the chromosome
    :return: full row from the dbSNP database corresponding to that position as a list
    """
    tb = tabix.open(dbsnp_path)
    query_str = f"{chr_to_RefSeq_dict[SNP_chr]}:{SNP_pos}-{SNP_pos}"  # For example: NC_000006.12:17100-17100
    matches = tb.querys(query_str)
    # The columns in the dbSNP file are: CHROM POS ID REF ALT QUAL FILTER INFO
    match_list = [x for x in matches]  # Convert the generator to a list
    if not match_list:
        logging.warning(f"No matches for {SNP_chr}:{SNP_pos}-{SNP_pos}")
        return None
    else:
        return match_list  # This will be a list with a variable number of elements.


def replace_rsid_column_with_dbsnp(df: pd.DataFrame, rsid_col_name: str, chr_col_name: str, pos_col_name: str) \
        -> pd.DataFrame:
    """
    For a df, replace its column with rsIDs with a new column obtained by getting the rsid listed dbsnp for that
    position.

    :param df: dataframe with a column with rsIDs, a column with chromosome numbers, and a column with positions
    :param rsid_col_name: name of the column with rsIDs
    :param chr_col_name: name of the column with chromosome numbers
    :param pos_col_name: name of the column with positions
    :return: dataframe with the same columns as the input, but with a new column with rsIDs obtained from dbsnp
    """
    df['rsid'] = df.apply(lambda row: dbsnp_single_position_query(row[chr_col_name], row[pos_col_name])[0][2], axis=1)

    return df
