from helpers.getpaths import get_paths
import config
import tabix

"""
Access the dbSNP file in cbio3.
This file was downloaded from the dbSNP database, release 155. It is in hg38 format.
The chromosomes are named by accession names instead of chromosome numbers. E.g. 'NC_000001.11' instead of 1.
For more information, check README file in cbio3/data/dbSNP/
"""

# TODO: add pytabix to requirements.txt
dbsnp_path = get_paths(config.cbio_root)['dbsnp']  # Path to dbSNP file

path_to_chr_RefSeq = '/'.join(dbsnp_path.split('/')[:-1] + ['chr_to_RefSeq.txt'])  # Path to file with chr to RefSeq mappings
with open(path_to_chr_RefSeq, 'r') as f:  # Load the file with the mappings into a dictionary
    chr_to_RefSeq_dict = {int(line.split('\t')[0]): line.split('\t')[1].rstrip() for line in f}


def dbsnp_single_position_query(SNP_chr : int, SNP_pos : int):
    """
    Query the dbSNP file for a single position.
    It converts the chromosome number to the RefSeq chromosome name before querying.

    :param SNP_chr: chromosome number or name
    :param SNP_pos: position on the chromosome
    :return: full row from the dbSNP database corresponding to that position as a list
    """
    tb = tabix.open(dbsnp_path)
    query_str = f"{chr_to_RefSeq_dict[SNP_chr]}:{SNP_pos}-{SNP_pos}"  # For example: NC_000006.12:17100-17100
    print("query_str ", query_str)
    matches = tb.querys(query_str)
    match_list = [x for x in matches]
    if not match_list:
        print(f"WARNING: No matches for {SNP_chr}:{SNP_pos}-{SNP_pos}")
        return None
    else:
        return match_list[0]  # Return the first match (there should only be one)
