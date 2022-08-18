from helpers.getpaths import get_paths
import config
import tabix

# TODO: add pytabix to requirements.txt
dbsnp_path = get_paths(config.cbio_root)['dbsnp']  # Path to dbSNP file


def dbsnp_single_position_query(SNP_chr, SNP_pos):
    """
    Query the dbSNP file for a single position. Returns full row from the dbSNP database corresponding to that position.
    :param SNP_chr: chromosome number or name
    :param SNP_pos: position on the chromosome
    """
    tb = tabix.open(dbsnp_path)
    matches = tb.querys(f"{SNP_chr}:{SNP_pos}-{SNP_pos}")  # For example: 1:17100-17100
    match_list = [x for x in matches]
    if not match_list:
        print(f"WARNING: No matches for {SNP_chr}:{SNP_pos}-{SNP_pos}")
        return None
    else:
        return match_list[0]  # Return the first match (there should only be one)
