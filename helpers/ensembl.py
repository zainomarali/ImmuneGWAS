import subprocess
from helpers.getpaths import get_paths
import config

alias_table_path = get_paths(config.cbio_root)['ensembl']


def get_gene_symbol(gene_id: str) -> str:
    """ Search for the gene id in the alias table and return the gene symbol. The search uses grep so it should also
    work for gene ids with version numbers.
    :param gene_id: gene id to search for
    :return: gene symbol
    # TODO: instead of grepping every time, create a  JSON and do a lookup in it. It will be much faster.
    """
    process = subprocess.Popen(["grep", "-w", gene_id, alias_table_path], stdout=subprocess.PIPE)  # grep the gene id
    sp_output, sp_error = process.communicate()  # Get the output of the process, both stdout and stderr
    if sp_error:  # I don't know if this would even be reached. Probably other errors would be raised first.
        raise ValueError(f"Error when trying to grep the rsid from the sumstats file:\n{sp_error}")
    if not sp_output:  # If the output is empty, the rsid is not in the file.
        raise ValueError(f"Gene ID {gene_id} not found in the ensembl alias file.")
    sp_output_list = sp_output.decode("utf-8").split('\t')  # decode bytes to string and split by tab
    if len(sp_output_list) > 5:
        print(f"WARNING: Multiple entries found with aliases for {gene_id}. Using the first one.")
        sp_output_list = sp_output_list[:5]
        sp_output_list[-1] = sp_output_list[-1].split('\n')[0]  # Remove the newline from the last element
    gene_stable_ID, gene_stable_ID_version, gene_name, gene_synonym, HGNC_symbol = sp_output_list  # Unpack the output

    return gene_name

