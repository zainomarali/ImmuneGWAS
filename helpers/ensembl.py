import json
from helpers.getpaths import get_paths
import config

alias_table_path = get_paths(config.cbio_root)['ensembl']


def load_json_file_to_dict(path):
    """ Load JSON file to dict."""
    with open(path, 'r') as in_file:
        return json.load(in_file)


def get_gene_symbol(gene_id: str) -> str:
    """ Search for the gene id in the alias JSON file. Return corresponding gene symbol.
    :param gene_id: gene id to search for
    :return: gene symbol
    """
    alias_dict = load_json_file_to_dict(alias_table_path)
    if gene_id in alias_dict:
        return alias_dict[gene_id]
    else:
        raise ValueError("Gene id not found in alias table.")
