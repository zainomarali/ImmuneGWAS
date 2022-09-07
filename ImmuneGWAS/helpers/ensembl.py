import json
import logging

from ImmuneGWAS.helpers.getpaths import get_paths
import ImmuneGWAS.config as config

alias_table_path = get_paths(config.cbio_root)['ensembl']


def load_json_file_to_dict(path: str):
    """ Load JSON file to dict.

    :param path: path to JSON file"""
    with open(path, 'r') as in_file:
        return json.load(in_file)


def get_gene_symbol(gene_id: str) -> str:
    """ Search for the gene id in the alias JSON file. Return corresponding gene symbol. If no gene symbol is found, it
    returns the gene id.

    :param gene_id: gene id to search for
    :return: gene symbol, or ensembl gene id if it cannot be found
    """
    alias_dict = load_json_file_to_dict(alias_table_path)
    if gene_id in alias_dict:
        return alias_dict[gene_id]
    else:
        logging.warning(f"Gene ID {gene_id} not found in alias table.")
        return gene_id
