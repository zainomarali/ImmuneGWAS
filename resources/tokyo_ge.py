import pandas as pd
from helpers.getpaths import get_paths
import config

tokyo_ge_path = get_paths(config.cbio_root)['ge_tokyo']


def get_ge_tokyo_dataframe():
    """
    Read the table of all the gene expression data from the Tokyo dataset.
    """
    df = pd.read_csv(tokyo_ge_path, sep='\t')
    return df


def tokyo_ge_query(gene_id_list):
    """
    Perform a lookup by Gene ID in the Tokyo gene expression table. Return the matching rows.
    Input should be a list, but it also works if it's a single gene ID.
    """
    full_df = get_ge_tokyo_dataframe()

    if type(gene_id_list) == str:
        gene_id_list = [gene_id_list]
    gene_df = full_df[full_df['Gene_id'].isin(gene_id_list)]
    return gene_df
