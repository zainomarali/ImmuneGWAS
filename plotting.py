import seaborn as sns
import matplotlib.pyplot as plt
from resources.tokyo_ge import tokyo_ge_query


def plot_tokyo_ge(gene_id_list):
    """
    Plot the gene expression data from the Tokyo dataset.
    """
    df = tokyo_ge_query(gene_id_list)
    df.set_index("Gene_name", inplace=True)
    values_df = df[df.columns.to_list()[2:]].copy()  # DF without the Gene_id and gene_name columns
    sns.heatmap(values_df)
    plt.show()
