import seaborn as sns
import matplotlib.pyplot as plt
from resources.tokyo_ge import tokyo_ge_query


def plot_tokyo_ge(gene_id_list):
    """
    Plot the gene expression data from the Tokyo dataset.
    """
    df = tokyo_ge_query(gene_id_list)
    ordered_cols = ['Plasmablast', 'LDG', 'Neu', 'pDC', 'mDC', 'CL_Mono', 'Int_Mono',
       'CD16p_Mono', 'NC_Mono', 'Naive_B', 'USM_B', 'DN_B', 'SM_B', 'NK',
       'Mem_CD8', 'EM_CD8', 'TEMRA_CD8', 'Naive_CD8', 'Fr_I_nTreg',
       'Naive_CD4', 'Fr_III_T', 'Fr_II_eTreg', 'CM_CD8', 'Th2', 'Tfh', 'Th17',
       'Mem_CD4', 'Th1']
    df.set_index("Gene_name", inplace=True)
    values_df = df[df.columns.to_list()[2:]].copy()  # DF without the Gene_id and gene_name columns
    sns.heatmap(values_df)
    plt.show()
