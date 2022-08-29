from helpers.ensembl import *


def test_get_gene_symbol():
    output = get_gene_symbol("ENSG00000210077")
    assert output == "MT-TV"
