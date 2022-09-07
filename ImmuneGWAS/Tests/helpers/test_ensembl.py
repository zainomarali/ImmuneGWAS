from ImmuneGWAS.helpers.ensembl import *


def test_get_gene_symbol():
    output = get_gene_symbol("ENSG00000210077")
    assert output == "MT-TV"


def test_get_gene_symbol_not_found():
    output = get_gene_symbol("ENSG00000_made_up_id")
    assert output == "ENSG00000_made_up_id"
