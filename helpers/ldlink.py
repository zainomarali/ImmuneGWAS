import requests
import pandas as pd
import warnings
from requests.structures import CaseInsensitiveDict


"""
Functions to access LDProxy and LDTrait API from LDLink, using Zain's access token: da0eb217dded

There are two functions: 'ldproxy' and 'ldtrait'
In both cases, the function takes an rsid as input, and returns pandas dataframes that correspond to API access.
By default, they use EUR populations.
To search for Japanese populations for ImmunexUT set "pop='JPT'" when calling function.

Downstream processing needed: 

1) Convert allele defs to EA/OA nomenclature.
2) Use dbSNP to make sure all rsids are always converted to their most recent version
"""

warnings.filterwarnings('ignore')


def ldtrait(rsid, pop='"CEU+FIN+GBR+TSI+IBS"'):
    """
    Access LDtrait and get all the available traits for the requested SNP.

    :param rsid: rsid you are interested in
    :param pop: Populations to use in the LDProxy query. Default is European.
    """

    headers = CaseInsensitiveDict()
    headers["Content-Type"] = "application/json"
    url = 'https://ldlink.nci.nih.gov/LDlinkRest/ldtrait?token=da0eb217dded'

    d1 = '{"snps": '
    d2 = ', "pop": ' + pop + ', "r2_d": "r2", "r2_d_threshold": "0.8", "window": "500000", "genome_build": ' \
                             '"grch38_high_coverage"} '
    rsid = '"' + rsid + '"'
    d = d1 + rsid + d2
    response = requests.post(url, data=d, headers=headers)

    if not response.ok:
        raise RuntimeError(f"LDLink HTTP POST response status code: {response.status_code}. Make sure arguments are "
                           "correct.")
    inputlist = []
    for row in response:
        inputlist.append(row)
    inputlist = [row.decode('UTF-8') for row in inputlist]
    inputlist = "".join(inputlist).split('\n')
    inputlist = [i for i in inputlist if i != '']
    cols = ['Query', 'GWAS Trait', 'RS Number', 'Position (GRCh38)', 'Alleles', 'R2', 'D', 'Risk Allele',
            'Effect Size (95% CI)', 'Beta or OR', 'P-value']
    df = pd.DataFrame(columns=cols)
    for i in range(1, len(inputlist)):
        row = inputlist[i].split('\t')
        if len(row) == len(df.columns):
            df.loc[i] = (inputlist[i].split('\t'))
    return df


def ldproxy(rsid, pop='CEU+FIN+GBR+TSI+IBS', threshold=0.8):
    """
    Access LDproxy and get all the SNPs in LD (R2 > 0.8) with the input SPN.

    :param rsid: rsid you are interested in
    :param pop: Populations to use in the LDProxy query. Default is European.
    :param threshold: R2 threshold. 0.8 by default
    """

    token = 'da0eb217dded'  # this is the token I got for API access can be different for other users
    # these populations were selected for closeness to Nordic population - can also be altered
    params = (
        ('var', rsid),
        ('pop', pop),
        ('r2_d', 'r2'),
        ('token', token),
        ('genome_build', 'grch38_high_coverage')
    )
    response = requests.get('https://ldlink.nci.nih.gov/LDlinkRest/ldproxy', params=params, verify=False)
    if not response.ok:
        raise RuntimeError(f"LDLink HTTP POST response status code: {response.status_code}. Make sure arguments are "
                           "correct.")
    # the requests library produces a response object that must be first converted into a single string
    # and then split by tabs and newlines

    inputlist = []
    for row in response:  # Fill up inputlist with the response from the API
        inputlist.append(row)
    inputlist = [row.decode('UTF-8') for row in inputlist]
    x = inputlist[0]
    for i in range(1, len(inputlist)):
        x = x + inputlist[i]
    x_split = x.split('\n')
    for i in range(len(x_split)):
        x_split[i] = x_split[i].split('\t')
    df = pd.DataFrame(x_split)
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df = df.drop(df.index[-1])
    if "R2" in df.columns:
        df.R2 = df.R2.astype(float)
        df = df[df.R2 > threshold]
    else:
        print("ERROR: df lacks 'R2' column. It might be empty.")
    return df
