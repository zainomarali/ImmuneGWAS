from helpers.getpaths import get_paths
import config
import os

eqtl_cat_path = get_paths(config.cbio_root)['eqtl_cat']+"/curated_crediblesets"  # Path to eQTL catalogue directory. Contains several files.


def get_eqtl_cat_file_list() -> list:
    return os.listdir(eqtl_cat_path)

# TODO: start lookup with the regular 'gene expression' studies (tagged as '_ge_' in the file names)
