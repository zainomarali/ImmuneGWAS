import subprocess
import pandas as pd
import logging

from ImmuneGWAS.helpers import getpaths, dbsnp, ldlink  # Import all the helper functions
import ImmuneGWAS.resources.immune_GWAS as immune_GWAS
import ImmuneGWAS.config as config


class Results:
    """
    Simple container class to store the results of a variant query.
    This class is always instantiated as a member of a Variant class object.
    """

    eqtl_cat_df = None
    eqtlgen_cis_df = None
    tokyo_eqtl_df = None
    ldtrait_df = None

    def __init__(self):
        pass

    def set_eqtl_cat_df(self, df):
        self.eqtl_cat_df = df

    def set_eqtlgen_cis_df(self, df):
        self.eqtlgen_cis_df = df

    def set_tokyo_eqtl_df(self, df):
        self.tokyo_eqtl_df = df

    def set_ldtrait_df(self, df):
        self.ldtrait_df = df

    def eqtl_cat(self):
        if self.eqtl_cat_df is None:
            logging.info('No eqtl_cat dataframe found.')
            return None
        else:
            return self.eqtl_cat_df

    def eqtlgen_cis(self):
        if self.eqtlgen_cis_df is None:
            logging.info('No eqtlgen_cis dataframe found.')
            return None
        else:
            return self.eqtlgen_cis_df

    def tokyo_eqtl(self):
        if self.tokyo_eqtl_df is None:
            logging.info('No tokyo dataframe found.')
            return None
        else:
            return self.tokyo_eqtl_df

    def ldtrait(self):
        if self.ldtrait_df is None:
            logging.info('No ldtrait dataframe found.')
            return None
        else:
            return self.ldtrait_df


class Variant:
    def __init__(self, rsid: str, chrom: int, pos: int, EA: str, OA: str):
        """
        Initialise a new Variant object.
        :param rsid: The rsID of the variant.
        :param chrom: The chromosome number of the variant.
        :param pos: The position of the variant.
        :param EA: The Effect Allele.
        :param OA: The Other Allele.
        """
        self.rsid = rsid
        self.chrom = chrom
        self.pos = pos
        self.EA = EA
        self.OA = OA
        self.LDblock = None  # Initialize an empty LDblock attribute. LDblock will be a pandas dataframe later.
        self.gwas_phenotypes = []  # Initialize an empty list of gwas_phenotypes.

        self.results = Results()  # Initialize a Results object.

        logging.info(f"Variant {self.rsid} initialised.")
        self.__cross_reference_dbsnp()  # Sanity check the variant against the dbSNP database.
        logging.info(f"Variant {self.rsid} successfully cross-referenced with dbSNP.")
        self.set_LDblock()  # Set the LDblock attribute for the Variant. No input means calculate using LDlink.
        self.__gwas_output_lookup()

    @classmethod
    def from_rsid(cls, rsid: str):
        """
        Alternative constructor for creating a Variant object from a rsid alone.
        NOTE: the current approach is to use grep, but it might be too unreliable

        :param rsid: The rsID of the variant.
        """
        path = getpaths.get_sumstats_path(config.cbio_root)
        process = subprocess.Popen(["grep", "-w", rsid, path], stdout=subprocess.PIPE)  # grep the rsid from the file
        sp_output, sp_error = process.communicate()  # Get the output of the process, both stdout and stderr
        if sp_error:  # I don't know if this would even be reached. Probably other errors would be raised first.
            raise ValueError(f"Error when trying to grep the rsid from the sumstats file:\n{sp_error}")
        if not sp_output:  # If the output is empty, the rsid is not in the file.
            raise ValueError(f"rsID {rsid} not found in the sumstats file.")
        sp_output_list = sp_output.decode("utf-8").split('\t')  # decode bytes to string and split by tab
        chrom = int(sp_output_list[2])
        pos = int(sp_output_list[3])
        OA = sp_output_list[4]  # Note the order. Second allele is the EA in hail output
        EA = sp_output_list[5]

        return cls(rsid, chrom, pos, EA, OA)

    def __cross_reference_dbsnp(self) -> None:
        """
        Cross-reference the variant with the dbSNP database to ensure the rsID-position assignment is correct.
        TODO : check alleles against the dbSNP database too?
        """
        dbsnp_row_list = dbsnp.dbsnp_single_position_query(self.chrom, self.pos)  # Full dbSNP row for the position
        if len(dbsnp_row_list) > 1:  # There's supposed to only be one row, but just in case...
            logging.warning(f"Multiple dbSNP rows found for {self.chrom}:{self.pos}."
                            "Using only the first one for cross-reference.")
        dbsnp_row = dbsnp_row_list[0]
        # The columns in the dbSNP file are: CHROM POS ID REF ALT QUAL FILTER INFO
        dbsnp_rsid = dbsnp_row[2]  # dbSNP rsID for the position
        if dbsnp_rsid != self.rsid:
            raise ValueError(
                f"rsID {self.rsid} does not match dbSNP rsID {dbsnp_rsid} for that position {self.chrom, self.pos}")

    @staticmethod
    def __map_alleles(full_string: str, EA: str, OA: str) -> dict:
        """
        Map the alleles in the LDblock to the EA and OA of the variant.
        The query to LDproxy returns a database with a column called 'Correlated_Alleles'.
        Each entry in column 'Correlated_Alleles' is a string of the form: 'G=C,A=T', meaning that the G allele is
        equivalent to C for the proxy SNP, and that and the A allele is equivalent to T for the proxy SNP.

        :param full_string: The string in the 'Correlated_Alleles' column of the LDproxy output. (E.g. 'G=C,A=T')
        :param EA: The Effect Allele in the Variant object.
        :param OA: The Other Allele in the Variant object.
        :return: A dictionary with the EA and OA as keys, and the corresponding alleles of the proxy as values.
        """
        allele_dict = {}
        a1 = full_string.split(",")[0].split("=")  # Get the first allele
        a2 = full_string.split(",")[1].split("=")  # Get the second allele
        if a1[0] == EA:
            allele_dict[EA] = a1[1]
            allele_dict[OA] = a2[1]
        elif a1[0] == OA:
            allele_dict[OA] = a1[1]
            allele_dict[EA] = a2[1]

        return allele_dict

    def __gwas_output_lookup(self) -> None:
        """
        Look up the GWAS output for the LD block of the variant. What gwas_phenotypes is the LDblock associated with?
        This method gets called by the constructor.
        """
        logging.info(f"Performing lookup of GWAS gwas_phenotypes for variant {self.rsid}'s LDblock.")
        gwas_hits_df = immune_GWAS.immuneGWAS_LDblock_query(self)
        if not gwas_hits_df.empty:
            self.gwas_phenotypes = gwas_hits_df["Phenotype"].unique().tolist()
        else:
            logging.warning(f"No GWAS hits found for the LDblock of {self.rsid}.")

    def get_rsid(self):
        return self.rsid

    def get_fullpos(self):
        return self.chrom, self.pos

    def get_pos(self):
        return self.pos

    def get_chrom(self):
        return self.chrom

    def get_EA(self):
        return self.EA

    def get_OA(self):
        return self.OA

    def get_LDblock(self):
        return self.LDblock

    def get_gwas_phenotypes(self):
        return self.gwas_phenotypes

    def set_LDblock(self, calculate: bool = True, new_df: pd.DataFrame = None) -> None:
        """
        Set the LDblock attribute for the Variant. If calculate is True, the LDblock is calculated by using LDproxy with
        the rsid attribute of the Variant object used for the query. If calculate is False, the LDblock is set to the
        one given as the new_df argument.

        :param calculate: If true, calculate the LDblock from the rsid of the Variant object.
        :param new_df: DataFrame to set as the new LDblock attribute. Only used if calculate is False.
        """
        if calculate:
            df = ldlink.ldproxy(self.rsid)
            if 'Correlated_Alleles' not in df.columns:
                raise ValueError(f"No LDblock found for {self.rsid}. Please choose a different rsID.")

            df['EA'] = df.Correlated_Alleles.apply(lambda x: self.__map_alleles(x, self.EA, self.OA)[self.EA])
            df['OA'] = df.Correlated_Alleles.apply(lambda x: self.__map_alleles(x, self.EA, self.OA)[self.OA])
            df['chrom'] = df.Coord.apply(lambda x: int(x.split(":")[0][-1]))
            df['hg38_pos'] = df.Coord.apply(lambda x: int(x.split(":")[1]))
            df = df[['RS_Number', 'chrom', 'hg38_pos', 'EA', 'OA', 'R2', 'MAF']]
            self.LDblock = df
        else:
            self.LDblock = new_df
