import subprocess
from helpers import getpaths, dbsnp, ldlink  # Import all the helper functions
import config


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
        self.LDblock = None  # Initialize an empty LDblock attribute. LDblock is a pandas dataframe

        self.__cross_reference_dbsnp()  # Sanity check the variant against the dbSNP database.
        self.set_LDblock()  # Set the LDblock attribute for the Variant.

    @classmethod
    def from_rsid(cls, rsid: str):
        """
        Alternative constructor for creating a Variant object from a rsid alone.
        NOTE: the current approach is to use grep, but it might be too unreliable. Maybe it's better to use pandas?
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
            print(f"Multiple dbSNP rows found for {self.chrom}:{self.pos}. Using the first one for cross-reference.")
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

    def set_LDblock(self) -> None:
        """
        Set the LDblock attribute for the Variant.
        """
        df = ldlink.ldproxy(self.rsid)
        if 'Correlated_Alleles' not in df.columns:
            raise ValueError(f"No LDblock found for {self.rsid}. Please choose a different rsID.")

        df['EA'] = df.Correlated_Alleles.apply(lambda x: self.__map_alleles(x, self.EA, self.OA)[self.EA])
        df['OA'] = df.Correlated_Alleles.apply(lambda x: self.__map_alleles(x, self.EA, self.OA)[self.OA])
        df['chrom'] = df.Coord.apply(lambda x: int(x.split(":")[0][-1]))
        df['hg38_pos'] = df.Coord.apply(lambda x: int(x.split(":")[1]))
        df = df[['RS_Number', 'chrom', 'hg38_pos', 'EA', 'OA', 'R2', 'MAF']]
        self.LDblock = df
