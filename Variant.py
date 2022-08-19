import subprocess
from helpers import getpaths, dbsnp  # Import all the helper functions
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
        self.__cross_reference_dbsnp()  # Sanity check the variant against the dbSNP database.

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
            raise ValueError(f"rsID {rsid} not found in sumstats file.")
        sp_output_list = sp_output.decode("utf-8").split('\t')  # decode bytes to string and split by tab
        chrom = int(sp_output_list[2])
        pos = int(sp_output_list[3])
        OA = sp_output_list[4]  # Note the order. Second allele is the EA in hail output
        EA = sp_output_list[5]

        return cls(rsid, chrom, pos, EA, OA)

    def __cross_reference_dbsnp(self) -> None:
        """
        Cross-reference the variant with the dbSNP database to ensure the rsID-position assignment is correct.
        TODO : check alleles against the dbSNP database too.
        """
        dbsnp_rsid = dbsnp.dbsnp_single_position_query(self.chrom, self.pos)[2]  # dbSNP rsID for the position
        if dbsnp_rsid != self.rsid:
            raise ValueError(f"rsID {self.rsid} does not match dbSNP rsID {dbsnp_rsid}")

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
