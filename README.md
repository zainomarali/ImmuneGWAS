
# ImmuneGWAS
Python package for analysis of our immune-cell GWAS results.

This package integrates several eQTL databases to look for genetic effects of SNPs identified by GWAS. I will
allow to elucidate the "story" behind many of these hits by showing the cis and trans eQTL signals that are
associated with the SNP's LD block.

## Getting Started

### Installation

You can install the package by running the following command from the package's root directory:

```bash
pip install -e .
```
The package can then be imported as 
    
    import ImmuneGWAS


### Configuration

After you have installed the package, you will need to update the configuration file. Copy the
`ImmuneGWAS/config.dist.py` file and rename it to `ImmuneGWAS/config.py`. Change the default paths in the file to
the correct paths to cbio3 and the desired output folder.

## Usage

Start your analysis by creating a new Variant object for the SNP you want to analyze. You can either provide the 
complete information about the SNP like so: `Variant(rsid, chromosome, position, effect allele, other allele)` or just
the rsid: `Variant.from_rsid(rsid)`. The latter will only work if the variant is present in our GWAS summary statistics.

Once you have created the Variant object, you can use it to query the eQTL databases. For example, to get all the
cis-eQTLs from eQTLGen for the variant, you can run:

    eqtlgen_cis_LDblock_query(Variant)

This will update the contents of the `.results` attribute of the object with the results from the query. You can get
those results by running:

    Variant.results.get_eqtlgen_cis()

## Resources

The following databases are available, listed here under the names used in the scripts:
* **Tokyo**: a database of eQTLs from the [ImmuneXUT](https://www.immunexut.org/) project, from the University of Tokyo.
[Cell. 2021; 184(11): 3006-3021](https://doi.org/10.1016/j.cell.2021.03.056)
* **eQTL cat**: a database of eQTLs from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) project, by EMBL-EBI.
* **eQTLGen**: a database of eQTLs by the [eQTLGen](https://www.eqtlgen.org/) Consortium. Both cis-eQTL and trans-eQTL 
datasets are available.

Additionally, our preliminary GWAS results, and the output of the LDlink tool LDtrait are also available for each
variant. 
