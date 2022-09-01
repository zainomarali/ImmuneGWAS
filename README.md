# ImmuneGWAS
Scripts for analysis of GWAS results.

These scripts integrate several eQTL databases to look for genetic effects of SNPs identified by GWAS.


## Installation

You can install the package with pip by running this command from inside the package directory:

```bash
pip install -e ./
```

The following databases are available, listed here under the names used in the scripts:
* **Tokyo**: a database of eQTLs from the [ImmuneXUT](https://www.immunexut.org/) project, from the University of Tokyo.
[Cell. 2021; 184(11): 3006-3021](https://doi.org/10.1016/j.cell.2021.03.056)
* **eQTL cat**: a database of eQTLs from the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) project, by EMBL-EBI.
* **eQTLGen**: a database of eQTLs by the [eQTLGen](https://www.eqtlgen.org/) Consortium. Both cis-eQTL and trans-eQTL 
datasets are available.

The package is installed by typing

    pip install -e .

In the root folder, the package can then be imported as 
    
    import ImmuneGWAS