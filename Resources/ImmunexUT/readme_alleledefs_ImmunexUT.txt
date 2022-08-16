In order to get full_allele_defs for the ImmunexUT hits I had to use 
the full summary statistics files that are in:

cbio3/projects/Zain_2021/ImmuNEXT_Japan_Cell2021/Full_summarystats/E-GEAD-420.processed/

As per ImmunexUT documentation all of the betas reported refer to the ALT allele. 

I simply pulled out the columns for REF and ALT allele in the full sumstats file. 

Then I renamed REF to OA and ALT to EA to make it consistent with our EA/OA nomenclature.

Should do a sanity check on a few to make sure that it is correct!
