import logging

cbio_root = "/media/ludvig/"  # Path to cbio3 without the cbio3 name
output_folder = "/media/ludvig/Project_Storage/BloodVariome/Immune_flow_GWAS/results/immuneGWAS_output/"  # Destination for output files

# Initialize logging
logging.basicConfig(filename=output_folder + "ImmuneGWAS.log",
                    filemode='w',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%Y %I:%M:%S',
                    level=logging.INFO)
