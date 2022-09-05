import logging

cbio_root = "/home/antton/"  # Path to cbio3
output_folder = "/home/antton/Desktop/ImmuneGWAS_post_output/"  # Destination for output files

if cbio_root == "/default/path/to/cbio/":
    raise ValueError("Default path to cbio3 is used. Please change this in config.py")

if output_folder == "/default/path/to/output/folder/":
    raise ValueError("Default path to output folder is used. Please change this in config.py")

# Initialize logging
logging.basicConfig(filename=output_folder + "ImmuneGWAS.log",
                    filemode='w',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%Y %I:%M:%S',
                    level=logging.INFO)
