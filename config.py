import logging

cbio_root = "/media/"  # Path to cbio3
output_folder = "/home/zain/Desktop/"  # Destination for output files

# Initialize logging
logging.basicConfig(filename=output_folder + "ImmuneGWAS.log",
                    filemode='w',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d-%Y %I:%M:%S',
                    level=logging.INFO)
