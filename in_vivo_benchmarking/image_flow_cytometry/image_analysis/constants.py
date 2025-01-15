"""
This module defines constants used across the project for file paths, segmentation parameters, and column names.

Authors: Luca Deininger, Angelo Jovin Yamachui Sitcheu
"""

# Directories
H5_DIR = '/mnt/lsdf_iai-aida/Daten_Deininger/projects/ifc_data_24_07_16'
CELL_MASK_DIR = '/mnt/lsdf_iai-aida/Daten_Deininger/projects/ifc_data_24_07_16/cellpose_masks'

# CSV files
CSV_FILE_FEATURES = 'image_features.csv'
CSV_FILE_DATA_OVERVIEW = 'data_overview.csv'
IMG_ID = 'img_id'
MASK_FILE = 'mask_file'

# Cell segmentation
BF_CHANNEL = 0  # Brightfield channel index
CELLPOSE_DIAMETER = 20  # Cell diameter for CellPose