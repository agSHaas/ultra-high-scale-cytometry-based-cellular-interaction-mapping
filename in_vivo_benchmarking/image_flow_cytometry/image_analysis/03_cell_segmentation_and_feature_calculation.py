"""
This module performs cell segmentation and subsequent cell feature extraction.

Authors: Luca Deininger, Angelo Jovin Yamachui Sitcheu
"""

import numpy as np
import math
from cellpose import models
import os
from os.path import join as opj
import pandas as pd
import h5py
from tqdm import tqdm
from skimage.measure import regionprops_table
from constants import (
    CSV_FILE_FEATURES, 
    CSV_FILE_DATA_OVERVIEW, 
    CELL_MASK_DIR, 
    CELLPOSE_DIAMETER, 
    BF_CHANNEL, 
    IMG_ID, 
    MASK_FILE
)

def cell_segmentation():
    """
    Perform cell segmentation using CellPose and save masks as .npy files.

    - Reads the overview CSV file for image data.
    - Processes each image using the CellPose model.
    - Saves the segmentation mask as a .npy file for each image.
    """
    df = pd.read_csv(CSV_FILE_DATA_OVERVIEW)
    model = models.Cellpose(gpu=True, model_type="cyto2")

    for _, row in tqdm(df.iterrows(), total=len(df), desc="CellPose", unit="image"):
        outfile = opj(CELL_MASK_DIR, f'{row["img_id"]}.npy')
        if os.path.isfile(outfile):
            continue  # Skip if the mask file already exists
        with h5py.File(row['file'], 'r') as hf:
            img = hf.get('image')[BF_CHANNEL, :, :]
            mask, _, _, _ = model.eval(img, diameter=CELLPOSE_DIAMETER, channels=[0, 0])
            np.save(outfile, mask.astype('uint8'))

    print('CellPose finished')

def feature_calculation():
    """
    Calculate cell features from segmentation masks and save them to a CSV file.

    - Reads the overview CSV file for image data.
    - Adds paths to the corresponding mask files.
    - Extracts features such as label, eccentricity, area, axis lengths, and perimeter.
    - Saves the aggregated features as a CSV file.
    """
    df = pd.read_csv(CSV_FILE_DATA_OVERVIEW)
    df[MASK_FILE] = df[IMG_ID].apply(lambda x: opj(CELL_MASK_DIR, f'{x}.npy'))

    features = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Feature Calculation", unit="image"):
        if not os.path.isfile(row[MASK_FILE]):
            continue  # Skip if the mask file does not exist
        mask = np.load(row[MASK_FILE])
        
        props = regionprops_table(
            mask, 
            properties=['label', 'eccentricity', 'area', 'axis_minor_length', 'axis_major_length', 'perimeter']
        )
        x = pd.DataFrame(props)
        if x.shape[0] == 0:
            x[IMG_ID] = [row[IMG_ID]]
        else:
            x[IMG_ID] = [row[IMG_ID]] * x.shape[0]
        features.append(x)
    df_feat = pd.concat(features)
    df_feat['circularity'] = df_feat.apply(lambda row: (4 * math.pi * row['area']) / (row['perimeter']**2), axis=1)
    df_feat.to_csv(CSV_FILE_FEATURES, index=False)
    print('Features written to CSV')

def main():
    """
    Main function to run cell segmentation and feature extraction.

    - Calls the cell_segmentation function to generate masks.
    - Calls the feature_calculation function to extract features.
    """
    os.makedirs(CELL_MASK_DIR, exist_ok=True)
    cell_segmentation()
    feature_calculation()
    
if __name__ == "__main__":
    main()