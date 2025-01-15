"""
This module creates an overview dataframe with one row per event.
The dataframe includes information such as file paths, sample types, image numbers, and unique image IDs.

Authors: Luca Deininger, Angelo Jovin Yamachui Sitcheu
"""

import os
import pandas as pd
from constants import H5_DIR, CSV_FILE_DATA_OVERVIEW

def find_h5_files(root_dir):
    """
    Traverse the given root directory to find all .h5 files.

    Args:
        root_dir (str): Path to the root directory to search for .h5 files.

    Returns:
        list: A list of file paths to all found .h5 files.
    """
    h5_files = []
    file_count = 0
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith('.h5'):
                filepath = os.path.join(dirpath, filename)
                h5_files.append(filepath)
                file_count += 1
    return h5_files

def get_sample_type(sample_name):
    """
    Map sample names to their respective types.

    Args:
        sample_name (str): Name of the sample.

    Returns:
        str: Sample type (e.g., 'infected', 'control', etc.) or 'Unknown sample type' if not found.
    """
    sample_map = {
        "A1": "infected",
        "A2": "infected",
        "A3": "infected",
        "B1": "control",
        "B2": "control",
        "B4": "control",
        "C1": "infected mix",
        "C4": "control mix"
    }
    
    return sample_map.get(sample_name, "Unknown sample type") 


def main():
    """
    Main function to generate an overview dataframe from .h5 files.

    - Finds all .h5 files in the specified directory.
    - Extracts metadata such as sample names, image file names, and image numbers.
    - Adds sample type and unique image ID columns to the dataframe.
    - Saves the dataframe to a CSV file for further analysis.
    """
    h5_files = find_h5_files(H5_DIR)
    df = pd.DataFrame(h5_files, columns=['file'])
    
    df['sample'] = df['file'].apply(lambda x: os.path.dirname(x).split('_')[-1])
    df['img_file'] = df['file'].apply(lambda x: os.path.basename(x))
    df['img_nr'] = df['img_file'].apply(lambda x: x.replace('.h5', '').split('_')[-1])
    df['sample_type'] = df['sample'].apply(lambda x: get_sample_type(x))
    df['img_id'] = df['sample'] + '_' + df['img_nr'].astype('str')

    df.to_csv(CSV_FILE_DATA_OVERVIEW, index=False)
    print('Created overview table')

if __name__ == "__main__":
    main()