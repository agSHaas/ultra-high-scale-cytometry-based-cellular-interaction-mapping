"""
This module reads image data from compensated image files (.cif files) into h5 files in python.
extract_CIF() needs to be run for each individual CIF file.

Inspired by: Minh Doan, Holger Hennig
Adapted by: Angelo Jovin Yamachui Sitcheu, Luca Deininger
"""

import bioformats
import bioformats.formatreader
import javabridge

import os
import h5py
import numpy as np
from tqdm import tqdm
import javabridge
import bioformats
 
from multiprocessing import Process


def save_image_as_h5(reader, id, file_path):
    """
    Save image and corresponding mask from a given series in the CIF file to an HDF5 file.

    Args:
        reader (bioformats.ImageReader): The ImageReader instance to read image data.
        id (int): Series ID of the image in the CIF file.
        file_path (str): Path to save the .h5 file.

    Returns:
        int: Always returns 0 upon successful saving.
    """
    with h5py.File(f'{file_path}/image_{id//2}.h5', 'w') as hf:
        image = np.transpose(reader.read(series=id, rescale=False), (2, 0, 1))
        mask = np.transpose(reader.read(series=id+1, rescale=False), (2, 0, 1))
        hf.create_dataset('image', data=image)
        hf.create_dataset('mask', data=mask)
    return 0
  
def my_function(path):
    """
    Process a CIF file to extract images and masks, and save them as .h5 files.

    Args:
        path (str): Path to the CIF file to process.

    Prints:
        Metadata about the file and progress during the processing.
    """
    print(path)
    try:
        javabridge.start_vm(class_path=bioformats.JARS, run_headless=True)
        reader = bioformats.ImageReader(path=path)
        metadata = bioformats.get_omexml_metadata(path)
        xml_reader = bioformats.OMEXML(metadata)
        image_count = xml_reader.get_image_count()
        print("Image count: ", image_count // 2)
        path_to_save = path.split(".")[0]
        print("Path to save:", path_to_save)
        os.makedirs(path_to_save, exist_ok=True)
        with tqdm(total=image_count//2, desc="Processing images", unit="image") as pbar:
            for i in range(0, image_count, 2):
                save_image_as_h5(reader, i, path_to_save)
                pbar.update(1)
        print("Dataset saved as .h5 files in :", path_to_save)
    finally:
        javabridge.kill_vm()
 
def extract_CIF(path): 
    """
    Extract images and masks from a CIF file by running `my_function` in a separate process.

    Args:
        path (str): Path to the CIF file to process.
    """
    p = Process(target=my_function, args=(path,))
    p.start()
    p.join()

def main():
    """
    Main function to demonstrate processing a CIF file.

    Adjust the path to the CIF file as needed.
    """
    cif_file = "/mnt/lsdf_iai-aida/Daten_Deininger/projects/ifc_data_24_07_16/LCMV_experiment_Sample_B1.cif"
    extract_CIF(path=cif_file)
    
if __name__ == "__main__":
    main()