
import os
import configparser
from astropy.io import fits
import numpy as np

from Files_and_folders_preparation import file_and_folder_preparation
from Master_image_creation import file_calibration

def image_calibration():
    
    '''LOADING CONFIG FILE'''
    print("LOADING CONFIG FILE")
    config = configparser.ConfigParser()
    config.read('config.ini')
    ''''''''''''''''''''''''''    
    
    '''GETTING THE PATH ARRAYS AND FILE ARRAYS FOR CALIBRATION FRAMES AND SCIENCE FRAMES'''
    image_file_array, path_array, \
    main_files_obj_raw, main_paths_obj_raw = file_and_folder_preparation(config)
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    print('CALIBRATING RAW IMAGES')
    '''GETTING MASTER FRAMES FOR CALIBRATION'''
    calibration_files = file_calibration(image_file_array, path_array)
    ''''''''''''''''''''''''''''''''''''''''''''
    
    counter = 0
    threshold = 3600
    '''PERFORMING THE CALIBRATION (WITH C4-16000 CAMERA SPECIALTY)'''
    for file in main_paths_obj_raw:
        
        with fits.open(file) as obj_raw:
            obj_raw_data, helper = obj_raw[0].data, obj_raw[0].data
         
            obj_where_less = obj_raw_data < threshold
            obj_where_more = obj_raw_data > threshold
            
            helper = helper.astype('float32') # otherwise type casting error (combined darks are also float)
            
            np.subtract(helper, calibration_files[0], out=helper, where=obj_where_less)
            np.subtract(helper, calibration_files[1], out=helper, where=obj_where_more)
            
            fneg = helper < 0
            helper[fneg] = 0    #setting the negative values after dark subtracting to zero
                            #unravel_index(a.argmax(), a.shape)
            helper = helper.astype('float32')
            
            np.divide(helper, calibration_files[2], out=helper, where=obj_where_less)
            np.divide(helper, calibration_files[3], out=helper, where=obj_where_more)
            foverflow= helper > 65535
            helper[foverflow] = 65535 #setting possible values that overflowed back into 16-bit range 
            helper = helper.astype('uint16')
            
            obj_raw[0].data = helper
            obj_raw[0].writeto(path_array[5] / main_files_obj_raw[counter],overwrite=True)
            counter += 1
            print(counter)
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''  
    
    '''CREATING PATH TO FOLDER WITH CALIBRATED RAW IMAGES FOR PHOTOMETRY'''
    main_files_obj_red = os.listdir(path_array[5])
    main_paths_obj_red = [path_array[5].__str__()+'//'+main_files_obj_red[i] for i in range(len(main_files_obj_red))]
    
    return main_paths_obj_red
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
