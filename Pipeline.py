import numpy as np
import os
from pathlib import Path
import configparser

from Image_calibration import image_calibration
from Photometry_parameters import phot_params_from_user
from Photometry import LightCurve
from Plotting import plotting
from Centroid_shifting_and_LC_making import centroid_shifting_and_lc_making
boolean_inputs = ["Y", "n"]

'''ASKING WHETHER CALIBRATION OF IMAGES HAS ALREADY BEEN COMPLETED'''

while True:
    after_calibration_or_not = input("Has calibration of raw images already been completed? If yes,"
                                   " the program will continue doing photometry with objects \nfound in"
                                   " calib_obj directory of working directory (it is therefore expected"
                                   " that calibration has been completed and \nthe working directory"
                                   " structure is unchanged), if no then calibration will be performed."
                                   " Y/n: ")
    if after_calibration_or_not in boolean_inputs:
        break
    else:
        print("Invalid input")
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 

'''RUNNING CALIBRATION SCRIPT AND CREATING DIRECTORIES'''

if after_calibration_or_not == "n":
        
        main_paths_obj_calib = image_calibration()
        
        config = configparser.ConfigParser()
        config.read('config.ini')     
        
        working_dir = config['PATHS']['workingdirectory']
        working_dir_path = Path(working_dir)
        
        reduced_lc_path = Path(working_dir +'\\calibrated\\Light_Curve')
        reduced_lc_path.mkdir(exist_ok=True)
        
        reduced_phot_path = Path(working_dir +'\\calibrated\\reduced_obj_phot' )
        reduced_phot_path.mkdir(exist_ok=True)

      
else:
        print("LOADING CALIBRATED IMAGES")
        config = configparser.ConfigParser()
        config.read('config.ini')        
        
        working_dir = config['PATHS']['workingdirectory']
        working_dir_path = Path(working_dir)
        
        reduced_lc_path = Path(working_dir +'\\calibrated\\Light_Curve')
        reduced_lc_path.mkdir(exist_ok=True)
        
        calib_obj_path = Path(working_dir +'\\calibrated\\reduced_obj' )
        calib_obj_path.mkdir(exist_ok=True)
        
        main_files_obj_calib = os.listdir(calib_obj_path)
        main_paths_obj_calib = [calib_obj_path.__str__()+'\\'+main_files_obj_calib[i] for i in range(len(main_files_obj_calib))]
        
        calib_phot_path = Path(working_dir +'\\calibrated\\reduced_obj_phot' )
        calib_phot_path.mkdir(exist_ok=True)
''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''ASKING WHETHER PHOTOMETRY HAS ALREADY BEEN COMPLETED'''

while True:
    after_phot_or_not = input("Has photometry already been completed? If yes the program will load"
                              " the appropriate light curve data and will continue \non to matching and plotting,"
                              " if no then photometry will be performed. Y/n: ")
    if after_phot_or_not in boolean_inputs:
        break
    else:
        print("Invalid input")
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''RUNNING PHOTOMETRY SCRIPTS'''    
    
if after_phot_or_not == "n":
           user_input_array = phot_params_from_user(main_paths_obj_calib)
           targ_lc, comp_lc, vali_lc, radii = LightCurve(main_paths_obj_calib, user_input_array, calib_phot_path, reduced_lc_path)                                  

else:
           while True:
               after_match_or_not = input("Has matching and light curve data generation already been"
                                          " completed? If yes the program will load the appropriate"
                                          " light curve data and will continue on to plotting, if no"
                                          " then matching will be performed. Y/n: ")
               if after_match_or_not in boolean_inputs:
                   break
               else:
                   print("Invalid input")
                   
           if after_match_or_not == "n":
                      radii = np.loadtxt(reduced_lc_path.__str__()+'\\radii.txt')                     
                      targ_lc, comp_lc, vali_lc, radii = centroid_shifting_and_lc_making(main_paths_obj_calib,radii)
                      
           else:
                       
                      targ_lc = np.loadtxt(reduced_lc_path.__str__()+'\\target_lc.txt') 
                      comp_lc = np.loadtxt(reduced_lc_path.__str__()+'\\comp_lc.txt')
                      vali_lc = np.loadtxt(reduced_lc_path.__str__()+'\\vali_lc.txt')
                      radii = np.loadtxt(reduced_lc_path.__str__()+'\\radii.txt')
''''''''''''''''''''''''''''''''


'''PLOTTING LIGHT CURVE'''   
       
plotting(targ_lc, comp_lc, vali_lc, radii)
''''''''''''''''''''''''''



