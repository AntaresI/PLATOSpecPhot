import numpy as np
import os
from pathlib import Path
import configparser

from ReductionPlatoSpecv0_3 import reduction
from photometry_parametersPlatoSpec import phot_params_from_user
from LightCurvePlatoSpec import LightCurve
from PlottingPlatoSpec import Plotting

boolean_inputs = ["yes", "no"]

'''PREPARING FILES AND OBTAINING PARAMETERS FROM USER FOR PHOTOMETRY'''
while True:
    after_reduction_or_not = input("Has reduction of raw images already been completed? If yes,\
                                   the program will continue doing photometry with objects found in\
                                   reduced_obj directory of working directory (it is therefore expected \
                                   that reduction has been completed and the working directory\
                                   structure is unchanged), if no then reduction will be performed.\
                                   yes or no: ")
    if after_reduction_or_not in boolean_inputs:
        break
    else:
        print("Invalid input")
        
if after_reduction_or_not == "no":
        
        main_paths_obj_red = reduction()
        
        config = configparser.ConfigParser()
        config.read('config.ini')     
        
        working_dir = config['PATHS']['workingdirectory']
        working_dir_path = Path(working_dir)
        
        reduced_lc_path = Path(working_dir +'\\reduced\\Light_Curve')
        reduced_lc_path.mkdir(exist_ok=True)
        
        reduced_phot_path = Path(working_dir +'\\reduced\\reduced_obj_phot' )
        reduced_phot_path.mkdir(exist_ok=True)
        
else:
        print("LOADING REDUCED IMAGES")
        config = configparser.ConfigParser()
        config.read('config.ini')        
        
        working_dir = config['PATHS']['workingdirectory']
        working_dir_path = Path(working_dir)
        
        reduced_lc_path = Path(working_dir +'\\reduced\\Light_Curve')
        reduced_lc_path.mkdir(exist_ok=True)
        
        reduced_obj_path = Path(working_dir +'\\reduced\\reduced_obj' )
        reduced_obj_path.mkdir(exist_ok=True)
        
        main_files_obj_red = os.listdir(reduced_obj_path)
        main_paths_obj_red = [reduced_obj_path.__str__()+'\\'+main_files_obj_red[i] for i in range(len(main_files_obj_red))]
        
        reduced_phot_path = Path(working_dir +'\\reduced\\reduced_obj_phot' )
        reduced_phot_path.mkdir(exist_ok=True)

while True:
    after_phot_or_not = input("Has photometry already been completed? If yes the program will load\
                                    the appropriate light curve data and will continue on to plotting,\
                                    if no then photometry will be performed, yes or no: ")
    if after_phot_or_not in boolean_inputs:
        break
    else:
        print("Invalid input")
        
if after_phot_or_not == "no":
           user_input_array = phot_params_from_user(main_paths_obj_red)
           targ_lc, comp_lc, vali_lc, radii = LightCurve(main_paths_obj_red, user_input_array, reduced_phot_path, reduced_lc_path)                                  

else:
           targ_lc = np.loadtxt(reduced_lc_path.__str__()+'\\target_lc.txt') 
           comp_lc = np.loadtxt(reduced_lc_path.__str__()+'\\comp_lc.txt')
           vali_lc = np.loadtxt(reduced_lc_path.__str__()+'\\vali_lc.txt')
           radii = np.loadtxt(reduced_lc_path.__str__()+'\\radii.txt')
''''''

'''DOING PHOTOMETRY'''

# while True:
#     after_phot_or_not = input("Has photometry already been completed? If yes the program will load\
#                                     the appropriate light curve data and will continue on to plotting,\
#                                     if no then photometry will be performed, yes or no: ")
                                    
#     if after_phot_or_not in boolean_inputs:
#         break
#     else:
#         print("Invalid input")

# if after_phot_or_not == "no":
    
#            targ_lc, comp_lc, vali_lc, radii = LightCurve(main_paths_obj_red, user_input_array, reduced_phot_path, reduced_lc_path)     

# else:
#            targ_lc = np.loadtxt(reduced_lc_path.__str__()+'\\targ_lc.txt') 
#            comp_lc = np.loadtxt(reduced_lc_path.__str__()+'\\comp_lc.txt')
#            vali_lc = np.loadtxt(reduced_lc_path.__str__()+'\\vali_lc.txt')
#            radii = np.loadtxt(reduced_lc_path.__str__()+'\\radii.txt')
''''''

'''PLOT LIGHT CURVE'''          
Plotting(targ_lc, comp_lc, vali_lc, radii)




