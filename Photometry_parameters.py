import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils.utils import calc_total_error
from Stars import getting_stars
from photutils import CircularAperture

def phot_params_from_user(main_paths_obj_red):
    while True:
        
        boolean_inputs = ["Y", "n"]
        '''GETTING INDEX OF THE CALIBRATED IMAGE FOR THE USER TO TEST DESIRED PHOTOMETRY PARAMETERS'''
        
        user_test_img_input = input("Enter index of an image which you want to display for photometry parameter"
                                   " setting and star detection. The index of an \nimage corresponds to its position"
                                  " in the reduced folder. Index of an image: ")
        data_test_img = fits.getdata(main_paths_obj_red[int(user_test_img_input)])     
        print("The path to the image is ",main_paths_obj_red[int(user_test_img_input)])                    
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''GETTING TESTING PHOTOMETRY PARAMATERS FROM THE USER'''
        user_parameter_input = input("Write parameters for the code to display an image with circled detected stars,"
                          " separated by a comma. For additional info to the proper setting of these parameters"
                          " see the photutils documentation to SigmaClip, detect_threshold, \ndetect_sources,"
                          " circular_footprint, make_source_mask, SExtractorBackground, Background2D and"
                          " DAOStarFinder.\nBox_size should be an integer which divides the number of pixels in "
                          " both dimensions. \nThreshold should be an integer which estimates the value above"
                          " which a pixel could be considered as part of a star.\nVmin and vmax define "
                          " the interval of pixels which will be brightened in the plot, allowing for"
                          " objects to be seen \nbetter, similar to SIPS histogram and stretch. \nThe"
                          " parameters to be entered in this order are - sigma (typically 5), nsigma"
                          " (typically 20-30), npixels (typically 5),\nbox_size "
                          " (typically 64), fwhm (typically 3 to 10), threshold, vmin (typically 0),"
                          " vmax (typically 50 to 400), radius of aperture to use in photometry (circle"
                          " of corresponding radius around every detected object will be displayed).\n"
                          "Then press Enter: ")
        
        user_input_array = user_parameter_input.split(",")
        user_input_array = [int(param) for param in user_input_array]
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''
        # mask = pht.make_source_mask(data, nsigma=30, npixels=5, sigclip_sigma=5, sigclip_iters=20, \
        #    dilate_size=11)  deprecated since photutils 1.5.0
        '''GETTING DETECTED STARS AND BACKGROUND OF AN IMAGE'''
        stars, bkg = getting_stars(user_input_array, data_test_img)
        ''''''''''''''''''''''''''''''
        
        '''GETTING POSITIONS OF CENTROIDS OF THE STARS'''
        positions = np.transpose((stars['xcentroid'], stars['ycentroid']))
        apertures = CircularAperture(positions, r=user_input_array[8])
        ''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''PLOTTING THE IMAGE WITH CIRCLED STARS'''
        plt.figure(figsize=(8,8))
        plt.imshow(data_test_img, cmap='Greys_r', origin='lower', vmin=user_input_array[6], \
                   vmax=user_input_array[7], interpolation='nearest')
        apertures.plot(color='red', lw=1.5, alpha=0.5)
        plt.show()
        ''''''''''''''''''''''''''''''''''''''''''''
        
        error = calc_total_error(data_test_img-bkg.background, bkg.background_rms, 0.85)
        print("Median of error is: ", np.median(error))
        
        '''ASKING USER WHETHER TO USE THESE PARAMETERS FOR THE WHOLE PHOTOMETRY'''
        while True:
            yes_or_no = input("Do you wish to use these parameters for the whole photometry? Y/n: ")
            if yes_or_no in boolean_inputs:
                break
            else:
                print("Invalid input")
            
        if yes_or_no == "Y":
            return user_input_array
  
