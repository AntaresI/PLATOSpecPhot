
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils.utils import calc_total_error
from aperturesPlatoSpec import getting_apertures
from photutils import CircularAperture

def phot_params_from_user(main_paths_obj_red):
    while True:
        
        boolean_inputs = ["yes", "no"]
        user_test_img_input = input("Enter index of an image which you want to display for photometry parameter\
                                   setting and star detection. The index of an image corresponds to its position\
                                   in the reduced folder. Index of an image: ")
        data_test_img = fits.getdata(main_paths_obj_red[int(user_test_img_input)])                         
        
        user_parameter_input = input("Write parameters for the code to display an image with circled detected stars,\
                           separated by a comma. For additional info to the proper setting of these parameters \
                           see the photutils documentation to SigmaClip, detect_threshold, detect_sources,\
                           circular_footprint, make_source_mask, SExtractorBackground, Background2D and\
                           DAOStarFinder. Box_size should be an integer which divides the number of pixels in \
                           both dimensions. Threshold should be an integer which estimates the value above\
                           which a pixel could be considered as part of a star. Vmin and vmax define \
                           the interval of pixels which will be brightened in the plot, allowing for \
                           objects to be seen better, similar to SIPS histogram and stretch. The \
                           parameters to be entered in this order are - sigma (typically 5), nsigma \
                           (typically 20-30), npixels (typically 5),\ radius (typically 11), box_size \
                           (typically 64), fwhm (typically 3 to 10), threshold, vmin (typically 0), \
                           vmax (typically 50 to 400), radius of aperture to use in photometry (circle\
                           of corresponding radius around every detected object will be displayed).\
                           Then press Enter: ")
        
        user_input_array = user_parameter_input.split(",")
        user_input_array = [int(param) for param in user_input_array]
        # mask = pht.make_source_mask(data, nsigma=30, npixels=5, sigclip_sigma=5, sigclip_iters=20, dilate_size=11)  deprecated since photutils 1.5.0
        sources, bkg = getting_apertures(user_input_array, data_test_img)
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=user_input_array[9])
        # fig, axes = plt.subplots()
        # axes.imshow(data_test_img, cmap='Greys_r', origin='lower', vmin=user_input_array[7], vmax=user_input_array[8], interpolation='nearest')
        # plt.show()
        plt.figure(figsize=(8,8))
        plt.imshow(data_test_img, cmap='Greys_r', origin='lower', vmin=user_input_array[7], vmax=user_input_array[8], interpolation='nearest')
        apertures.plot(color='red', lw=1.5, alpha=0.5)
        plt.show()
        error = calc_total_error(data_test_img-bkg.background, bkg.background_rms, 0.85) #0.85 je gain
        print("Median of error is: ", np.median(error))
        
        while True:
            yes_or_no = input("Do you wish to use these parameters for the whole photometry? Type yes, or type no: ")
            if yes_or_no in boolean_inputs:
                break
            else:
                print("Invalid input")
            
        if yes_or_no == "yes":
            return user_input_array
  