
import photutils as pht
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint


def getting_apertures(user_input_array, data_img):
    
    sigma_clip = SigmaClip(sigma=user_input_array[0], maxiters=20)      #apparently increasing this to 4 helped with one image that looked completely alright but was still getting the box error, how does this work?
    threshold = detect_threshold(data_img, nsigma=user_input_array[1], sigma_clip=sigma_clip) 
    segment_img = detect_sources(data_img, threshold, npixels=user_input_array[2])  
    footprint = circular_footprint(radius=user_input_array[3])   
    mask = segment_img.make_source_mask(footprint=footprint)    #https://photutils.readthedocs.io/en/stable/background.html#masking-sources
    
    mean, median, std = sigma_clipped_stats(data_img,mask=mask, sigma=user_input_array[0], maxiters=20) #this to display for user the possible optimal threshold value
    
    bkg_estimator = pht.SExtractorBackground(sigma_clip=sigma_clip)
    bkg = pht.Background2D(data_img, (user_input_array[4], user_input_array[4]), mask=mask, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator) #https://photutils.readthedocs.io/en/stable/background.html#d-background-and-noise-estimation
    print("Mean of background is: ", mean," Median of background is: ", bkg.background_median, " Median of standard deviation of background is: ", bkg.background_rms_median)
    #ValueError at pht.Background2D usually indicates something fundamentally wrong with the reduced images (some weird shapes or zero pixel circles due to bad reduction)
    
    # iraffind = pht.IRAFStarFinder(fwhm=5, threshold=17.*bkg.background_rms_median,exclude_border=True, sharplo=
    # 0.2, sharphi=1.0, roundlo=-1, roundhi=1)    #To-DO: Threshold according to photutils doc could probably be set to 5*std (this tuned to muniwin), but testing showed that it's wrong
    # sources = iraffind(data - bkg.background_median)
    
    daofind = pht.DAOStarFinder(fwhm=user_input_array[5], threshold=user_input_array[6],exclude_border=True, sharplo=
    0.2, sharphi=1.0, roundlo=-1, roundhi=1)        
    sources = daofind(data_img - bkg.background_median)
    if sources is None:
        print('Using IRAFStarFinder instead')
        iraffind = pht.IRAFStarFinder(fwhm=user_input_array[5], threshold=user_input_array[6],exclude_border=True, sharplo=
        0.2, sharphi=1.0, roundlo=-1, roundhi=1)
        sources = iraffind(data_img - bkg.background)
        
    return sources, bkg