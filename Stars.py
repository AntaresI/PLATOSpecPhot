
import photutils as pht
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint


def getting_stars(user_input_array, data_img):
    
    '''MASKING OF STARS'''
    sigma_clip = SigmaClip(sigma=user_input_array[0], maxiters=20)      #apparently increasing this to 4 helped with one image that looked completely alright but was still getting the box error, how does this work?
    threshold = detect_threshold(data_img, nsigma=user_input_array[1], sigma_clip=sigma_clip) 
    segment_img = detect_sources(data_img, threshold, npixels=user_input_array[2])  
    footprint = circular_footprint(radius=11)   #this value worked well for the data I obtained from testing run
    mask = segment_img.make_source_mask(footprint=footprint)    #https://photutils.readthedocs.io/en/stable/background.html#masking-sources
    ''''''''''''''''''''''''
     
    mean, median, std = sigma_clipped_stats(data_img,mask=mask, sigma=user_input_array[0],
                                            maxiters=20) 
    
    #https://photutils.readthedocs.io/en/stable/background.html#d-background-and-noise-estimation
    
    '''BACKGROUND COMPUTATION''' 
    bkg_estimator = pht.SExtractorBackground(sigma_clip=sigma_clip)
    bkg = pht.Background2D(data_img, (user_input_array[3], user_input_array[3]), mask=mask,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    print("Mean of background is: ", mean," Median of background is: ", bkg.background_median,
          " Median of standard deviation of background is: ", bkg.background_rms_median)
    ''''''''''''''''''''''''''
    
    '''STAR DETECTION VIA DAOPHOT ALGORITHM'''
    daofind = pht.DAOStarFinder(fwhm=user_input_array[4], threshold=user_input_array[5],
                                exclude_border=True, sharplo=0.2, sharphi=1.0, roundlo=-1, 
                                roundhi=1)        
    
    sources = daofind(data_img - bkg.background)
    
    if sources is None:
        print('Using IRAFStarFinder instead')
        iraffind = pht.IRAFStarFinder(fwhm=user_input_array[4], threshold=user_input_array[5],
                                      exclude_border=True, sharplo=0.2, sharphi=1.0, roundlo=-1,
                                      roundhi=1)
        sources = iraffind(data_img - bkg.background)
    ''''''''''''''''''''''''''''''''''''''''''''    
    
    return sources, bkg