import numpy as np
import photutils as pht
from astropy.io import fits
from photutils.utils import calc_total_error
import os
from photutils.aperture import CircularAnnulus, ApertureStats

from Stars import getting_stars
from Centroid_shifting_and_LC_making import centroid_shifting_and_lc_making

def LightCurve(main_paths_obj_red, user_input_array, reduced_phot_path, reduced_lc_path):
    
    '''GETTING METHOD OF BACKGROUND COMPUTATION FROM USER'''
    while True:
        one_or_two = [1,2]
        user_background_method_input = input("Write which method you wish to use as background computation."
                                         " 1 corresponds to the iteration sigma method. 2 corresponds to"
                                         " the annulus method. 1 or 2: ")
    
        user_background_method_input = int(user_background_method_input)
        if user_background_method_input in one_or_two:
            break
        else:
            print("Invalid input")
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''   
    '''IF BACKGROUND COMPUTATION METHOD IS THE ANNULUS METHOD THEN GETTING THE ANNULUS RADII FROM USER'''
    if user_background_method_input == 2:
        user_annulus_radii_input = input("Write the circular annulus inner and outer radius separated by"
                                        " a comma which will be used to estimate the background \naround an"
                                         " object. Note that the numbers you type in will be added to the"
                                         " aperture radius (so if you type in 5,7 then \naperture radius + 5"
                                         " will be used as a radius of the inner circular annulus and"
                                         " aperture radius + 7 will be used as a \nradius for outer circular"
                                        " annulus. Then press Enter: ")
                                          
        user_annulus_radii_input = user_annulus_radii_input.split(",")
        annulus_radii = [int(float(param)) for param in user_annulus_radii_input]                                  
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    '''GETTING APERTURE RADII FROM USER'''
    user_aperture_radii_input = input("Write the aperture radii that you wish to use for photometry"
                                      " separated by a comma, then press Enter: ")
                                       
    user_aperture_radii_input = user_aperture_radii_input.split(",")
    radii = [int(float(param)) for param in user_aperture_radii_input]
    np.savetxt(reduced_lc_path.__str__()+'\\radii.txt',radii) 
    ''''''''''''''''''''''''''''''''''''    
        
    '''FINDING ALL STARS IN SEQUENCE AND GATHERING INFO ABOUT THEIR POSITIONS AND MAGNITUDES'''
    if user_background_method_input == 2:
        for i,ifile in enumerate(main_paths_obj_red):
          
            print("aperture photometry :", i+1,"/",len(main_paths_obj_red),ifile)
            
            '''CREATING PHOT METADATA FILE AND LOADING IMAGE'''
            rootname,_ = os.path.splitext(ifile)
            photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
            data_img = fits.getdata(ifile)                                 
            ''''''''''''''''''''''''''''''''''''''''''''''''''
            
            '''GETTING DETECTED STARS IN THE IMAGE'''
            sources, bkg = getting_stars(user_input_array, data_img)
            ''''''''''''''''''''''''''''''''''''''''''
            
            # https://photutils.readthedocs.io/en/stable/aperture.html#aperture-photometry-with
            #-multiple-apertures-at-each-position
          
            '''GETTING POSITIONS OF CENTROIDS OF THE STARS'''
            positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
            apertures = [pht.CircularAperture(positions, r=r) for r in radii] 
            ''''''''''''''''''''''''''''''''''''''''''''''''
            
            '''PERFORMING APERTURE PHOTOMETRY'''
            error = calc_total_error(data_img-bkg.background, bkg.background_rms, 1/0.85)
            aper_phot = pht.aperture_photometry(data_img, apertures, error=error)
            ''''''''''''''''''''''''''''''''''''
            print("Number of objects detected is: ", len(aper_phot))
            print("X coordinates of centroids on the image are: ",aper_phot['xcenter']," Y coordinates of centroids on the image are: ",aper_phot['ycenter'])
            
            '''SUBTRACTING FLUX FROM CIRCULAR ANNULUS AND CONVERTING FLUX TO MAGNITUDE''' 
            for j in range(len(radii)):
                
                annulus_aperture = CircularAnnulus(positions, r_in=radii[j]+annulus_radii[0], r_out=radii[j]+annulus_radii[1])
                aperstats = ApertureStats(data_img, annulus_aperture)
                annulus_bkg = apertures[j].area*aperstats.mean
                
                fcol = 'aperture_sum_'+str(j)
                ecol = 'aperture_sum_err_'+str(j)
                
                flux = aper_phot[fcol]-annulus_bkg     
                fluxerr = np.sqrt(aper_phot[ecol]**2+annulus_bkg) 
                
                mag = -2.5*np.log10(flux)
                magerr = 2.5/(flux*np.log(10))*fluxerr
                
                aper_phot[fcol] = mag
                aper_phot[ecol] = magerr
                
                aper_phot.rename_column(fcol,'mag_'+str(j))
                aper_phot.rename_column(ecol,'magerr_'+str(j))
                aper_phot.write(photfile,overwrite=True)
            ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''   
            
    else:
        for i,ifile in enumerate(main_paths_obj_red):
          
            print("aperture photometry :", i+1,"/",len(main_paths_obj_red),ifile)
            
            '''CREATING PHOT METADATA FILE AND LOADING IMAGE'''
            rootname,_ = os.path.splitext(ifile)
            photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
            data_img = fits.getdata(ifile)                             
            ''''''''''''''''''''''''''''''''''''''''''''''''''
            
            '''GETTING DETECTED STARS IN THE IMAGE'''
            sources, bkg = getting_stars(user_input_array, data_img)
            ''''''''''''''''''''''''''''''''''''''''''
            
            # https://photutils.readthedocs.io/en/stable/aperture.html#aperture-photometry-with-multiple-apertures-at-each-position
            
            '''GETTING POSITIONS OF CENTROIDS OF THE STARS'''
            positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
            apertures = [pht.CircularAperture(positions, r=r) for r in radii]  
            ''''''''''''''''''''''''''''''''''''''''''''''''
            
            '''PERFORMING APERTURE PHOTOMETRY WITH SUBTRACTED BACKGROUND'''
            error = calc_total_error(data_img-bkg.background, bkg.background_rms, 1/0.85)
            aper_phot = pht.aperture_photometry(data_img-bkg.background, apertures, error=error)
            print("Number of objects detected is: ", len(aper_phot))
            print("X coordinates of centroids on the image are: ",aper_phot['xcenter']," Y coordinates of centroids on the image are: ",aper_phot['ycenter'])
            ''''''''''''''''''''''''''''''''''''
            
            '''CONVERTING FLUX TO MAGNITUDE''' 
            for j in range(len(radii)):
                
                fcol = 'aperture_sum_'+str(j)
                ecol = 'aperture_sum_err_'+str(j)
                
                flux = aper_phot[fcol] 
                fluxerr = aper_phot[ecol] 
                
                mag = -2.5*np.log10(flux)
                magerr = 2.5/(flux*np.log(10))*fluxerr
                
                aper_phot[fcol] = mag
                aper_phot[ecol] = magerr
                
                aper_phot.rename_column(fcol,'mag_'+str(j))
                aper_phot.rename_column(ecol,'magerr_'+str(j))
                aper_phot.write(photfile,overwrite=True)
            ''''''''''''''''''''''''''''''''''''   
            
            
    target_lc, comp_lc, vali_lc, radii = centroid_shifting_and_lc_making(main_paths_obj_red, radii)
    
    return target_lc, comp_lc, vali_lc, radii
    
