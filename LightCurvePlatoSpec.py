
import numpy as np
import photutils as pht
from astropy.io import fits
from photutils.utils import calc_total_error
from astropy.table import Table
import os
from astropy.time import Time
from photutils.aperture import CircularAnnulus, ApertureStats

from aperturesPlatoSpec import getting_apertures


def LightCurve(main_paths_obj_red, user_input_array, reduced_phot_path, reduced_lc_path):
    
    user_aperture_radii_input = input("Write the aperture radii that you wish to use for photometry \
                                       separated by a comma, then press Enter: ")
                                       
    user_aperture_radii_input = user_aperture_radii_input.split(",")
    radii = [int(param) for param in user_aperture_radii_input]
    
    user_annulus_radii_input = input("Write the circular annulus inner and outer radius separated by \
                                      a comma which will be used to estimate the background around an \
                                      object. Note that the numbers you type in will be added to the \
                                      aperture radius (so if you type in 5,7 then aperture radius + 5 \
                                      will be used as a radius of the inner circular annulus and \
                                      aperture radius + 7 will be used as a radius for outer circular \
                                      annulus. Then press Enter: ")
                                      
    user_annulus_radii_input = user_annulus_radii_input.split(",")
    annulus_radii = [int(param) for param in user_annulus_radii_input]                                  
    
    '''FINDING ALL STARS IN IMAGE AND GATHERING INFO ABOUT THEIR POSITIONS AND MAGNITUDES'''
    
    for i,ifile in enumerate(main_paths_obj_red):
      
        print("aperture photometry :", i+1,"/",len(main_paths_obj_red),ifile)
        
        rootname,_ = os.path.splitext(ifile)
        photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
        data_img = fits.getdata(ifile)                                   # photfile = rootname+'-phot.fits'
        
        sources, bkg = getting_apertures(user_input_array, data_img)
        positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
        apertures = [pht.CircularAperture(positions, r=r) for r in radii]  # https://photutils.readthedocs.io/en/stable/aperture.html#aperture-photometry-with-multiple-apertures-at-each-position
        error = calc_total_error(data_img-bkg.background, bkg.background_rms, 0.85)
        aper_phot = pht.aperture_photometry(data_img, apertures, error=error)
        print(len(aper_phot))
        print(aper_phot['xcenter'],aper_phot['ycenter'])
        
        '''SUBTRACTING FLUX FROM CIRCULAR ANNULUS AND CONVERTING FLUX TO MAGNITUDE''' 
        for j in range(len(radii)):
            
            annulus_aperture = CircularAnnulus(positions, r_in=radii[j]+annulus_radii[0], r_out=radii[j]+annulus_radii[1])
            aperstats = ApertureStats(data_img, annulus_aperture)
            annulus_bkg = apertures[j].area*aperstats.mean
            
            fcol = 'aperture_sum_'+str(j)
            ecol = 'aperture_sum_err_'+str(j)
            
            flux = aper_phot[fcol]-annulus_bkg
            fluxerr = np.sqrt(aper_phot[ecol]**2+annulus_bkg) #annulus_bkg je square root na druhou
            
            mag = -2.5*np.log10(flux)
            magerr = 2.5/(flux*np.log(10))*fluxerr
            
            aper_phot[fcol] = mag
            aper_phot[ecol] = magerr
            
            aper_phot.rename_column(fcol,'mag_'+str(j))
            aper_phot.rename_column(ecol,'magerr_'+str(j))
            aper_phot.write(photfile,overwrite=True)
    ''''''
          
    '''GETTING USER INPUT FOR MATCHING AND LIGHTCURVE COMPUTATION'''
    user_matching_index_input = input("Enter index of an image which you want to use for matching. The index of an image corresponds to its position\
                                        in the reduced folder. Index of the image: ")
    user_matching_index_input = int(user_matching_index_input)  
                                  
    user_coordinates_targ_input = input("Write down the coordinates of target star in form of xcoord, ycoord. Then press enter: ")
    user_coordinates_comp_input = input("Write down the coordinates of comparison star in form of xcoord, ycoord. Then press enter:   ")
    user_coordinates_vali_input = input("Write down the coordinates of validation star in form of xcoord, ycoord. Then press enter:   ")
    user_coordinates_targ_input = user_coordinates_targ_input.split(",")
    user_coordinates_comp_input = user_coordinates_comp_input.split(",")
    user_coordinates_vali_input = user_coordinates_vali_input.split(",")
    x_targ, y_targ = (int(user_coordinates_targ_input[0]), int(user_coordinates_targ_input[1]))
    x_comp, y_comp = (int(user_coordinates_comp_input[0]), int(user_coordinates_comp_input[1]))
    x_vali, y_vali = (int(user_coordinates_vali_input[0]), int(user_coordinates_vali_input[1]))
    ''''''
    

    '''GETTING COORDINATES OF MATCHING IMAGE AND PREPARING COLUMNS ACCORDINGLY'''
    rootname,_=os.path.splitext(main_paths_obj_red[user_matching_index_input])
    photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
    phot1 = Table.read(photfile)
    x1 = phot1['xcenter']
    y1 = phot1['ycenter']
    if 'x_sht' not in phot1.colnames:
        xcol = Table.Column(x1,name='x_sht')
        ycol = Table.Column(y1,name='y_sht')
        phot1.add_columns([xcol,ycol])
    else:
        phot1['x_sht'] = x1
        phot1['y_sht'] = y1
    phot1.write(photfile,overwrite=True)
    ''''''
    
    
    for i,ifile in enumerate(main_paths_obj_red):
        
         rootname,_=os.path.splitext(ifile)
         photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
         print("calculate shifts :", i+1,"/",len(main_paths_obj_red),ifile)
         
         if i == user_matching_index_input:
             continue
             
         else:
             phot2 = Table.read(photfile)
             nphot2 = len(phot2)
             
             x2 = phot2['xcenter']
             y2 = phot2['ycenter']
             
             xx = []
             yy = []
             
             for j in range(nphot2):
                 xx.extend((x1-x2[j]))
                 yy.extend((y1-y2[j]))
                 
             xx = np.array(xx)
             yy = np.array(yy)
             
             xhist,xbins = np.histogram(xx,bins=500)
             yhist,ybins = np.histogram(yy,bins=500)
             # xhist,xbins = np.histogram(XX,range=[-200,200],bins=401)
             # yhist,ybins = np.histogram(YY,range=[-200,200],bins=401)
             
             idx = np.argmax(xhist)
             xsht0 = (xbins[idx]+xbins[idx+1])/2.0
             idx = np.argmax(yhist)
             ysht0 = (ybins[idx]+ybins[idx+1])/2.0
             print("initial shift:",xsht0,ysht0)
             # print(np.abs(XX-xsht0),np.abs(YY-ysht0))
             
             mask = (np.abs(xx-xsht0)<np.abs(2*xsht0)) & (np.abs(yy-ysht0)<np.abs(2*ysht0)) #To-DO: how to set this properly??
             print(mask.sum())
             if mask.sum() == 0:
                 constant = 2
                 while True:
                     constant += 1
                     mask = (np.abs(xx-xsht0)<np.abs(constant*xsht0)) & (np.abs(yy-ysht0)<np.abs(constant*ysht0))
                     if mask.sum() != 0:
                         print('tuned sum is:',mask.sum())
                         break
             xsht1 = np.nanmedian(xx[mask])
             ysht1 = np.nanmedian(yy[mask])
             print("finetuned shift:",xsht1,ysht1)
             
             if 'x_sht' not in phot1.colnames:
                 xcol = Table.Column(x2+xsht1,name='x_sht')
                 ycol = Table.Column(y2+ysht1,name='y_sht')
                 phot2.add_columns([xcol,ycol])
             else:
                 phot2['x_sht'] = x2+xsht1
                 phot2['y_sht'] = y2+ysht1
             phot2.write(photfile,overwrite=True)
    ''''''    
    
    '''3.krok'''
    num_of_apertures = len(radii)
    num_of_files = len(main_paths_obj_red)
    
    target_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    comp_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    vali_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    
    # target_lc = np.zeros((2+2*num_of_apertures,num_of_files))
    # comp_lc = np.zeros((2+2*num_of_apertures,num_of_files))
    # vali_lc = np.zeros((2+2*num_of_apertures,num_of_files))
    
    # lc_targ = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 targ.txt',lc_targ)
    # lc_comp = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 comp.txt',lc_comp)
    # lc_vali = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 vali.txt',lc_vali)
    print("calculating light curves...")
    
    for i,ifile in enumerate(main_paths_obj_red):
        
        rootname,_=os.path.splitext(ifile) 
        photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
        head = fits.getheader(ifile)                        # photfile=rootname+'-phot.fits'
        
        
        datestr = head['DATE-OBS']
        # timestr = head['UT']
        # datetime = datestr+'T'+timestr.strip()
        t = Time(datestr,format='isot',scale='utc')
        mjd = t.mjd
        target_lc[0,i] = mjd
        comp_lc[0,i] = mjd
        vali_lc[0,i] = mjd
     
        print("MJD: ",datestr,mjd)
        print("reading:", i+1,"/",len(main_paths_obj_red),ifile)
     
        phot = fits.getdata(photfile)
        x = phot['x_sht']
        y = phot['y_sht']
        
        # get target star
        d = np.sqrt((x-x_targ)**2+(y-y_targ)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dt = d[idx]
        # print(x,y,d,idx,icat,dt,icat.)
        if d[idx]<60:   # TO-DO - how to set this adaptively?
            for j in range(num_of_apertures):
                target_lc[j+1,i] = iphot['mag_'+str(j)]    #aperture_sum_'+str(j)
                target_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]  #aperture_sum_err_'+str(j)
        else:
            target_lc[1:,i] = np.nan
     
     # get comparison star
        d = np.sqrt((x-x_comp)**2+(y-y_comp)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dc = d[idx]
        if d[idx]<60:
            for j in range(num_of_apertures):
                comp_lc[j+1,i] = iphot['mag_'+str(j)]
                comp_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]
        else:
            comp_lc[1:,i] = np.nan
     
        # get validation star
        d = np.sqrt((x-x_vali)**2+(y-y_vali)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dv = d[idx]
        if d[idx]<60:
            for j in range(num_of_apertures):
     
                vali_lc[j+1,i] = iphot['mag_'+str(j)]
                vali_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]
        else:
            vali_lc[1:,i] = np.nan
     
        print(dt,dc,dv)
    
    '''SAVING LIGHTCURVES AND RADII TO TEXT FILE IN LIGHTCURVES FOLDER'''    
    np.savetxt(reduced_lc_path.__str__()+'\\target_lc.txt',target_lc) 
    np.savetxt(reduced_lc_path.__str__()+'\\comp_lc.txt',comp_lc)
    np.savetxt(reduced_lc_path.__str__()+'\\vali_lc.txt',vali_lc) 
    np.savetxt(reduced_lc_path.__str__()+'\\radii.txt',radii) 
    ''''''
    
    return target_lc, comp_lc, vali_lc, radii