import numpy as np
from astropy.io import fits

from astropy.table import Table
import os
import configparser
from pathlib import Path
from astropy.time import Time

def centroid_shifting_and_lc_making(main_paths_obj_red,radii):
    
    config = configparser.ConfigParser()
    config.read('config.ini')     
    
    working_dir = config['PATHS']['workingdirectory']
    
    reduced_lc_path = Path(working_dir +'//calibrated//Light_Curve')
    reduced_lc_path.mkdir(exist_ok=True)
    
    reduced_phot_path = Path(working_dir +'//calibrated//reduced_obj_phot' )
    reduced_phot_path.mkdir(exist_ok=True) 
    
    '''GETTING INDEX OF IMAGE USED FOR CENTROID SHIFTING FROM USER'''
    user_matching_index_input = input("Enter index of an image which you want to use for centroid "
                                       "shifting. The index of an image corresponds to its position"
                                       " in the reduced folder. Index of the image:  ")
    user_matching_index_input = int(user_matching_index_input)  
    print("The path to the referenced image is ",main_paths_obj_red[int(user_matching_index_input)])    
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    '''GETTING COORDINATES OF TARGET,COMPARISON AND VALIDATION STAR'''                          
    user_coordinates_targ_input = input("Write down the coordinates of target star (visual) in form"
                                        " of xcoord, ycoord. Then press enter:   ")
    user_coordinates_comp_input = input("Write down the coordinates of comparison star in form"
                                        " of xcoord, ycoord. Then press enter:   ")
    user_coordinates_vali_input = input("Write down the coordinates of validation star in form"
                                        " of xcoord, ycoord. Then press enter:   ")
    user_coordinates_targ_input = user_coordinates_targ_input.split(",")
    user_coordinates_comp_input = user_coordinates_comp_input.split(",")
    user_coordinates_vali_input = user_coordinates_vali_input.split(",")
    x_targ, y_targ = (int(user_coordinates_targ_input[0]), int(user_coordinates_targ_input[1]))
    x_comp, y_comp = (int(user_coordinates_comp_input[0]), int(user_coordinates_comp_input[1]))
    x_vali, y_vali = (int(user_coordinates_vali_input[0]), int(user_coordinates_vali_input[1]))
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    '''GETTING COORDINATES OF THE STARS FROM MATCHING IMAGE METADATA FILE'''
    rootname,_=os.path.splitext(main_paths_obj_red[user_matching_index_input])
    photfile = str(reduced_phot_path)+"//"+rootname.split("//")[-1]+"-phot.fits"
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
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    '''CENTROID SHIFTING AND LIGHT CURVE MAKING FOR THE WHOLE SEQUENCE'''
    for i,ifile in enumerate(main_paths_obj_red):
        
         '''LOADING PHOT METADATA FILE AND LOADING IMAGE'''
         rootname,_=os.path.splitext(ifile)
         photfile = str(reduced_phot_path)+"//"+rootname.split("//")[-1]+"-phot.fits"
         print("calculate shifts :", i+1,"/",len(main_paths_obj_red),ifile)
         ''''''''''''''''''''''''''''''''''''''''''''''''
         
         if i == user_matching_index_input:
             continue
         
         
         else:
             '''GETTING COORDINATES OF THE STARS FROM IMAGE METADATA FILE'''     
             phot2 = Table.read(photfile)
             nphot2 = len(phot2)
             
             x2 = phot2['xcenter']
             y2 = phot2['ycenter']
             ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
             
             '''INITIALIZING ARRAYS OF DIFFERENCES'''
             xx = []
             yy = []
             ''''''''''''''''''''''''''''''''''''''
             
             '''COMPUTING DIFFERENCES AND APPENDING THEM TO THE ARRAYS OF DIFFERENCES'''
             for j in range(nphot2):
                 xx.extend((x1-x2[j]))
                 yy.extend((y1-y2[j]))
                 
             xx = np.array(xx)
             yy = np.array(yy)
             ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
             
             '''COMPUTING HISTOGRAM OF THE ARRAY AND TAKING THE MOST FREQUENT RESULT'''
             xhist,xbins = np.histogram(xx,bins=500)
             yhist,ybins = np.histogram(yy,bins=500)
             
             idx = np.argmax(xhist)
             xsht0 = (xbins[idx]+xbins[idx+1])/2.0
             
             idx = np.argmax(yhist)
             ysht0 = (ybins[idx]+ybins[idx+1])/2.0
             
             print("Initially computed shift in x direction: ",xsht0)
             print("Initially computed shift in y direction: ",ysht0)
             ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
             
             '''ADDITIONAL FINETUNING(NOT REALLY NECESSARY BECAUSE OF THE LARGE NUMBER OF BINS)'''
             mask = (np.abs(xx-xsht0)<np.abs(2*xsht0)) & (np.abs(yy-ysht0)<np.abs(2*ysht0))
             # print(mask.sum())
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
             print("Finetuned shift in x direction: ",xsht1)
             print("Finetuned shift in y direction: ",ysht1)
             ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
             
             '''PERFORMING SHIFT OF CENTROIDS OF STARS'''
             if 'x_sht' not in phot1.colnames:
                 xcol = Table.Column(x2+xsht0,name='x_sht')
                 ycol = Table.Column(y2+ysht0,name='y_sht')
                 phot2.add_columns([xcol,ycol])
             else:
                 phot2['x_sht'] = x2+xsht0
                 phot2['y_sht'] = y2+ysht0
             phot2.write(photfile,overwrite=True)
             ''''''''''''''''''''''''''''''''''''''''''
             
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''    
    
    '''INITIALIZING ARRAYS FOR VARIABLE, COMPARISON AND VALIDATION STAR MAGNITUDES'''
    num_of_apertures = len(radii)
    num_of_files = len(main_paths_obj_red)
    
    target_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    comp_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    vali_lc = np.zeros((1+2*num_of_apertures,num_of_files))
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    '''CREATING LIGHT CURVES OF VARIABLE, COMPARISON AND VALIDATION STAR'''
    print("calculating light curves...")
    
    for i,ifile in enumerate(main_paths_obj_red):
        
        '''LOADING PHOT METADATA FILE'''
        rootname,_=os.path.splitext(ifile) 
        photfile = str(reduced_phot_path)+"//"+rootname.split("//")[-1]+"-phot.fits"
        head = fits.getheader(ifile)                    
        ''''''''''''''''''''''''''''''
        
        '''GETTING THE MJD DATE OF OBSERVATION'''
        datestr = head['DATE-OBS']
        timestr = head['UT']
        datetime = datestr+'T'+timestr.strip()    #this is tuned to the fits headers time format from camera C4-16000
        # t = Time(datestr,format='isot',scale='utc')
        t = Time(datetime,format='isot',scale='utc')
        mjd = t.mjd
        target_lc[0,i] = mjd
        comp_lc[0,i] = mjd
        vali_lc[0,i] = mjd
        ''''''''''''''''''''''''''''''''''''''''''
        
        # print("MJD: ",datestr,mjd)
        print("Reading image metadata:", i+1,"/",len(main_paths_obj_red),ifile)
        
        phot = fits.getdata(photfile)
        x = phot['x_sht']
        y = phot['y_sht']
        
        
        '''TAKING THE CLOSEST CENTROID TO INPUT VARIABLE STAR COORDINATES'''
        d = np.sqrt((x-x_targ)**2+(y-y_targ)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dt = d[idx]
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''APPEND THE MAGNITUDE OF VARIABLE STAR TO THE RESPECTIVE ARRAY'''
        if d[idx]<60:   
            for j in range(num_of_apertures):
                target_lc[j+1,i] = iphot['mag_'+str(j)]  
                target_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]  
        else:
            target_lc[1:,i] = np.nan
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''TAKING THE CLOSEST CENTROID TO INPUT COMPARISON STAR COORDINATES'''
        d = np.sqrt((x-x_comp)**2+(y-y_comp)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dc = d[idx]
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''APPEND THE MAGNITUDE OF COMPARISON STAR TO THE RESPECTIVE ARRAY'''
        if d[idx]<60:
            for j in range(num_of_apertures):
                comp_lc[j+1,i] = iphot['mag_'+str(j)]
                comp_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]
        else:
            comp_lc[1:,i] = np.nan
     
        '''TAKING THE CLOSEST CENTROID TO INPUT VALIDATION STAR COORDINATES'''
        d = np.sqrt((x-x_vali)**2+(y-y_vali)**2)
        idx = np.argmin(d)
        iphot = phot[idx]
        dv = d[idx]
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''APPEND THE MAGNITUDE OF VALIDATION STAR TO THE RESPECTIVE ARRAY'''
        if d[idx]<60:
            for j in range(num_of_apertures):
     
                vali_lc[j+1,i] = iphot['mag_'+str(j)]
                vali_lc[num_of_apertures+j+1,i] = iphot['magerr_'+str(j)]
        else:
            vali_lc[1:,i] = np.nan
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''    
        print("Distance between input coordinates of target star and closest object is: ",dt)
        print("Distance between input coordinates of comparison star and closest object is: ",dc)
        print("Distance between input coordinates of validation star and closest object is: ",dv)
        
    '''SAVING LIGHTCURVES AND RADII TO TEXT FILE IN LIGHTCURVES FOLDER'''    
    np.savetxt(reduced_lc_path.__str__()+'//target_lc.txt',target_lc) 
    np.savetxt(reduced_lc_path.__str__()+'//comp_lc.txt',comp_lc)
    np.savetxt(reduced_lc_path.__str__()+'//vali_lc.txt',vali_lc) 
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    return target_lc, comp_lc, vali_lc,radii
