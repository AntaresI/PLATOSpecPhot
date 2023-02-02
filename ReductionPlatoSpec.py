
import os
from pathlib import Path
import configparser
import ccdproc as ccdp
from astropy.io import fits
from astropy import units as u
from astropy.nddata import CCDData
from convenience_functions import show_image
import matplotlib.pyplot as plt
from astropy.stats import mad_std
import numpy as np


config = configparser.ConfigParser()
config.read('config.ini')

'''LOADING FILES'''
# To-do: organize path, input and naming variables system (probably done)
print('LOADING FILES')
main_path_obj_raw = Path(config['PATHS']['rawimages'])
main_images_obj_raw = ccdp.ImageFileCollection(main_path_obj_raw)

main_path_dark_hi = Path(config['PATHS']['imagedarkshi'])
main_darks_raw_hi = ccdp.ImageFileCollection(main_path_dark_hi)

main_path_dark_lo = Path(config['PATHS']['imagedarkslo'])
main_darks_raw_lo = ccdp.ImageFileCollection(main_path_dark_lo)

main_path_flat_hi = Path(config['PATHS']['flatshi'])
main_flats_raw_hi = ccdp.ImageFileCollection(main_path_flat_hi)

main_path_flat_lo = Path(config['PATHS']['flatslo'])
main_flats_raw_lo = ccdp.ImageFileCollection(main_path_flat_lo)

flat_path_dark_hi = Path(config['PATHS']['flatdarkshi'])
flat_darks_raw_hi = ccdp.ImageFileCollection(flat_path_dark_hi)

flat_path_dark_lo = Path(config['PATHS']['flatdarkslo'])
flat_darks_raw_lo = ccdp.ImageFileCollection(flat_path_dark_lo)

reduced = config['PATHS']['reduced']
reduced_path = Path(reduced)
''''''

'''CREATING FOLDERS TO SAVE REDUCED IMAGES'''
print('CREATING FOLDERS')
calibrated_data_path = Path(reduced +'\\reduced')
calibrated_data_path.mkdir(exist_ok=True)   # To-do : try reduction again, if the files will be overwritten (they are?)
subtracted_flats_hi_path = Path(reduced +'\\reduced\\subtracted_flats_hi')
subtracted_flats_hi_path.mkdir(exist_ok=True)
subtracted_flats_lo_path = Path(reduced +'\\reduced\\subtracted_flats_lo')
subtracted_flats_lo_path.mkdir(exist_ok=True)
reduced_obj_path = Path(reduced +'\\reduced\\reduced_obj' )
reduced_obj_path.mkdir(exist_ok=True)
''''''

'''ADDING ADU UNIT TO FITS HEADERS OF IMAGES (NOT OBJECT)''' # To-do: is this going to be necessary for main images even?
                                                             # To-do: check for existence of bunit column beforehand maybe?
adu_files = []
zeromax_files = []
print('ADDING ADU UNIT AND CHANGING ZEROS AND MAX VALUES')

main_files_obj_raw = os.listdir(main_path_obj_raw)
main_paths_obj_raw = [main_path_obj_raw.__str__()+'\\'+main_files_obj_raw[i] for i in range(len(main_files_obj_raw))]
adu_files.extend(main_paths_obj_raw)
zeromax_files.extend(main_paths_obj_raw)

main_files_dark_hi = os.listdir(main_path_dark_hi)
main_paths_dark_hi = [main_path_dark_hi.__str__()+'\\'+main_files_dark_hi[i] for i in range(len(main_files_dark_hi))]
adu_files.extend(main_paths_dark_hi)

main_files_dark_lo = os.listdir(main_path_dark_lo)
main_paths_dark_lo = [main_path_dark_lo.__str__()+'\\'+main_files_dark_lo[i] for i in range(len(main_files_dark_lo))]
adu_files.extend(main_paths_dark_lo)

main_files_flat_hi = os.listdir(main_path_flat_hi)
main_paths_flat_hi = [main_path_flat_hi.__str__()+'\\'+main_files_flat_hi[i] for i in range(len(main_files_flat_hi))]
adu_files.extend(main_paths_flat_hi)
zeromax_files.extend(main_paths_flat_hi)

main_files_flat_lo  = os.listdir(main_path_flat_lo)
main_paths_flat_lo  = [main_path_flat_lo.__str__()+'\\'+main_files_flat_lo[i] for i in range(len(main_files_flat_lo))]
adu_files.extend(main_paths_flat_lo)
zeromax_files.extend(main_paths_flat_lo)

flat_files_dark_hi = os.listdir(flat_path_dark_hi)
flat_paths_dark_hi = [flat_path_dark_hi.__str__()+'\\'+flat_files_dark_hi[i] for i in range(len(flat_files_dark_hi))]
adu_files.extend(flat_paths_dark_hi)

flat_files_dark_lo = os.listdir(flat_path_dark_lo)
flat_paths_dark_lo  = [flat_path_dark_lo.__str__()+'\\'+flat_files_dark_lo[i] for i in range(len(flat_files_dark_lo))]
adu_files.extend(flat_paths_dark_lo)

for file in adu_files:
    with fits.open(file, 'update') as f:
        for hdu in f:
            
            '''BUNIT ADDITION'''
            
            hdu.header['bunit']='adu'     
            
            '''FINDING ZEROS/65535 IN ALL IMAGES AND CHANGING THEM TO AVERAGE OF NON-ZERO/NON-65535 VALUES AROUND IT'''
            if file in zeromax_files:
                data = hdu.data 
                
                zeroinds = np.argwhere((data==0) | (data==65535))
                zeroinds = tuple(map(tuple, zeroinds))       #convert array to tuple
                
                if len(zeroinds) != 0:
                 
                    for i in range(len(zeroinds)):
                        value = 0 
                        counterval = 0
                        for j in range(zeroinds[i][0]-1, zeroinds[i][0]+2):
                            for k in range(zeroinds[i][1]-1, zeroinds[i][1]+2):
                                
                                if j >= data.shape[1]:
                                    j -= 5
                                if k >= data.shape[1]:
                                    k -= 5
                                    
                                if data[j,k] != 0:
                                    value += data[j,k]
                                    counterval += 1
                                    
                        value = int(value/counterval)
                        data[zeroinds[i]] = value 
           
            f.close()  
''''''

# '''GETTING EXPTIME FOR OUTPUT NAMES'''
# # To-do: maybe group it together?

'''MASTER DARK FOR IMAGES CREATION'''

print('CREATING MASTER DARKS')
# to-do: možnost grafu mediánu do pdf/png, soubor, se zobrazenými darky zasebou, aby bylo hned vidět ustřelení
raw_darks_hi = main_darks_raw_hi.files_filtered(imagetyp='dark', include_path=True)
raw_darks_lo = main_darks_raw_lo.files_filtered(imagetyp='dark', include_path=True)

exptimerawdarks = main_darks_raw_hi.summary['exptime'][0]

combined_raw_darks_hi = ccdp.combine(raw_darks_hi, method='median',mem_limit=2000e6,dtype='int16')
combined_raw_darks_lo = ccdp.combine(raw_darks_lo, method='median',mem_limit=2000e6,dtype='int16')

combined_raw_darks_hi.meta['combined'] = True
combined_raw_darks_lo.meta['combined'] = True

hi_dark_file_name = 'master_raw_dark_hi_{:6.2f}.fit'.format(exptimerawdarks)
lo_dark_file_name = 'master_raw_dark_lo_{:6.2f}.fit'.format(exptimerawdarks)

combined_raw_darks_hi.write(calibrated_data_path / hi_dark_file_name, overwrite=True)
combined_raw_darks_lo.write(calibrated_data_path / lo_dark_file_name, overwrite=True)
''''''

'''MASTER DARK FOR FLATS CREATION'''

print('CREATING MASTER DARKS FOR FLATS')

flat_darks_hi = flat_darks_raw_hi.files_filtered(imagetyp='dark', include_path=True)
flat_darks_lo = flat_darks_raw_lo.files_filtered(imagetyp='dark', include_path=True)

exptimeflatdarks = flat_darks_raw_hi.summary['exptime'][0]

combined_flat_darks_hi = ccdp.combine(flat_darks_hi, method='median',mem_limit=2000e6,dtype='int16')
combined_flat_darks_lo = ccdp.combine(flat_darks_lo, method='median',mem_limit=2000e6,dtype='int16')

combined_flat_darks_hi.meta['combined'] = True
combined_flat_darks_lo.meta['combined'] = True

hi_flat_dark_file_name = 'master_flat_dark_hi_{:6.2f}.fit'.format(exptimeflatdarks)
lo_flat_dark_file_name = 'master_flat_dark_lo_{:6.2f}.fit'.format(exptimeflatdarks)

combined_flat_darks_hi.write(calibrated_data_path / hi_flat_dark_file_name, overwrite=True)
combined_flat_darks_lo.write(calibrated_data_path / lo_flat_dark_file_name, overwrite=True)
''''''

'''MASTER FLAT FOR IMAGES CREATION'''
#To-do: somehow make more elegant the saving and loading on flats according to hi or lo gain

print('SUBTRACTING MASTER DARKS FROM FLATS')

for ccd, file_name in main_flats_raw_hi.ccds(imagetyp='flat', return_fname=True): 
    ccd.data = ccd.data.astype('int16')
    ccd = ccdp.subtract_dark(ccd, combined_flat_darks_hi, exposure_time='exptime', exposure_unit=u.second)
    ccd.data = ccd.data.astype('int16')
    ccd.write(subtracted_flats_hi_path / ('flat-subtracted-hi' + file_name),overwrite=True)    
    
for ccd, file_name in main_flats_raw_lo.ccds(imagetyp='flat', return_fname=True): 
    ccd.data = ccd.data.astype('int16')
    ccd = ccdp.subtract_dark(ccd, combined_flat_darks_lo, exposure_time='exptime', exposure_unit=u.second)
    ccd.data = ccd.data.astype('int16')
    ccd.write(subtracted_flats_lo_path / ('flat-subtracted-lo' + file_name),overwrite=True)    

subtracted_flats_hi = ccdp.ImageFileCollection(subtracted_flats_hi_path)
subtracted_flats_lo = ccdp.ImageFileCollection(subtracted_flats_lo_path)



subtractedflats_hi = subtracted_flats_hi.files_filtered(imagetyp='flat', include_path=True)
subtractedflats_lo = subtracted_flats_lo.files_filtered(imagetyp='flat', include_path=True)

exptimeflatshi = subtracted_flats_hi.summary['exptime'][0]
exptimeflatslo = subtracted_flats_lo.summary['exptime'][0]

print('CREATING MASTER FLATS')

def inv_median(image):
    return 1/np.median(image) # scaling function for flats

combined_flats_hi = ccdp.combine(subtractedflats_hi, method='median',scale=inv_median, mem_limit=2000e6,dtype='float32')
combined_flats_lo = ccdp.combine(subtractedflats_lo, method='median',scale=inv_median, mem_limit=2000e6,dtype='float32')

combined_flats_hi.meta['combined'] = True
combined_flats_lo.meta['combined'] = True


hi_flat_file_name = 'master_flat_main_hi{:6.2f}.fit'.format(exptimeflatshi) 
lo_flat_file_name = 'master_flat_main_lo{:6.2f}.fit'.format(exptimeflatslo) 

combined_flats_hi.write(calibrated_data_path / hi_flat_file_name,overwrite=True)
combined_flats_lo.write(calibrated_data_path / lo_flat_file_name,overwrite=True)
''''''

'''RAW IMAGES CALIBRATION'''

print('CALIBRATING RAW IMAGES')
# to-do: fix the flat problem with the team
counter = 0
threshold = 3600
for file in main_paths_obj_raw:
    
    with fits.open(file) as f:
        f1, f2 = f[0], f[0]
        f4 = f1.data
        f2 = f2.data
        # dead pixel map, vytvořit mapu true/false hodnot, nahradit ho nějakými průměrnými hodnotami kolem (zjistit), už na začátku, podívat do dat a zjistit pattern těch nul
        # vynést je, otevřít snímek, najít nuly
        # linearity check, je detektor lineární? a kam? 
        # u flatů hodnoty, graf, zkusit vynásobení koeficientem ze SIPS, vynést
        # signal to noise calculator, zadat magnitudu, exptime a vypadne signal to noise
        #zkusti přidat offset, nepřidat offset
        f3 = np.zeros((f4.shape[0],f4.shape[1]))
        fwhereless = f4 < threshold
        fwheremore = f4 > threshold
        f3 = f3 + f2
        np.subtract(f3, combined_raw_darks_hi, out=f3, where=fwhereless)
        np.subtract(f3, combined_raw_darks_lo, out=f3, where=fwheremore)
        
        fneg = f3 < 0
        f3[fneg] = 0    #setting the negative values after dark subtracting to zero (what to do? ask team)
                        #unravel_index(a.argmax(), a.shape)
        f3=f3.astype('float16')
        
        np.divide(f3, combined_flats_hi, out=f3, where=fwhereless)
        np.divide(f3, combined_flats_lo, out=f3, where=fwhereless)
        
        f3 = f3.astype('uint16')
        # f1.data[f1.data<threshold] = (f1.data[f1.data<3600]-combined_raw_darks_hi.data[f1.data<3600]) / combined_flats_hi[f1.data<3600]
        # f1.data[f2.data>3600] = (f1.data[f2.data>3600]-combined_raw_darks_lo.data[f2.data>3600]) / combined_flats_lo[f2.data>3600]
        f1.data = f3
        f1.writeto(reduced_obj_path / main_files_obj_raw[counter],overwrite=True)
        counter += 1
''''''    
'''CREATING PATH TO FOLDER WITH REDUCED RAW IMAGES FOR PHOTOMETRY'''
main_files_obj_red = os.listdir(reduced_obj_path)
main_paths_obj_red = [reduced_obj_path.__str__()+'\\'+main_files_obj_red[i] for i in range(len(main_files_obj_red))]
