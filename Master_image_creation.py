import ccdproc as ccdp
from astropy import units as u
import numpy as np
import configparser

def file_calibration(image_file_array, path_array):
     
     '''INITIALIZING ARRAY TO STORE THE CALIBRATION FILES INTO'''
     
     calibration_files = []
     ''''''''''''''''''''''''''''''''''''''''''''''''''''''''
     '''LOADING CONFIG FILE FOR MEMORY LIMIT AND CREATING VARIABLE WITH MEMORY LIMIT USER INPUT'''
     config = configparser.ConfigParser()
     config.read('config.ini')
     memory_limit = float(int(config['OTHERS']['memorylimit'])*10**9)
     ''''''''''''''''''''''''''''''''''''''''''
     
     '''MASTER DARK FOR IMAGES CREATION'''
     
     print('CREATING MASTER DARKS')
     
     raw_darks_hi = image_file_array[1].files_filtered(imagetyp='dark', include_path=True)
     raw_darks_lo = image_file_array[2].files_filtered(imagetyp='dark', include_path=True)
     
     exptimerawdarks = image_file_array[1].summary['exptime'][0]
     
     combined_raw_darks_hi = ccdp.combine(raw_darks_hi, method='median',mem_limit=memory_limit,dtype='int16')
     combined_raw_darks_lo = ccdp.combine(raw_darks_lo, method='median',mem_limit=memory_limit,dtype='int16')
     
     combined_raw_darks_hi.meta['combined'] = True
     combined_raw_darks_lo.meta['combined'] = True
     
     hi_dark_file_name = 'master_raw_dark_hi_{:4.2f}s.fit'.format(exptimerawdarks)
     lo_dark_file_name = 'master_raw_dark_lo_{:4.2f}s.fit'.format(exptimerawdarks)
     
     combined_raw_darks_hi.write(path_array[1] / hi_dark_file_name, overwrite=True)
     combined_raw_darks_lo.write(path_array[1] / lo_dark_file_name, overwrite=True)
     
     calibration_files.append(combined_raw_darks_hi)
     calibration_files.append(combined_raw_darks_lo)
     ''''''''''''''''''''''''''''''''''''
     
     '''MASTER DARK FOR FLATS CREATION'''
     
     print('CREATING MASTER DARKS FOR FLATS')
     # to-do: možnost grafu mediánu do pdf/png, soubor, se zobrazenými darky zasebou, aby bylo hned vidět ustřelení
     flat_darks_hi = image_file_array[5].files_filtered(imagetyp='dark', include_path=True)
     flat_darks_lo = image_file_array[6].files_filtered(imagetyp='dark', include_path=True)
     
     exptimeflatdarks = image_file_array[5].summary['exptime'][0]
     
     combined_flat_darks_hi = ccdp.combine(flat_darks_hi, method='median',mem_limit=memory_limit,dtype='int16')
     combined_flat_darks_lo = ccdp.combine(flat_darks_lo, method='median',mem_limit=memory_limit,dtype='int16')
     
     combined_flat_darks_hi.meta['combined'] = True
     combined_flat_darks_lo.meta['combined'] = True
     
     hi_flat_dark_file_name = 'master_flat_dark_hi_{:4.2f}s.fit'.format(exptimeflatdarks)
     lo_flat_dark_file_name = 'master_flat_dark_lo_{:4.2f}s.fit'.format(exptimeflatdarks)
     
     # combined_flat_darks_hi.data = changing_zeroes(combined_flat_darks_hi.data)
     # combined_flat_darks_lo.data = changing_zeroes(combined_flat_darks_lo.data)
     
     combined_flat_darks_hi.write(path_array[1] / hi_flat_dark_file_name, overwrite=True)
     combined_flat_darks_lo.write(path_array[1] / lo_flat_dark_file_name, overwrite=True)
     ''''''''''''''''''''''''''''''''''''
     
     '''MASTER FLAT FOR IMAGES CREATION'''
     
     print('SUBTRACTING MASTER DARKS FROM FLATS')
     
     for ccd, file_name in image_file_array[3].ccds(imagetyp='flat', return_fname=True): 
         ccd.data = ccd.data.astype('int16')
         ccd = ccdp.subtract_dark(ccd, combined_flat_darks_hi, exposure_time='exptime', exposure_unit=u.second)
         ccd.data = ccd.data.astype('int16')
         ccd.write(path_array[3] / ('flat-subtracted-hi' + file_name),overwrite=True)    
         
     for ccd, file_name in image_file_array[4].ccds(imagetyp='flat', return_fname=True): 
         ccd.data = ccd.data.astype('int16')
         ccd = ccdp.subtract_dark(ccd, combined_flat_darks_lo, exposure_time='exptime', exposure_unit=u.second)
         ccd.data = ccd.data.astype('int16')
         ccd.write(path_array[4] / ('flat-subtracted-lo' + file_name),overwrite=True)    
     
     subtracted_flats_hi = ccdp.ImageFileCollection(path_array[3])
     subtracted_flats_lo = ccdp.ImageFileCollection(path_array[4])
     
     subtractedflats_hi = subtracted_flats_hi.files_filtered(imagetyp='flat', include_path=True)
     subtractedflats_lo = subtracted_flats_lo.files_filtered(imagetyp='flat', include_path=True)
     
     exptimeflatshi = subtracted_flats_hi.summary['exptime'][0]
     exptimeflatslo = subtracted_flats_lo.summary['exptime'][0]
     
     print('CREATING SCALED AND UNSCALED MASTER FLATS')
     
     def inv_median(image):
         return 1/np.median(image) # scaling function for flats
     
     combined_flats_scaled_hi = ccdp.combine(subtractedflats_hi, method='median',scale=inv_median, mem_limit=2000e6,dtype='float32')
     combined_flats_scaled_lo = ccdp.combine(subtractedflats_lo, method='median',scale=inv_median, mem_limit=2000e6,dtype='float32')
     
     combined_flats_hi = ccdp.combine(subtractedflats_hi, method='median',mem_limit=memory_limit,dtype='int16')
     combined_flats_lo = ccdp.combine(subtractedflats_lo, method='median',mem_limit=memory_limit,dtype='int16')
     
     combined_flats_scaled_hi.meta['combined'] = True
     combined_flats_scaled_lo.meta['combined'] = True
     combined_flats_hi.meta['combined'] = True
     combined_flats_lo.meta['combined'] = True
     
     hi_flat_scaled_file_name = 'master_flat_scaled_main_hi_{:4.2f}s.fit'.format(exptimeflatshi) 
     lo_flat_scaled_file_name = 'master_flat_scaled_main_lo_{:4.2f}s.fit'.format(exptimeflatslo) 
     hi_flat_file_name = 'master_flat_main_hi_{:4.2f}s.fit'.format(exptimeflatshi) 
     lo_flat_file_name = 'master_flat_main_lo_{:4.2f}s.fit'.format(exptimeflatslo) 
     
     combined_flats_scaled_hi.write(path_array[2] / hi_flat_scaled_file_name,overwrite=True)
     combined_flats_scaled_lo.write(path_array[2] / lo_flat_scaled_file_name,overwrite=True)
     
     combined_flats_hi.write(path_array[2] / hi_flat_file_name,overwrite=True)
     combined_flats_lo.write(path_array[2] / lo_flat_file_name,overwrite=True)
     
     calibration_files.append(combined_flats_scaled_hi)
     calibration_files.append(combined_flats_scaled_lo)
     
     return calibration_files
