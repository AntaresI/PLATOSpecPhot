import os
from pathlib import Path
import ccdproc as ccdp
from astropy.io import fits

from changing_zeroesPlatoSpec import changing_zeroes

def file_and_folder_preparation(configfile):
    
    '''INITIALIZATION OF ARRAY WITH IMAGEFILECOLLECTIONS AND PATHS TO THEM'''
    image_file_array = []
    path_array = []
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    '''LOADING AND STORING FILES'''
    print('LOADING REDUCTION FILES')
    
    main_path_obj_raw = Path(configfile['PATHS']['rawimages'])       #To-do: how to check for same exposure time and filters?
    main_images_obj_raw = ccdp.ImageFileCollection(main_path_obj_raw)
    image_file_array.append(main_images_obj_raw)
    
    main_path_dark_hi = Path(configfile['PATHS']['imagedarkshi'])
    main_darks_raw_hi = ccdp.ImageFileCollection(main_path_dark_hi)
    image_file_array.append(main_darks_raw_hi)
    
    main_path_dark_lo = Path(configfile['PATHS']['imagedarkslo'])
    main_darks_raw_lo = ccdp.ImageFileCollection(main_path_dark_lo)
    image_file_array.append(main_darks_raw_lo)
    
    main_path_flat_hi = Path(configfile['PATHS']['flatshi'])
    main_flats_raw_hi = ccdp.ImageFileCollection(main_path_flat_hi)
    image_file_array.append(main_flats_raw_hi)
    
    main_path_flat_lo = Path(configfile['PATHS']['flatslo'])
    main_flats_raw_lo = ccdp.ImageFileCollection(main_path_flat_lo)
    image_file_array.append(main_flats_raw_lo)
    
    flat_path_dark_hi = Path(configfile['PATHS']['flatdarkshi'])
    flat_darks_raw_hi = ccdp.ImageFileCollection(flat_path_dark_hi)
    image_file_array.append(flat_darks_raw_hi)
    
    flat_path_dark_lo = Path(configfile['PATHS']['flatdarkslo'])
    flat_darks_raw_lo = ccdp.ImageFileCollection(flat_path_dark_lo)
    image_file_array.append(flat_darks_raw_lo)   
    ''''''''''''''''''''''''''''''
    
    '''GROUPING TOGETHER FILES WHICH HAVE ADU UNIT IN THE HEADER SO THAT THEY WON'T GET CYCLED OVER'''
    
    main_images_obj_raw_filtered = main_images_obj_raw.files_filtered(bunit="adu",include_path=True)
    main_darks_raw_hi_filtered = main_darks_raw_hi.files_filtered(bunit="adu",include_path=True)
    main_darks_raw_lo_filtered = main_darks_raw_lo.files_filtered(bunit="adu",include_path=True)
    main_flats_raw_hi_filtered = main_flats_raw_hi.files_filtered(bunit="adu",include_path=True)
    main_flats_raw_lo_filtered = main_flats_raw_lo.files_filtered(bunit="adu",include_path=True)
    flat_darks_raw_hi_filtered = flat_darks_raw_hi.files_filtered(bunit="adu",include_path=True)
    flat_darks_raw_lo_filtered = flat_darks_raw_lo.files_filtered(bunit="adu",include_path=True)
    
    files_with_adu = main_images_obj_raw_filtered+main_darks_raw_hi_filtered+main_darks_raw_lo_filtered\
                     + main_flats_raw_hi_filtered+main_flats_raw_lo_filtered+flat_darks_raw_hi_filtered\
                     + flat_darks_raw_lo_filtered  
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    
    '''INITIALIZATION OF ARRAYS THAT WILL CONTAIN ALL FILES AND FLAT FILES'''
    all_files = []
    flat_files = []
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    print('ADDING ADU UNIT AND CHANGING ZERO ADU PIXELS FOR FLATS')
    
    '''GROUPING UP ALL THE FILES, GROUPING UP FLATS, AND THEN FILTERING THEM TO GET A FINAL GROUP
        OF IMAGES THAT NEED TO BE CYCLED OVER TO ADD ADU UNIT INTO THE HEADER AND CHANGE ZERO ADU PIXELS'''
    
    main_files_obj_raw = os.listdir(main_path_obj_raw)
    main_paths_obj_raw = [main_path_obj_raw.__str__()+'\\'+main_files_obj_raw[i] for i in range(len(main_files_obj_raw))]
    all_files.extend(main_paths_obj_raw)
  
    main_files_dark_hi = os.listdir(main_path_dark_hi)
    main_paths_dark_hi = [main_path_dark_hi.__str__()+'\\'+main_files_dark_hi[i] for i in range(len(main_files_dark_hi))]
    all_files.extend(main_paths_dark_hi)
    
    main_files_dark_lo = os.listdir(main_path_dark_lo)
    main_paths_dark_lo = [main_path_dark_lo.__str__()+'\\'+main_files_dark_lo[i] for i in range(len(main_files_dark_lo))]
    all_files.extend(main_paths_dark_lo)
    
    main_files_flat_hi = os.listdir(main_path_flat_hi)
    main_paths_flat_hi = [main_path_flat_hi.__str__()+'\\'+main_files_flat_hi[i] for i in range(len(main_files_flat_hi))]
    all_files.extend(main_paths_flat_hi)
    flat_files.extend(main_paths_flat_hi)
    
    main_files_flat_lo  = os.listdir(main_path_flat_lo)
    main_paths_flat_lo  = [main_path_flat_lo.__str__()+'\\'+main_files_flat_lo[i] for i in range(len(main_files_flat_lo))]
    all_files.extend(main_paths_flat_lo)
    flat_files.extend(main_paths_flat_lo)
    
    flat_files_dark_hi = os.listdir(flat_path_dark_hi)
    flat_paths_dark_hi = [flat_path_dark_hi.__str__()+'\\'+flat_files_dark_hi[i]  \
                          for i in range(len(flat_files_dark_hi))]
    all_files.extend(flat_paths_dark_hi)

    flat_files_dark_lo = os.listdir(flat_path_dark_lo)
    flat_paths_dark_lo  = [flat_path_dark_lo.__str__()+'\\'+flat_files_dark_lo[i] \
                           for i in range(len(flat_files_dark_lo))]
        
    all_files.extend(flat_paths_dark_lo)
   
    
    files_to_add_adu = list(set(all_files)-set(files_with_adu))
    files_to_cycle_over = list(set().union(files_to_add_adu, flat_files))
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    counter = 0
    for file in files_to_cycle_over:
        with fits.open(file, 'update') as f:
            for hdu in f:
                
                '''BUNIT ADDITION'''  
                
                hdu.header['bunit']='adu'     
                ''''''''''''''''''''
                
                
                # '''FINDING ZEROS(formerly 65535 values) IN ALL IMAGES AND CHANGING THEM TO AVERAGE OF NON-ZERO/(formerly NON-65535) VALUES AROUND IT'''
                # '''ZERO IS PROBLEM BECAUSE OF FLAT DIVISION, IT IS A DEAD PIXEL
                '''CHANGING ZERO ADU PIXELS FOR FLATS'''
                
                if file in flat_files: 
                    data = hdu.data 
                    data = changing_zeroes(data)
                '''''' ''''''''''''''''''''''''''''''   
                
                counter += 1
                print(counter)
                f.close()  
    
       
    '''CREATING FOLDERS TO SAVE REDUCED IMAGES AND APPENDING THEIR PATHS TO PATH_ARRAY'''
    
    print('CREATING FOLDERS')
    reduced = configfile['PATHS']['workingdirectory']
    # reduced_path = Path(reduced)
    
    calibrated_data_path = Path(reduced +'\\calibrated')
    calibrated_data_path.mkdir(exist_ok=True)
    path_array.append(calibrated_data_path)
        
    reduced_darks_path = Path(reduced +'\\calibrated\\reduced_darks')
    reduced_darks_path.mkdir(exist_ok=True)
    path_array.append(reduced_darks_path)
    
    reduced_flats_path = Path(reduced +'\\calibrated\\reduced_flats')
    reduced_flats_path.mkdir(exist_ok=True)
    path_array.append(reduced_flats_path)
    
    subtracted_flats_hi_path = Path(reduced +'\\calibrated\\subtracted_flats_hi')
    subtracted_flats_hi_path.mkdir(exist_ok=True)
    path_array.append(subtracted_flats_hi_path)
    
    subtracted_flats_lo_path = Path(reduced +'\\calibrated\\subtracted_flats_lo')
    subtracted_flats_lo_path.mkdir(exist_ok=True)
    path_array.append(subtracted_flats_lo_path)

    reduced_obj_path = Path(reduced +'\\calibrated\\reduced_obj' )
    reduced_obj_path.mkdir(exist_ok=True)
    path_array.append(reduced_obj_path)
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''            
            
    return image_file_array, path_array, main_files_obj_raw, main_paths_obj_raw     

