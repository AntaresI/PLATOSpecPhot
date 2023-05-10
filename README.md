# PLATOSpecPhot

PLATOSpecPhot is a pipeline which performs differential photometry on sky images coming from the C4-16000 photometric CCD camera. It consists of 10 scripts and one configuration file.
## Instalation
1) Setup your virtual Python environment in such a way, that the `Python` version is >= 3.9, and that you have these libraries installed (the versions listed are those for which PLATOSpecPhot pipeline has been developed and works, that doesn't necessarily mean that the pipeline doesn't work if you have different versions installed)

    `numpy>=1.24 `
   
    `matplotlib>=3.7.0 `
    
    `astropy>=5.2.1 `
    
    `ccdproc>=2.4.0 `
    
    `photutils>=1.5.0 `
    
2) Download the 10 Python scripts (either manually or by cloning the repository) and place them in one single folder. 
3) Now if you run the script `pipeline.py` (meaning you navigate into the folder with all the scripts and then enter `python pipeline.py` (or python3) into command line) and you get asked the first question about performing calibration, you should have everything ready correctly (you can then simply press `ctrl+c` so that the pipeline stops)

## Preparing the FITS images
The pipeline works with FITS images (files with extension `.fit` or `.fits`) coming from CCD cameras. FITS files have headers (see image below) which can be viewed in a program working with scientific images such as SIPS. These headers contain additional information about the images. For the pipeline to perform the light-curve generating properly, the header MUST contain a row named `UT` with the hour of observation, and a row `DATE-OBS` with the date o observation. If your FITS header doesn't contain these two rows and the time is instead written in some other way, you will have to edit the header and split the information manually. 

Momentálně je zde outdated (a dost možná i ne dobře funkční) verze všech skriptů a žádný návod. Všechny skripty 100% funkční pipeline a návod budou dodány v týdnu mezi 8.5. a 14.5. 2023
