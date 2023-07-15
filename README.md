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
The pipeline works with FITS images (files with extension `.fits`, AstroPy doesn't work with `.fit` files!!!) coming from CCD cameras. FITS images have headers (see image below) which can be viewed in a program working with scientific images such as SIPS. These headers contain additional information about the images. For the pipeline to perform the light-curve generating properly, the header MUST contain a row named `UT` with the hour of observation, and a row `DATE-OBS` with the date o observation. If your FITS header doesn't contain these two rows and the time is instead written in some other way, you will have to edit the header and split the information manually. Also, please ensure that there is a row called `imagetyp` which contains the information about the type of the image, whether it is a dark image or a flat image (lowercased). Otherwise the pipeline also won't work.  
![hi](/img/picture1.png)
Other than that, nothing else has to be done to the FITS images.
## Preparing the workplace and editing the configuration file
### For C4-16000 camera 
Before the pipeline can run, it needs to know 8 paths to folders - 7 of them contain different kinds of images used for differential photometry of FITS images coming from C4-16000 camera. Open the `config.ini` file (you can edit it in a text editor) and you will see where to write the paths to (see image below). 
![hi](/img/picture2.png)
Warning: do not edit anything else other than the paths themselves!! Otherwise that will break the configuration file and therefore the pipeline processing.

Enter the respective paths to the folders containing the appropriate images. Meaning before you use the pipeline, ensure that: 
1) You have moved the images of an object to a folder whose path you will then enter into the `rawimages` row.
2) You have moved the two sets of dark images obtained through hi-gain channel and lo-gain channel into separate folders whose paths you will then enter into the `imagedarkshi` and `imagedarkslo` row.
3) You have moved the two sets of dark images to be used for flat images correction obtained through hi-gain channel and lo-gain channel into separate folders whose paths you will then enter into the `flatdarkshi` and `flatdarkslo` row.
4) You have moved the two sets of flat images obtained through hi-gain channel and lo-gain channel into separate folders whose paths you will then enter into the `imageflatshi` and `imageflatslo` row.

Then enter the path to a folder where everything will be stored - the calibrated images, the metadata files, the light-curve text files etc.
### For images that don't use the special double channel calibration
In that case there is only one set of dark images, one set of flat images and one set of dark images used for flat calibration. So simply put these into one folder respectively and then enter that respective folder into both `hi` and `lo` rows and the pipeline will calibrate the object images accordingly.
## Working with the pipeline
Once all the folders have been prepared and added to the configuration file, save the configuration file, then navigate in the command line to the folder where you have all 10 scripts stored and enter `python pipeline.py`. 

It is assumed that you are running the pipeline for the first time on this setup of images. As you are completing each step of calibration, photometry and light-curve computation, then assuming you do not move any file from the `calibrated` folder created in the folder you have input into the `workingdirectory` row of configuration file, it will be possible to skip the steps of the pipeline if you for example want to jump right into the plotting. You can think of this as checkpoints that will be granted for completing each step of the calibration and photometry process.
### Step 1: calibration (first question of the pipeline)
The pipeline script will ask you whether calibration of the object images has already been completed. As stated above this is a checkpoint question (if the calibration has already been completed, the pipeline will load the respective calibrated images saved in the `calibrated` folder and jump to the next step, which is photometry). 
Simply enter `n` and calibration of images will start. First, the pipeline needs to manually input a bunit row into the header of each file for the purpose of being able to use the `ccdproc` library which requires this. Then it goes through flat images and removes all zeros in order to not divide by it at the end of calibration process. Then the calibration will start, which means the pipeline will create master images from the set of calibration images provided and then perform the typical arithmetic operation used for calibrating object images in astronomy (see image below with an example). The numbers outputted by the pipeline after creation of the respective master frames are stating which object image is being processed.
![hi](/img/picture3.png)
### Step 2: photometry (second to sixth question of the pipeline)
After calibration of the object images, the pipeline will ask whether photometry has already been completed. This is a second checkpoint. If it hasn't been completed already, simply enter `n`. You will then be asked to enter the index of an image which you want to use for testing the photometry parameters. Indices of images should correspond to the ordering in the `calibrated` folder but in case it doesn't the pipeline outputs for every entered index the path to the chosen image (and if you don't like the image, then after inputting the testing photometry parameters properly, you can then go back to choosing and index after closing the plot). Pick the image that you want to use for testing and enter its index. You will then be prompted to enter the photometry parameters. A brief description of what each parameter means is outputted by the pipeline, but in case you would like a more detailed description, then visit https://photutils.readthedocs.io/en/stable/aperture.html#aperture-photometry-photutils-aperture
https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html#photutils.detection.DAOStarFinder

After entering the parameters separated by a comma, you will then get a plot displayed with circled stars in the image. It is recommended to test on multiple images whether the `DAOStarfinder` really finds the stars that you wish to use for differential photometry. If not then fiddle with the parameters. If you are satisfied with the circled stars, then enter `Y` in the next question. See images below for examples of the inputs and outputs of this step.
![hi](/img/picture4.png)
![hi](/img/picture5.png)

When you are satisfied with the photometry parameters, you will then be asked which method of background computation to use. If you choose annulus method, then you will be asked to input the inner and outer aperture radii, separated by a comma.
Afterwards you need to input all the aperture radii for which you want to compute photometry, separated by a comma. These will then be saved in a text file, so that it can be used as the final checkpoint, when you just want to see the light-curves. 
Once you input the radii, photometry process starts, so wait for all the images to be processed (see image below). For every processed image a metadata file with the extension `-phot` is created which contains all the information output from https://photutils.readthedocs.io/en/stable/api/photutils.aperture.aperture_photometry.html#photutils.aperture.aperture_photometry

![hi](/img/picture6.png)
### Step 3: differential photometry (seventh and eight question)
After finishing photometry, the centroid coordinates have to be shifted relatively to the position of centroids on one reference image, so that they exactly match. For this purpose, the pipeline will ask you to select a reference image. The task of the `Centroid_shifting_and_LC_making.py` script is to load all the centroid coordinates and then with a histogram method determine how much are the centroids of an image shifted in the x and y direction relative to the centroids on the reference image and move them accordingly so that their positions are ideally the same or at least very close (1 to 3 pixels) to the centroids of reference image.

Once you have selected the referenced image, the pipeline will right away ask you to input the x and y coordinates of target star, comparison star and validation star on the reference image, separated by a comma. Then after inputting these coordinates, centroid shifting will start. 

The pipeline will print out the amount of pixels by which it shifted all the centroids in the image (see image below. It is therefore highly recommended to open the reference image in a program that lets you visualize the FITS image, SIPS is very recommended for this, and then another image for which you see the appropriate shift and see whether it matches. If it doesn't, then you either need to redo the photometry with better parameters (since some stars might not have always been detected in the images), or enter the script itself and play with the amount of bins on lines 101 and 102.
![hi](/img/picture7.png)
If the matching was performed correctly, then right away the light-curve generating will start. That means that the pipeline opens each phot file, it finds out which of the centroid coordinates is the closest to the ones you input and it will take the magnitude and append it to an array. This is repeated for all images and the resulting arrays are then saved into text files so that plotting can immediately be done. 

Since the method for detecting which centroid is the proper one consists of finding the minimal distance, it should be verified whether the pipeline is really detecting the correct image. This can be estimated by the output of the pipeline (see image below) where the are 3 numbers that state the distance from the closest centroid coordinates and the ones you input for target star, comparison star and validation star respectively. Since comparison and validation stars are usually dim, they may not always be detected so it can sometimes happen that these numbers will reach hundreds or thousands, which is wrong. The best number to look at is the first one corresponding to target star, if it is somewhere around 1 to 20, then the centroid was identified correctly.
![hi](/img/picture8.png)
If you are not satisfied, you can always quite the pipeline and navigate with the usage of question checkpoints by the pipeline immediately to the centroid shifting and light curve making and try it again. 
### Step 4: plotting (final questions)
Once data for light curves are gathered and saved in an array, the pipeline will proceed to plotting. The pipeline will print out the standard deviations of the difference of magnitudes of target and comparison star (V-C), and comparison and validation star (C-K) along with the mean error of the target star magnitudes for every aperture radius (see image below). 
![hi](/img/picture9.png)
Then the pipeline will ask for a number of the aperture radius for which it loads the magnitudes and plots them versus time in Julian date format which is acquired from the FITS header of each image during light curve array generating (see image below)
![hi](/img/picture10.png)
A good approach to picking the best aperture radius for plotting is to look for the aperture which has the lowest mean error and if there are multiple, then select the one with lower standard deviation (V-C). Then the pipeline will ask whether you wish to plot also errors and then it will plot the light curve. Final question is whether you wish to plot the light curve for a different aperture. If not then the pipeline has finished its task and will exit.

