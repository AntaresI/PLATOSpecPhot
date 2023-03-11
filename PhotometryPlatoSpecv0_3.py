import numpy as np
import photutils as pht
from astropy.io import fits
# from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils.utils import calc_total_error
from astropy.table import Table
# from astropy import table
# import glob
import os
import configparser
from pathlib import Path
from astropy.time import Time
# from photutils.segmentation import SegmentationImage
# from photutils.utils import circular_footprint
import ccdproc as ccdp
from astropy import units as u
from ReductionPlatoSpecv0_3 import reduction

main_paths_obj_red, config = reduction()

'''For testing one image, not part of overall script'''

data = fits.getdata('G:\\reduced\\20220903\\AUMic 20220903\\reduced\\reduced_obj\\z20220903001203.fit')
data = fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\DSTuc 20220823\\reduced\\reduced_obj\\z20220823001015.fit')
data = fits.getdata('G:\\20220830\\DS_Tuc\\z20220830000370.fit')
# To-DO: What to do with the sigma and sigclip_sigma? both have so far affected the desired display of stars
# ...in C4-16000 images
# mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
# print((mean, median, std)) 
# f,axs=plt.subplots(1,2,figsize=(16,8))
# axs[0].imshow(data,vmin=300,vmax=600,origin='lower')
# axs[0].set_title("data")
# axs[1].imshow(mask,origin='lower')
# axs[1].set_title("mask")

mask = pht.make_source_mask(data, nsigma=30, npixels=5, sigclip_sigma=5, sigclip_iters=20, dilate_size=11) 
sigma_clip = SigmaClip(sigma=3.) #apparently increasing this to 4 helped with one image that looked completely alright but was still getting the box error, how does this work?
bkg_estimator = pht.SExtractorBackground()
bkg = pht.Background2D(data, (64, 64), mask=mask,filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
print(bkg.background_median,bkg.background_rms_median)
#ValueError at pht.Background2D usually indicates something fundamentally wrong with the reduced images (some weird shapes or zero pixel circles due to bad reduction)

# iraffind = pht.IRAFStarFinder(fwhm=5, threshold=17.*bkg.background_rms_median,exclude_border=True, sharplo=
# 0.2, sharphi=1.0, roundlo=-1, roundhi=1)    #To-DO: How to adaptively set the threshold? (this tuned to muniwin)
# sources = iraffind(data - bkg.background_median)

daofind = pht.DAOStarFinder(fwhm=10, threshold=35,exclude_border=True, sharplo=
0.2, sharphi=1.0, roundlo=-1, roundhi=1)        
sources = daofind(data - bkg.background_median)

positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=30.)   # #To-DO: How to properly set r to display big enough circles? 
plt.figure(figsize=(8,8))
plt.imshow(data, cmap='Greys_r', origin='lower', vmin=900,vmax=1000, interpolation='nearest')
apertures.plot(color='red', lw=1.5, alpha=0.5)


error = calc_total_error(data-bkg.background, bkg.background_rms, 0.85) #0.85 je gain
print(np.median(error))
'''For testing one image, not part of the overall script'''

#To-Do: automatize this script, meaning that parameters must be input by the user, not into the code itself
'''1.krok'''
reduced = config['PATHS']['reduced']
reduced_path = Path(reduced)
# calibrated_data_path = Path(reduced +'\\reduced')
# calibrated_data_path.mkdir(exist_ok=True) 
# reduced_obj_path = Path(reduced +'\\reduced\\reduced_obj' )
# reduced_obj_path.mkdir(exist_ok=True)
reduced_phot_path = Path(reduced +'\\reduced\\reduced_obj_phot' )
reduced_phot_path.mkdir(exist_ok=True)
# main_files_obj_red = os.listdir(reduced_obj_path)
# main_paths_obj_red = [reduced_obj_path.__str__()+'\\'+main_files_obj_red[i] for i in range(len(main_files_obj_red))]
radii=[15,20,25,27,30,35,40,42,45,47,50,55,60] ## aperture radii in pixels
# radii = [10,11,12,13,14,15,16,17,18,20,22] 
for i,ifile in enumerate(main_paths_obj_red):
  
    print("aperture photometry :", i+1,"/",len(main_paths_obj_red),ifile)
    
    rootname,_ = os.path.splitext(ifile)
    photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
    # photfile = rootname+'-phot.fits'
    data = fits.getdata(ifile)
    
    ## or first mask sources then estimate the sky background
    mask = pht.make_source_mask(data, nsigma=30, npixels=5, sigclip_sigma=5, sigclip_iters=20, dilate_size=11) 
    sigma_clip = SigmaClip(sigma=4.)
    bkg_estimator = pht.SExtractorBackground()
    bkg = pht.Background2D(data, (64, 64), mask=mask,filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    print(bkg.background_median,bkg.background_rms_median)
    
    daofind = pht.DAOStarFinder(fwhm=10, threshold=35,exclude_border=True, sharplo=
    0.2, sharphi=1.0, roundlo=-1, roundhi=1)
    sources = daofind(data - bkg.background)
    if sources is None:
        print('Using IRAFStarFinder instead')
        iraffind = pht.IRAFStarFinder(fwhm=10, threshold=35,exclude_border=True, sharplo=
        0.2, sharphi=1.0, roundlo=-1, roundhi=1)
        sources = iraffind(data - bkg.background)
    positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
    apertures = [pht.CircularAperture(positions, r=r) for r in radii]
    error = calc_total_error(data-bkg.background, bkg.background_rms, 0.85)
    aper_phot = pht.aperture_photometry(data - bkg.background, apertures, error=error)
    print(len(aper_phot))
    print(aper_phot['xcenter'],aper_phot['ycenter'])
    #convert flux to magnitude
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
        
'''1.krok'''       
# f,axs = plt.subplots(1,2,figsize=(16,8))

# data1 = fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000438.fit')
# data2 = fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000460.fit')

# axs[0].imshow(data1,vmin=300,vmax=600,origin='lower')
# axs[0].set_title("image 1")
# axs[1].imshow(data2,vmin=300,vmax=600,origin='lower')
# axs[1].set_title("image 2")

# cat1 = Table.read('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000438-cat.fits')
# cat2 = Table.read('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000460-cat.fits')

# x1=cat1['xcenter']
# y1=cat1['ycenter']
# x2=cat2['xcenter']
# y2=cat2['ycenter']

# ncat1=len(cat1)
# ncat2=len(cat2)

# XX=[]
# YY=[]

# for i in range(ncat2):
#     XX.extend((x1-x2[i]))
#     YY.extend((y1-y2[i]))

# XX = np.array(XX)
# YY = np.array(YY)

# xhist,xbins = np.histogram(XX,range=[-200,200],bins=401)
# yhist,ybins = np.histogram(YY,range=[-200,200],bins=401)
# print(np.median(xhist),np.median(yhist))

# f,axs = plt.subplots(1,2,figsize=(16,8))
# axs[0].hist(XX,range=[-200,200],bins=401)
# axs[0].set_title("x shift")
# axs[1].hist(YY,range=[-200,200],bins=401)
# axs[1].set_title("y shift")

'''2.krok''' 
for i,ifile in enumerate(main_paths_obj_red):
 rootname,_=os.path.splitext(ifile)
 photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
 # photfile1 = rootname+'-phot.fits'
 print("calculate shifts :", i+1,"/",len(main_paths_obj_red),ifile)
 
 if i == 0:  #To-Do: user should be able to choose the matching image, for now it's always the first one
     cat1 = Table.read(photfile)
     x1 = cat1['xcenter']
     y1 = cat1['ycenter']
     if 'x_sht' not in cat1.colnames:
         xcol = Table.Column(x1,name='x_sht')
         ycol = Table.Column(y1,name='y_sht')
         cat1.add_columns([xcol,ycol])
     else:
         cat1['x_sht'] = x1
         cat1['y_sht'] = y1
     cat1.write(photfile,overwrite=True)
     
 else:
     cat2 = Table.read(photfile)
     ncat2 = len(cat2)
     
     x2 = cat2['xcenter']
     y2 = cat2['ycenter']
     
     xx = []
     yy = []
     
     for j in range(ncat2):
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
     if mask.sum()==0:
         constant = 2
         while True:
             constant += 1
             mask = (np.abs(xx-xsht0)<np.abs(constant*xsht0)) & (np.abs(yy-ysht0)<np.abs(constant*ysht0))
             if mask.sum()!=0:
                 print('Tuned sum is:',mask.sum())
                 break
     xsht1 = np.nanmedian(xx[mask])
     ysht1 = np.nanmedian(yy[mask])
     print("finetuned shift:",xsht1,ysht1)
     
     if 'x_sht' not in cat1.colnames:
         xcol=Table.Column(x2+xsht1,name='x_sht')
         ycol=Table.Column(y2+ysht1,name='y_sht')
         cat2.add_columns([xcol,ycol])
     else:
         cat2['x_sht']=x2+xsht1
         cat2['y_sht']=y2+ysht1
     cat2.write(photfile,overwrite=True)
'''2.krok'''     

# data=fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000438.fit')
'''před 3.krokem'''
x_targ,y_targ=(794,681)  # To-Do: this must be input by user
x_comp,y_comp=(1517,1424)
x_vali,y_vali=(1154,1631)
'''před 3.krokem'''
# aper_targ = CircularAperture((x_targ,y_targ), r=40.)
# aper_comp = CircularAperture((x_comp,y_comp), r=40.)
# aper_vali = CircularAperture((x_vali,y_vali), r=40.)
# plt.figure(figsize=(8,8))
# plt.imshow(data, cmap='Greys_r', origin='lower', vmin=300,vmax=600, interpolation='nearest')
# aper_targ.plot(color='red', lw=1.5, alpha=0.5)
# aper_comp.plot(color='cyan', lw=1.5, alpha=0.5)
# aper_vali.plot(color='yellow', lw=1.5, alpha=0.5)
# plt.title('red: target, cyan: comparison, yellow: validation')     

'''3.krok'''
naper = len(radii)
nfiles = len(main_paths_obj_red)

lc_targ = np.zeros((1+2*naper,nfiles))
lc_comp = np.zeros((1+2*naper,nfiles))
lc_vali = np.zeros((1+2*naper,nfiles))

# lc_targ = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 targ.txt',lc_targ)
# lc_comp = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 comp.txt',lc_comp)
# lc_vali = np.loadtxt('G:\\reduced\\20220822\\target_DS_Tuc-filter0\\LightCurves\\DSTuc 20220822 vali.txt',lc_vali)
print("calculating light curves...")

for i,ifile in enumerate(main_paths_obj_red):
    
    rootname,_=os.path.splitext(ifile) 
    photfile = str(reduced_phot_path)+"\\"+rootname.split("\\")[-1]+"-phot.fits"
    # photfile=rootname+'-phot.fits'
    head=fits.getheader(ifile)
    
    datestr=head['DATE-OBS']
    timestr=head['UT']
    datetime=datestr+'T'+timestr.strip()
    t=Time(datetime,format='isot',scale='utc')
    jd=t.mjd
    lc_targ[0,i]=jd
    lc_comp[0,i]=jd
    lc_vali[0,i]=jd
 
    print("MJD: ",datestr,jd)
    print("reading:", i+1,"/",len(main_paths_obj_red),ifile)
 
    phot=fits.getdata(photfile)
    x=phot['x_sht']
    y=phot['y_sht']
    
    # get target star
    d=np.sqrt((x-x_targ)**2+(y-y_targ)**2)
    idx=np.argmin(d)
    icat=phot[idx]
    dt=d[idx]
    # print(x,y,d,idx,icat,dt,icat.)
    if d[idx]<60:   # TO-DO - how to set this adaptively?
        for j in range(naper):
            lc_targ[j+1,i]=icat['mag_'+str(j)]    #aperture_sum_'+str(j)
            lc_targ[naper+j+1,i]=icat['magerr_'+str(j)]  #aperture_sum_err_'+str(j)
    else:
        lc_targ[1:,i]=np.nan
 
 # get comparison star
    d=np.sqrt((x-x_comp)**2+(y-y_comp)**2)
    idx=np.argmin(d)
    icat=phot[idx]
    dc=d[idx]
    if d[idx]<60:
        for j in range(naper):
            lc_comp[j+1,i]=icat['mag_'+str(j)]
            lc_comp[naper+j+1,i]=icat['magerr_'+str(j)]
    else:
        lc_comp[1:,i]=np.nan
 
    # get validation star
    d=np.sqrt((x-x_vali)**2+(y-y_vali)**2)
    idx=np.argmin(d)
    icat=phot[idx]
    dv=d[idx]
    if d[idx]<60:
        for j in range(naper):
 
            lc_vali[j+1,i]=icat['mag_'+str(j)]
            lc_vali[naper+j+1,i]=icat['magerr_'+str(j)]
    else:
        lc_vali[1:,i]=np.nan
 
    print(dt,dc,dv)

'''3.krok'''

'''4.krok upgrade''' #To-Do: play with the plotting options, what to display, what not to etc.
print('WRITING OUT ALL STANDARD DEVIATIONS')
for i in range(len(radii)):
    iaper = i
    rlc_targ = lc_targ[iaper+1,:]-lc_comp[iaper+1,:]
    rlc_vali = lc_comp[iaper+1,:]-lc_vali[iaper+1,:]
    print("Aperture number ",i," radius ",radii[i],"px"," STD C-K: ",round(np.nanstd(rlc_vali),4),"  STD V-C ",round(np.nanstd(rlc_targ),4))
for i in range(len(radii)):
    iaper = 10
     # for iaper aperture #To-Do: this must be chosen by the user
    rlc_targ = lc_targ[iaper+1,:]-lc_comp[iaper+1,:]
    rlc_vali = lc_comp[iaper+1,:]-lc_vali[iaper+1,:]
    print(np.nanstd(rlc_vali),np.nanstd(rlc_targ))
    a1 = 1.0/lc_comp[iaper+1,:]; e1 = lc_targ[iaper+naper+1,:]
    a2 = lc_targ[iaper+1,:]/lc_comp[iaper+1,:]**2; e2 = lc_comp[iaper+naper+1,:]
    rlcerr_targ = np.sqrt(a1**2*e1**2+a2**2*e2**2)
    a1 = 1.0/lc_comp[iaper+1,:]; e1 = lc_vali[iaper+naper+1,:]
    a2 = lc_vali[iaper+1,:]/lc_comp[iaper+1,:]**2; e2 = lc_comp[iaper+naper+1,:]
    rlcerr_vali = np.sqrt(a1**2*e1**2+a2**2*e2**2)
    print('photerr for target/comparison:',np.nanmedian(rlcerr_targ))
    print('photerr for validation/comparison:',np.nanmedian(rlcerr_vali))
    
    # idx = np.argmin(np.abs(lc_targ[0,:])) #To Do: how to set the subtract adaptively
    # norm_targ = np.nanmedian(rlc_targ[idx:])
    # norm_vali = np.nanmedian(rlc_vali[idx:])
    # tmpx = [np.min(lc_targ[0,:]),np.max(lc_targ[0,:])]
    
    plt.figure(figsize=(16,8))
    # plt.plot(lc_targ[0,:],rlc_vali,'b.')
    plt.plot(lc_targ[0,:],rlc_vali,'r.')
    # plt.plot(lc_targ[0,:],rlc_vali+0.08,'b.')
    # plt.plot(tmpx,[1.0,1.0],'g-',linewidth=2)
    # plt.plot(tmpx,[1.08,1.08],'g-',linewidth=2)
    # plt.xlim([59818.0040,59818.1836])
    # # plt.ylim([1.05,1.17])
    plt.ylim([-1.52,-1.10])
    plt.xlabel('MJD',fontsize=20)
    plt.ylabel('$\Delta m$')
    plt.gca().invert_yaxis()
    plt.title('AUMic 20220903 C-K '+str(radii[iaper])+" pxl")  #To-Do: do the plotting properly, automatize titles and somehow limits

'''SAVING THE ARRAYS FOR MAGNITUDES'''

np.savetxt('G:\\reduced\\20220903\\AUMic 20220903\\LightCurves\\AUMic 20220903 targ.txt',lc_targ) #To-Do: organize this
np.savetxt('G:\\reduced\\20220903\\AUMic 20220903\\LightCurves\\AUMic 20220903 comp.txt',lc_comp)
np.savetxt('G:\\reduced\\20220903\\AUMic 20220903\\LightCurves\\AUMic 20220903 vali.txt',lc_vali)



