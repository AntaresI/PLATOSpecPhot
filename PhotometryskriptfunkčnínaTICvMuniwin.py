import numpy as np
import photutils as pht
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
import matplotlib.pyplot as plt
from photutils import CircularAperture
from photutils.utils import calc_total_error
from astropy.table import Table
from astropy import table
import glob
import os
import configparser
from pathlib import Path
from astropy.table import Table
from astropy import table
from astropy.time import Time
config = configparser.ConfigParser()
config.read('config.ini')

'''For testing one image, not part of overall script'''
# data = fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Praktikum\\TIC12036_cast1\\tic12036_R_0074.fits')

# mask = pht.make_source_mask(data, nsigma=3, npixels=5, dilate_size=11) # To-DO: What to do with the sigma?
# mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
# print((mean, median, std)) 
# f,axs=plt.subplots(1,2,figsize=(16,8))
# axs[0].imshow(data,vmin=300,vmax=6000,origin='lower')
# axs[0].set_title("data")
# axs[1].imshow(mask,origin='lower')
# axs[1].set_title("mask")


# sigma_clip = SigmaClip(sigma=3.)
# bkg_estimator = pht.SExtractorBackground()
# bkg = pht.Background2D(data, (64, 64), mask=mask,filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
# print(bkg.background_median,bkg.background_rms_median)

# daofind = pht.IRAFStarFinder(fwhm=3, threshold=4.*bkg.background_rms_median,exclude_border=True, sharplo=
# 0.2, sharphi=1.0, roundlo=-1, roundhi=1)    #To-DO: How to adaptively set the threshold? (this tuned to muniwin)
# sources = daofind(data - bkg.background_median)


# positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
# apertures = CircularAperture(positions, r=20.)   # #To-DO: How to properly set r to display big enough circles? 
# plt.figure(figsize=(8,8))
# plt.imshow(data, cmap='Greys_r', origin='lower', vmin=900,vmax=1000, interpolation='nearest')
# apertures.plot(color='red', lw=1.5, alpha=0.5)


# error = calc_total_error(data-bkg.background, bkg.background_rms, 0.85) #0.85 je gain
# print(np.median(error))


#To-Do: automatize this script, meaning : parameters must be input by the user, not into the code itself
'''1.krok'''
reduced = config['PATHS']['reduced']
reduced_path = Path(reduced)
calibrated_data_path = Path(reduced +'\\reduced')
calibrated_data_path.mkdir(exist_ok=True) 
reduced_obj_path = Path(reduced +'\\reduced\\reduced_obj' )
reduced_obj_path.mkdir(exist_ok=True)
main_files_obj_red = os.listdir(reduced_obj_path)
main_paths_obj_red = [reduced_obj_path.__str__()+'\\'+main_files_obj_red[i] for i in range(len(main_files_obj_red))]
radii=[3,4,5,6,7.09,8,10,12,15,20,25,30,35,40,45,50] ## aperture radii in pixels

for i,ifile in enumerate(main_paths_obj_red):
    print("aperture photometry :", i+1,len(main_paths_obj_red),ifile)
    rootname,_ = os.path.splitext(ifile)
    catfile = rootname+'-cat.fits'
    data = fits.getdata(ifile)
    
    ## or first mask sources then estimate the sky background
    mask = pht.make_source_mask(data, nsigma=3, npixels=5, dilate_size=11)
    sigma_clip = SigmaClip(sigma=4.)
    bkg_estimator = pht.SExtractorBackground()
    bkg = pht.Background2D(data, (64, 64), mask=mask,filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    print(bkg.background_median,bkg.background_rms_median)
    
    daofind = pht.IRAFStarFinder(fwhm=3.0, threshold=4.*bkg.background_rms_median,exclude_border=True, sharplo=
    0.2, sharphi=1.0, roundlo=-1, roundhi=1)
    sources = daofind(data - bkg.background)
    positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
    apertures = [pht.CircularAperture(positions, r=r) for r in radii]
    error = calc_total_error(data-bkg.background, bkg.background_rms, 1.07)
    aper_phot = pht.aperture_photometry(data - bkg.background, apertures, error=error)
    print(len(aper_phot))
    print(aper_phot.colnames)
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
        aper_phot.write(catfile,overwrite=True)
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

'''2.krok''' #něco jako matching? 
for i,ifile in enumerate(main_paths_obj_red):
 rootname,_=os.path.splitext(ifile)
 catfile = rootname+'-cat.fits'
 print("calculate shifts :", i+1,len(main_paths_obj_red),ifile)
 
 if i == 0:  #To-Do: user should be able to choose the matching star, for now it's always the first one
     cat1 = Table.read(catfile)
     x1 = cat1['xcenter']
     y1 = cat1['ycenter']
     if 'x_sht' not in cat1.colnames:
         xcol = Table.Column(x1,name='x_sht')
         ycol = Table.Column(y1,name='y_sht')
         cat1.add_columns([xcol,ycol])
     else:
         cat1['x_sht'] = x1
         cat1['y_sht'] = y1
     cat1.write(catfile,overwrite=True)
     
 else:
     cat2 = Table.read(catfile)
     ncat2 = len(cat2)
     
     x2 = cat2['xcenter']
     y2 = cat2['ycenter']
     
     XX = []; YY = []
     
     for j in range(ncat2):
         XX.extend((x1-x2[j]))
         YY.extend((y1-y2[j]))
         
     XX = np.array(XX)
     YY = np.array(YY)
     
     xhist,xbins = np.histogram(XX,range=[-200,200],bins=401)
     yhist,ybins = np.histogram(YY,range=[-200,200],bins=401)
     idx = np.argmax(xhist)
     xsht0 = (xbins[idx]+xbins[idx+1])/2.0
     idx = np.argmax(yhist)
     ysht0 = (ybins[idx]+ybins[idx+1])/2.0
     print("initial shift:",xsht0,ysht0)
     # print(np.abs(XX-xsht0),np.abs(YY-ysht0))
     mask = (np.abs(XX-xsht0)<300) & (np.abs(YY-ysht0)<300) #To-DO: how to set this properly??
     print(mask.sum())
     xsht1 = np.median(XX[mask])
     ysht1 = np.median(YY[mask])
     print("finetuned shift:",xsht1,ysht1)
     
     if 'x_sht' not in cat1.colnames:
         xcol=Table.Column(x2+xsht1,name='x_sht')
         ycol=Table.Column(y2+ysht1,name='y_sht')
         cat2.add_columns([xcol,ycol])
     else:
         cat2['x_sht']=x2+xsht1
         cat2['y_sht']=y2+ysht1
     cat2.write(catfile,overwrite=True)
'''2.krok'''     

data=fits.getdata('D:\\Johnny\\Documents\\Matfyz\\astro\\Diplomová práce\\data-Fryda - kopie2\\data-Fryda\\reduced\\reduced_obj\\z20220905000438.fit')
'''před 3.krokem'''
x_targ,y_targ=(587.2,367.3)  # To-Do: this must be input by user
#x_comp,y_comp=(159.54-1,336.61-1)
#x_vali,y_vali=(111.89-1,358.47-1)
x_comp,y_comp=(562.2,296.3)
x_vali,y_vali=(362.3,312.8)
'''před3.krokem'''
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
print("calculating light curves...")

for i,ifile in enumerate(main_paths_obj_red):
    rootname,_=os.path.splitext(ifile)
    head=fits.getheader(ifile)
    datestr=head['DATE-OBS']
    # timestr=head['UT']
    # datetime=datestr+'T'+timestr.strip()
    t=Time(datestr,format='isot',scale='utc')
    jd=t.mjd
    lc_targ[0,i]=jd
    lc_comp[0,i]=jd
    lc_vali[0,i]=jd
 
    print("MJD: ",datestr,jd)
    catfile=rootname+'-cat.fits'
    print("reading:", i+1,len(main_paths_obj_red),ifile)
 
    cat=fits.getdata(catfile)
    x=cat['x_sht']
    y=cat['y_sht']
    # get target star
    d=np.sqrt((x-x_targ)**2+(y-y_targ)**2)
    idx=np.argmin(d)
    icat=cat[idx]
    dt=d[idx]
    # print(x,y,d,idx,icat,dt,icat.)
    if d[idx]<50:   # TO-DO - how to set this adaptively?
        for j in range(naper):
            lc_targ[j+1,i]=icat['mag_'+str(j)]    #aperture_sum_'+str(j)
            lc_targ[naper+j+1,i]=icat['magerr_'+str(j)]  #aperture_sum_err_'+str(j)
    else:
        lc_targ[1:,i]=np.nan
 
 # get comparison star
    d=np.sqrt((x-x_comp)**2+(y-y_comp)**2)
    idx=np.argmin(d)
    icat=cat[idx]
    dc=d[idx]
    if d[idx]<50:
        for j in range(naper):
            lc_comp[j+1,i]=icat['mag_'+str(j)]
            lc_comp[naper+j+1,i]=icat['magerr_'+str(j)]
    else:
        lc_comp[1:,i]=np.nan
 
    # get validation star
    d=np.sqrt((x-x_vali)**2+(y-y_vali)**2)
    idx=np.argmin(d)
    icat=cat[idx]
    dv=d[idx]
    if d[idx]<50:
        for j in range(naper):
 
            lc_vali[j+1,i]=icat['mag_'+str(j)]
            lc_vali[naper+j+1,i]=icat['magerr_'+str(j)]
    else:
        lc_vali[1:,i]=np.nan
 
    print(dt,dc,dv)

'''3.krok'''

'''4.krok upgrade''' #To-Do: play with the plotting options, what to display, what not to etc.
iaper = 4 # for iaper aperture #To-Do: this must be chosen by the user
rlc_targ = lc_targ[iaper+1,:]-lc_comp[iaper+1,:]
rlc_vali = lc_comp[iaper+1,:]-lc_vali[iaper+1,:]

a1 = 1.0/lc_comp[iaper+1,:]; e1 = lc_targ[iaper+naper+1,:]
a2 = lc_targ[iaper+1,:]/lc_comp[iaper+1,:]**2; e2 = lc_comp[iaper+naper+1,:]
rlcerr_targ = np.sqrt(a1**2*e1**2+a2**2*e2**2)
a1 = 1.0/lc_comp[iaper+1,:]; e1 = lc_vali[iaper+naper+1,:]
a2 = lc_vali[iaper+1,:]/lc_comp[iaper+1,:]**2; e2 = lc_comp[iaper+naper+1,:]
rlcerr_vali = np.sqrt(a1**2*e1**2+a2**2*e2**2)
print('photerr for target/comparison:',np.nanmedian(rlcerr_targ))
print('photerr for validation/comparison:',np.nanmedian(rlcerr_vali))

idx = np.argmin(np.abs(lc_targ[0,:]-59661.9)) #To Do: how to set the subtract adaptively
norm_targ = np.nanmedian(rlc_targ[idx:])
norm_vali = np.nanmedian(rlc_vali[idx:])
tmpx = [np.min(lc_targ[0,:]),np.max(lc_targ[0,:])]

plt.figure(figsize=(16,8))
plt.plot(lc_targ[0,:],rlc_vali,'r.')
# plt.plot(lc_targ[0,:],rlc_vali+0.08,'b.')
plt.plot(tmpx,[1.0,1.0],'g-',linewidth=2)
# plt.plot(tmpx,[1.08,1.08],'g-',linewidth=2)
plt.ylim([-0.45,-0.28])
plt.xlabel('MJD',fontsize=20)
plt.ylabel('$\Delta m$')
plt.gca().invert_yaxis()
plt.title("red: exoplanet transit, blue: validation star")
print(sigma_clipped_stats(2.5*np.log10(rlc_vali),sigma=3,maxiters=3))