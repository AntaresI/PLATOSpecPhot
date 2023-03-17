
import numpy as np
import matplotlib.pyplot as plt


def Plotting(targ_lc, comp_lc, vali_lc, radii):
    
    while True:
        boolean_inputs = ["yes", "no"]
        num_of_apertures = len(radii)
        
        print('WRITING OUT ALL STANDARD DEVIATIONS')
        
        for i in range(len(radii)):
            iaper = i
            targ_lc_diff = targ_lc[iaper+1,:]-comp_lc[iaper+1,:]
            vali_lc_diff = comp_lc[iaper+1,:]-vali_lc[iaper+1,:]
            print("Aperture number ",i," radius ",radii[i],"px"," STD C-K: ",round(np.nanstd(vali_lc_diff),4),"  STD V-C ",round(np.nanstd(targ_lc_diff),4))
       
        while True:
            which_aper_num = input("Enter the aperture number which you want to use for plotting: ")
            
            if int(which_aper_num) in range(len(radii)):
                break
            else:
                print("This number doesn't correspond to any aperture.")
                
        aper_num = int(which_aper_num)
             
        targ_lc_diff = targ_lc[aper_num+1,:]-comp_lc[aper_num+1,:]
        vali_lc_diff = comp_lc[aper_num+1,:]-vali_lc[aper_num+1,:]
        
        print(np.nanstd(vali_lc_diff),np.nanstd(targ_lc_diff))
        
        err_targ = targ_lc[aper_num+num_of_apertures+1,:]
        err_comp = comp_lc[aper_num+num_of_apertures+1,:]
        err_vali = vali_lc[aper_num+num_of_apertures+1,:]
        diff_lc_targ_err = np.sqrt(err_targ**2+err_comp**2)
        diff_lc_vali_err = np.sqrt(err_comp**2+err_vali**2)
        
        # print('photerr for target/comparison:',np.nanmedian(diff_lc_targ_err))
        # print('photerr for validation/comparison:',np.nanmedian(diff_lc_vali_err))
        
        # idx = np.argmin(np.abs(lc_targ[0,:]))
        # norm_targ = np.nanmedian(rlc_targ[idx:])
        # norm_vali = np.nanmedian(rlc_vali[idx:])
        # tmpx = [np.min(lc_targ[0,:]),np.max(lc_targ[0,:])]
        while True:
            errorbars = input("Do you wish to plot errorbars? Yes or no:")
            
            if errorbars in boolean_inputs:
                break
            else:
                print("Invalid input.")
        
        # while True:
        #         title_yes_no = input("Do you want to include a title of your choice? Yes or no: ")
                
        #         if title_yes_no in boolean_inputs:
        #             break
        #         else:
        #             print("Invalid input.")
                    
        if errorbars == "yes" :
            
            # if title_yes_no == "yes" :
            #     title = input("Enter title of graph: ")
            plt.figure(1,figsize=(16,8))
            plt.errorbar(vali_lc[0,:],vali_lc_diff,yerr=diff_lc_vali_err,fmt='b.',ecolor='black',capsize=2)
            plt.xlabel('MJD',fontsize=20)
            plt.ylabel('$\Delta m$')
            plt.gca().invert_yaxis()
            plt.title('C-K '+str(int(radii[aper_num]))+" pxl")
                
            # plt.title('AUMic 20220826 V-C '+str(radii[aper_num])+" pxl") 
            
            plt.figure(2,figsize=(16,8))
            plt.errorbar(targ_lc[0,:],targ_lc_diff,yerr=diff_lc_targ_err,fmt='r.',ecolor='black',capsize=2)
            plt.xlabel('MJD',fontsize=20)
            plt.ylabel('$\Delta m$')
            plt.gca().invert_yaxis()
            plt.title('V-C '+str(int(radii[aper_num]))+" pxl")
            # plt.title('AUMic 20220826 V-C '+str(radii[aper_num])+" pxl") 
            # plt.xlim([59818.0040,59818.1836])
            # plt.ylim([1.05,1.17])
            # plt.ylim([-3.22,-3.10])
            plt.show()
            
        else:
            plt.figure(1,figsize=(16,8))
            plt.errorbar(vali_lc[0,:],vali_lc_diff,yerr=None,fmt='b.',ecolor='black',capsize=2)
            plt.xlabel('MJD',fontsize=20)
            plt.ylabel('$\Delta m$')
            plt.gca().invert_yaxis()
            plt.title('C-K '+str(int(radii[aper_num]))+" pxl")
                
            # plt.title('AUMic 20220826 V-C '+str(radii[aper_num])+" pxl") 
            
            plt.figure(2,figsize=(16,8))
            plt.errorbar(targ_lc[0,:],targ_lc_diff,yerr=None,fmt='r.',ecolor='black',capsize=2)
            plt.xlabel('MJD',fontsize=20)
            plt.ylabel('$\Delta m$')
            plt.gca().invert_yaxis()
            plt.title('V-C '+str(int(radii[aper_num]))+" pxl")
            # plt.figure(2,figsize=(16,8))
            # plt.plot(targ_lc[0,:],targ_lc_diff,'r')
            # plt.xlabel('MJD',fontsize=20)
            # plt.ylabel('$\Delta m$')
            # plt.gca().invert_yaxis()
            # plt.title('AUMic 20220826 V-C '+str(radii[aper_num])+" pxl")   
            # plt.xlim([59818.0040,59818.1836])
            # # plt.ylim([1.05,1.17])
            # plt.ylim([-3.22,-3.10])
            plt.show()
        while True:
             yes_or_no = input("Do you wish to plot light curve for a different aperture? Yes or no: ")
             if yes_or_no in boolean_inputs:
                 break
             else:
                 print("Invalid input") 
                 
        if yes_or_no == 'no' :
            break