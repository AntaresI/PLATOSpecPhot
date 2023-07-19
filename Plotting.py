import numpy as np
import matplotlib.pyplot as plt
import math

def plotting(targ_lc, comp_lc, vali_lc, radii):
    
    while True:
        boolean_inputs = ["Y", "n"]
        num_of_apertures = len(radii)
        
        '''COMPUTING DIFFERENCES OF MAGNITUDES AND COMPUTING STANDARD DEVIATIONS'''
        print('WRITING OUT ALL STANDARD DEVIATIONS')
        
        for i in range(num_of_apertures):
            iaper = i
            targ_lc_diff = targ_lc[iaper+1,:]-comp_lc[iaper+1,:]
            vali_lc_diff = comp_lc[iaper+1,:]-vali_lc[iaper+1,:]
            err_targ = targ_lc[iaper+num_of_apertures+1,:]
            # err_comp = comp_lc[iaper+num_of_apertures+1,:]
            # err_vali = vali_lc[iaper+num_of_apertures+1,:]
            print("Aperture number ",i," radius ",radii[i],"px"," STD C-K: ",round(np.nanstd(vali_lc_diff),4),"  STD V-C: ",round(np.nanstd(targ_lc_diff),4),"  Mean error: ",round(np.nanmean(err_targ),5))
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''GETTING FROM USER WHICH APERTURE TO PLOT'''
        while True:
            which_aper_num = input("Enter the aperture number which you want to use for plotting: ")
            
            if int(which_aper_num) in range(len(radii)):
                break
            else:
                print("This number doesn't correspond to any aperture.")
                
        aper_num = int(which_aper_num)
        ''''''''''''''''''''''''''''''''''''''''''''
        
        targ_lc_diff = targ_lc[aper_num+1,:]-comp_lc[aper_num+1,:]
        vali_lc_diff = comp_lc[aper_num+1,:]-vali_lc[aper_num+1,:]
        
        '''COMPUTING ERRORS OF DIFFERENCES OF MAGNITUDES FOR THE CORRESPONDING CHOSEN APERTURE'''
        err_targ = targ_lc[aper_num+num_of_apertures+1,:]
        
        err_comp = comp_lc[aper_num+num_of_apertures+1,:]
        err_vali = vali_lc[aper_num+num_of_apertures+1,:]
        diff_lc_targ_err = np.sqrt(err_targ**2+err_comp**2)
        diff_lc_vali_err = np.sqrt(err_comp**2+err_vali**2)
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        
        '''GETTING INPUT FROM USER WHETHER TO PLOT ERRORBARS'''
        while True:
            errorbars = input("Do you wish to plot errorbars? Y/n: ")
            
            if errorbars in boolean_inputs:
                break
            else:
                print("Invalid input.")
        ''''''''''''''''''''''''''''''''''''''''''''''''''''''
                
        '''PLOTTING'''
        if errorbars == "Y" :
            # fig, axs = plt.subplots(2,gridspec_kw={'height_ratios': [3, 1]},sharex=True)
            # axs[0].scatter(times2x2,meansalllo2, s=20,color='blue')

            # popt, pcov = opt.curve_fit(lin, times2x2, meansalllo2)
            # axs[0].plot(times2x2, lin(times2x2, *popt), color='blue', linestyle='dashed')

            # axs[1].set_xlabel('Expoziční doba [s]')
            # axs[0].set_ylabel('Průměrné ADU [ADU]')
            # axs[1].set_ylabel('Rozptyl [ADU]')
            # axs[1].scatter(times2x2,stdsalllo2, s=20,color='blue')
            # axs[0].set_xticks(np.arange(0,210,30))
            # vali_lc[0,:]-math.floor(vali_lc[0,0])+2400000.5
            # print(vali_lc[0,:]-math.floor(vali_lc[0,0])+0.5)
            print(f"JD date of the first value is {vali_lc[0,0]}")
            # x = vali_lc[0,:]-math.floor(vali_lc[0,0])-0.5
            fig, axs = plt.subplots(2,gridspec_kw={'height_ratios': [3, 1]},sharex=True)
            axs[0].scatter(vali_lc[0,:]-math.floor(vali_lc[0,0]),vali_lc_diff,s=2,color='blue')
            axs[1].set_xlabel('JD - '+str(math.floor(vali_lc[0,0])),fontsize=10)
            axs[0].set_ylabel('Relative magnitude [mag]', fontsize=10)
            axs[1].set_ylabel('Error [mag]')
            axs[1].scatter(vali_lc[0,:]-math.floor(vali_lc[0,0]),diff_lc_vali_err, s=2,color='blue')
            axs[0].invert_yaxis()
            axs[0].set_title('Comparison - Validation      Aperture '+str(int(radii[aper_num]))+" pxl")
            
            figa, axos = plt.subplots(2,gridspec_kw={'height_ratios': [3, 1]},sharex=True)
            axos[0].scatter(targ_lc[0,:]-math.floor(targ_lc[0,0]),targ_lc_diff,s=2,color='red')
            axos[1].set_xlabel('JD - '+str(math.floor(targ_lc[0,0])),fontsize=10)
            axos[0].set_ylabel('Relative magnitude [mag]', fontsize=10)
            axos[1].set_ylabel('Error [mag]')
            axos[1].scatter(targ_lc[0,:]-math.floor(targ_lc[0,0]),diff_lc_targ_err, s=2,color='red')
            axos[0].invert_yaxis()
            axos[0].set_title('Visual - Comparison         Aperture '+str(int(radii[aper_num]))+" pxl")
            
            plt.show()
            
        else:
            print(f"JD date of the first value is {vali_lc[0,0]}")
            plt.figure(1,figsize=(16,8))
            plt.errorbar(vali_lc[0,:]-math.floor(vali_lc[0,0]),vali_lc_diff,yerr=None,fmt='b.',ecolor='black',capsize=2)
            plt.xlabel('JD - '+str(math.floor(vali_lc[0,0])),fontsize=10)
            plt.ylabel('Relative magnitude [mag]', fontsize=10)
            plt.gca().invert_yaxis()
            plt.title('Comparison - Validation      Aperture'+str(int(radii[aper_num]))+" pxl")
            
            plt.figure(2,figsize=(16,8))
            plt.errorbar(targ_lc[0,:]-math.floor(targ_lc[0,0]),targ_lc_diff,yerr=None,fmt='r.',ecolor='black',capsize=2)
            plt.xlabel('JD - '+str(math.floor(targ_lc[0,0])),fontsize=10)
            plt.ylabel('Relative magnitude [mag]', fontsize=10)
            plt.gca().invert_yaxis()
            plt.title('Visual - Comparison         Aperture '+str(int(radii[aper_num]))+" pxl")
            plt.show()
        ''''''''''''   
        
        while True:
             yes_or_no = input("Do you wish to plot light curve for a different aperture? Yes or no: ")
             if yes_or_no in boolean_inputs:
                 break
             else:
                 print("Invalid input") 
                 
        if yes_or_no == 'n' :
            break
        
