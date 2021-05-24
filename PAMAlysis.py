# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Fri May  7 13:36:59 2021
Software for analysing tiff image stacks created by pacman
@author: Olle
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import csv   
import pathlib



def perform_Analysis(fp,work_name, batch = False,debug = True,AOI_mode = "Projection"):
    """
    

    Parameters
    ----------
    fp : str
        Filepath.
    directory : bool, optional
        Determines whether the analysis is a batch(then fp should point to a directory
        containing tiff-stacks). The default is False.

    Returns
    -------
    None.

    """
    cv2.destroyAllWindows()
    tifs = []
    create_Plots = False
    minsize = 25
    maxsize = 125
    if(batch):
        fl = os.listdir(fp)
        tifs = [fp+ "/" + file for file in fl if ".tif" in file]
    else:
        tifs = [fp]
    try:
        os.makedirs('Output')
    except OSError:
        print("Could not create output folder")
    for idx, i in enumerate(tifs):
        if('/' in i):
            fn = i.split('/')[-1][:-3]
        else:
            fn = i
        tif = cv2.imreadmulti(i)[1]
        #Remove first 4 images         
        yields = make_Yield_Images(tif[4:])
        cv2.imshow("Random yield image", np.asarray(yields[np.random.randint(low=1,high=len(yields))]*255,dtype=np.uint8))
        mask = 0
        if("Projection" in AOI_mode):
            mask = create_Masks(yields)
            cv2.imshow("Mask",mask)
        cnts,hrs = cv2.findContours(mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        contimg = cv2.cvtColor(np.zeros_like(mask),cv2.COLOR_GRAY2BGR)
        drawn = cv2.drawContours(contimg,cnts,-1,(0,0,255),-1)
        cv2.imshow("Conts",drawn)
        #Remove all contours below a given size
        size_filt_cnts = []
        for cnt in cnts:
            mu = cv2.moments(cnt)
            if(mu['m00'] > minsize and mu['m00'] < maxsize):
                size_filt_cnts.append(cnt)
        #Draw a filtered mask
        filtered_mask = cv2.drawContours(cv2.cvtColor(np.zeros_like(mask,dtype=np.uint8),cv2.COLOR_GRAY2BGR),size_filt_cnts,-1,(255,255,255),-1)
        cv2.imshow("Filtered conts",filtered_mask)
        #Output table
        mean_yields = np.zeros(shape=(len(size_filt_cnts),len(yields)+1))
        for cellidx, cnt in enumerate(size_filt_cnts):
            #Get minimum bounding rect
            rect = cv2.boundingRect(cnt)
            cell_minis = yields[:,rect[1]:rect[1]+rect[3],rect[0]:rect[0]+rect[2]]
            filtered_mask = cv2.putText(filtered_mask,str(cellidx), (rect[0],rect[1]+int(rect[3]/2)), cv2.FONT_HERSHEY_PLAIN, 0.5,(0,255,0),thickness = 1)
            for timeidx,img in enumerate(cell_minis):                   
                mean_yield = np.nanmean(np.where(img!=0,img,np.nan))   
                if(timeidx == 0):
                    mean_yields[cellidx,0] = int(cellidx)
                mean_yields[cellidx,timeidx+1] = mean_yield  
        cv2.imshow("Numbered masks",filtered_mask)
        cv2.imwrite('Output/' + fn + "numbered masks.tif", filtered_mask)       
        #cv2.imshow("Cell_mini", cell_minis[np.random.randint(low=0,high=len(cell_minis))])
        with open('Output/' + work_name +'_'+ fn + 'csv', mode = 'w',newline="") as pos_file:
            yield_writer = csv.writer(pos_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
            for idx,part in enumerate(mean_yields):
                yield_writer.writerow(part)
            
        if(create_Plots):
            subs, names = subdivide_Yield(mean_yields[:,1:])
            plot_Values(subs,names,[0.0])
  
def plot_Values(yields, names, ylim, intervall = 300, rows = -1, columns = -1):
    #Assumes that yields is formatted as yields.shape = [n(subplots),n(samples),n(values)]
    if(rows == -1 or columns == -1):
        rows = int(len(yields)/2)
        columns = int(len(yields)/rows)
        
    fig,axs = plt.subplots(rows,columns, sharey = 'all',sharex='all')
    #axs[:].set_ylim(ylim)
    i = 0
    for row in axs:
        for col in row:
            print(len(yields[i]))
            if(len(yields[i]) > 1):
                print(yields[i][:][0])
                avg_line = np.mean(yields[i][:][0],axis=0)
                col.plot(range(0,len(yields[i][0])*intervall,intervall),avg_line)
                col.title.set_text(names[i])
            i += 1
              
def subdivide_Yield(cellyields, method = "Static Bins",threshold_size = 0.1, disc_pos = 1):
    #Several modes for finding subpopulations
    #Static Value = Divides up in subpopulations based on the yield value in disc pos 
    #using static bins of threshold_size
    #Distribution = Subpopulates based on distributions based on threshold size so that each
    #subpopulation contains sample_size*threshold_size samples
    #disc_pos selects what time-position to use when comparing values
    subpops = []  
    names = []
    if(method == "Static Bins"):
        subs = int(1/threshold_size)
        #We can't know the sizes of the populations from the start
        #subpops = np.zeros(shape = (subs,cellyields.shape[0],cellyields.shape[1]))
        subpops = []
        for idx in range(0,subs):
            temp = []
            for cell in cellyields:
                if(idx*threshold_size < cell[disc_pos] and cell[disc_pos] <= (idx+1)*threshold_size):          
                    temp.append(cell)
            subpops.append(temp)
            names.append(f"Subpopulation: {idx}. Static threshold: {(idx*threshold_size):.3f}-{((idx+1)*threshold_size):.3f}")          
    if(method == "Distribution"):
        #Create equally sized populations
        subs = len(cellyields)*threshold_size
        #We can know the sizes of the populations from the start.
        subpops = np.nan(shape = (1/threshold_size,subs,cellyields.shape[2]))
        #Sort based on disc_pos value
        ntile_size = len(cellyields*threshold_size)
        cellyields[cellyields[:,disc_pos].argsort()]
        for idx in range(0,10):
            #Grab percentile between idx to idx + threshold            
            subpops[idx] = cellyields[idx*ntile_size:(idx+1)*ntile_size]
            names.append(f"Subpopulation: {idx}. Percentage range: {idx*10}-{(idx+1)*10}")
           
    return subpops, names
    
    
def make_Yield_Images(img_stack):
    """
    Parameters
    ----------
    img_stack : Numpy, shape = 2x, 480, 640
        Input image stack of only Fo/Fm.

    Returns
    -------
    Yield image stack: x,480,640. Range of values are 0.0 to 1.0

    """
    #Assumes the initial 4 images have been removed
    Fo = img_stack[::2]
    Fm = img_stack[1::2]
    #Yield is defined as Fv/Fm or (Fm-Fo)/Fm
    Yield = []
    for i in range(len(Fo)):
        Mask = np.where(Fo[i] > 8,1,0)
        Fv = np.subtract(Fm[i],Fo[i],dtype = np.float32)
        #Floor to zero
        Fv = np.clip(Fv,0,255)*Mask
        cYield = np.divide(Fv,Fm[i],out=np.zeros_like(Fv),where=Fm[i]!=0)
        Yield.append(cYield)
    return np.asarray(Yield)
    
def create_Masks(imgstack):
    """
    Creates masks based on z-projection of yield values and a thresholding operation

    Parameters
    ----------
    imgstack : should be a MxNxP stack of yield images. Threshold is based on P
        DESCRIPTION.

    Returns
    -------
    Mask image.

    """
    summed = np.sum(imgstack,axis = 0) 
    #summed = np.asarray(summed,dtype=np.uint32)
    #Only keep pixels as part of a cell if they have an average of 1 pixel
    #intensity over the entire stack, i.e. the length of the stack.
    threshold = imgstack.shape[0]/20
    #Simple threshold
    summed[summed <= threshold] = 0
    summed = summed.astype(np.uint8)
    np.clip(summed,0,255,out=summed)
    return summed
    
    
def extract_Grid():
    #Extracting grid
    #Return binary mask of interest
    return 0
    
def template_matching():
    #For using a template image of the microwell instead
    #Return binary mask of interest
    return 0


if __name__ == '__main__':
    batch_flag = False
    fp = ""
    job_name = ""
    args = sys.argv[1:]
    if("/help" or "/h" in args):
        print("This is the PAMalysis software which computes quantum yields of tiffstacks output"
              + " from ImagingWinGigE v2.51d.\n")
        print("Following options are available:\n")
        print("[/batch /b] Sets the analysis to be in batch mode(analyses all files in a directory). Directory must then be provided with directory flag argument")
        print("[/dir /d ] followed by DIRECTORY. string argument specificying relative path.\n ")
        print("[/file /f] Only useable if batch mode is not on. Analyses a single stack. \n")
        print("[/job /j] Sets the output folder name. \n")
        print("Made by Olle Pontén, Uppsala University, 2021")
    if ("/Batch" in args or "/batch" in args or "/b" in args):
        batch_flag = True
        if("/dir" or "/d" in args):
            try:
                findex = args.index("/dir") + 1
            except ValueError:
                try:
                    findex = args.index("/d") + 1     
                except:
                    input("No directory flag found. Enter key to exit")
                    exit()
            fp = args[findex]
        else:
            #Assume we are just grabbing all tif stacks in current folder
            fp = str(pathlib.Path().absolute())
    else:
         if("/file" or "/f" in args):
            try:
                findex = args.index("/file") + 1
            except ValueError:
                try:
                    findex = args.index("/f") + 1     
                except:
                    input("No filename flag found. Enter key to exit")
                    exit()
            fp = args[findex]
    try:
        jindex = args.index("/job") + 1
    except ValueError:
        try:
            jindex = args.index("/j") + 1
        except:
            input("No jobname found. Enter key to exit")
            sys.quit()
        job_name = args[jindex]
    perform_Analysis(fp,job_name, batch = batch_flag)