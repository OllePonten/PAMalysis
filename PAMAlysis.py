# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Fri May 7 13:36:59 2021
Software for analysing tiff image stacks created by pacman
@author: Olle
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import cv2
import tifffile
import os
import csv   
import pathlib
import warnings

DEBUG = False

def load_PAM_Params(fp = "PAMSet.txt"):
    try:
        with open(fp, mode='r') as file:
            lines = file.readlines()
            lines = [line.rstrip('\n') for line in lines]
            lines = [line.split('=') for line in lines]
            params = dict(lines)
            return params
    except FileNotFoundError:
        print("No pamset text file found, proceeding with defaults")
        return []

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
    global DEBUG
    
    cv2.destroyAllWindows()
    plt.close('all')
    tifs = []
    border = 100
    create_Plots = True
    intervall = 5
    minsize = 5
    maxsize = 60
    subpopthreshold_size = 0.3
    subpopfloor = 0.1
    threshold = 25
    Hists = False
    filterMethods = {"SYD":0.2}
    if(batch):
        settings = load_PAM_Params(fp+"/PAMSet.txt")
    else:
        settings = load_PAM_Params()
    output_folder = ""
    if(len(settings) > 0):
        keys = settings.keys()
        if('minsize' in keys):
            try:
                minsize = int(settings['minsize'])
            except:
                print("Minsize badly formatted")
        if('maxsize' in keys):
            try:
                maxsize = int(settings['maxsize'])
            except:
                print("maxsize badly formatted")
        if('floor' in keys):
            try:
                subpopfloor = float(settings['floor'])
            except:
                print("floor badly formatted")
        if('subpop_size' in keys):
            try:
                subpopthreshold_size = float(settings['subpop_size'])
            except:
                print("subpopthreshold badly formatted")
        if('border' in keys):
            try:
                border = int(settings['border'])
            except:
                print("border badly formatted")
        if('threshold' in keys):
            try:
                threshold = float(settings['threshold'])
            except:
                print("Threshold badle formatted")
        if('intervall' in keys):
            try:
                intervall = int(settings['intervall'])
            except:
                print("Intervall badly formatted")
        if('SYD' in keys):
            try:
                SYD = bool(int(settings['SYD']))
                if(SYD == False):
                    try:
                        filterMethods.pop("SYD")
                    except:
                        pass
            except:
                print("SYD badly formatted")
        if('Histogram' in keys):
            try:
                Hists = bool(settings['Histogram'])
            except:
                pass
                
    outYields = dict()
    if(batch):      
        fl = os.listdir(fp)
        tifs = [fp+ "/" + file for file in fl if ".tif" in file]
        print("Running batch mode")
        print(f"Analysing following files: {str(tifs)}")
    else:
        tifs = [fp]
    try:
        output_folder = f'Output/{work_name}'
        os.makedirs(output_folder)
    except OSError:
        print("Could not create output folder")
    for fovidx, i in enumerate(tifs):
        if('/' in i):
            fn = i.split('/')[-1][:-4]
        else:
            fn = i
		try:
            tif = tifffile.imread(i)
        except:
            input(f"Could not load image: {fn} Please check filename")
            sys.exit()
        if(len(tif) > 6):
            tif = np.concatenate((tif[0:2],tif[4:]))
        else:
            tif = tif[0:2]
        tif_tags = {}
        #Get all tags
        with tifffile.TiffFile(i) as tif_file:
            for tag in tif_file.pages[0].tags.values():
                name, value = tag.name,tag.value
                tif_tags[name] = value
        imgwidth = 640
        imgheight = 480
        if(border > 0):
            yields = [frame[border:imgheight-border,border*2:imgwidth-border*2]for frame in tif]
            yields = make_Yield_Images(yields)
        else:
            yields = make_Yield_Images(tif)
        cv2.imshow("First yield image", np.asarray(yields[0]*255,dtype=np.uint8))
        yields_for_img = (yields*255).astype(dtype = np.uint8)
        tifffile.imwrite(f'{output_folder}/' + f"Yields_{work_name}_{fn}.tif",data = yields_for_img)
        mask = 0
        if("Projection" in AOI_mode):
            mask = create_Masks(yields, threshold)
            cv2.imshow("Mask",mask*255)
        cnts,hrs = cv2.findContours(mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        contimg = cv2.cvtColor(np.zeros_like(mask),cv2.COLOR_GRAY2BGR)
        drawn = cv2.drawContours(contimg,cnts,-1,(0,0,255),-1)
        cv2.imshow("Conts",drawn)
        print(f"Raw contours found: {len(cnts)}")
        #Remove all contours below a given size
        size_filt_cnts = []
        for cnt in cnts:
            mu = cv2.moments(cnt)
            if(mu['m00'] > minsize and mu['m00'] < maxsize):
                size_filt_cnts.append(cnt)
        print(f"Contours remaining filtering with size thresholds: {minsize}-{maxsize} pixels are {len(size_filt_cnts)}")
        #Draw a filtered mask
        filteredMask = cv2.drawContours(cv2.cvtColor(np.zeros_like(mask,dtype=np.uint8),cv2.COLOR_GRAY2BGR),size_filt_cnts,-1,(255,255,255),-1)
        cv2.imshow("Filtered conts",filteredMask)  
        #Output table
        meanYields = np.zeros(shape=(len(size_filt_cnts),len(yields)+1))
        for cellidx, cnt in enumerate(size_filt_cnts):
            #Get minimum bounding rect
            rect = cv2.boundingRect(cnt)
            cellMinis = yields[:,rect[1]-1:rect[1]+rect[3]+1,rect[0]-1:rect[0]+rect[2]+1]
            filteredMask = cv2.putText(filteredMask,str(cellidx), (rect[0],rect[1]+int(rect[3]/2)), cv2.FONT_HERSHEY_PLAIN, 0.5,(0,255,0),thickness = 1)
            for timeidx,img in enumerate(cellMinis):  
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        #First we make all zeros into np.nan Np.nanmean then allows us to sum
                        #only over nonnan numbers, giving an accurate average yield
                        #We do this within warnings to catch completely nan-filled areas
                        #(Likely cells who have wandered off)
                        meanYield = np.nanmean(np.where(img!=0,img,np.nan))  
                    except Warning:
                        #print(f"Nan/Zero yield enc. Timeindex: {timeidx}. Cellindex: {cellidx}. Setting zero")
                        meanYield = 0
                    if(timeidx == 0):
                        meanYields[cellidx,0] = cellidx
                meanYields[cellidx,timeidx+1] = meanYield  
        
        #THRESHOLD FILTER
        filteredYields = filter_Yields(meanYields, filterMethods)
        print(f"Yields remaining after filter: {len(filteredYields)}")
        cv2.imshow("Numbered masks",filteredMask)
        cv2.imwrite(f'{output_folder}/' + work_name + '_' + fn + "numbered masks.tif", filteredMask)       
        with open(f'{output_folder}/' + work_name +'_'+ fn + '.csv', mode = 'w',newline="") as pos_file:
            yield_writer = csv.writer(pos_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
            if(isinstance(settings,dict)):
                yield_writer.writerow(settings.items() )
            for idx,part in enumerate(filteredYields):
                yield_writer.writerow(part)     
        if(create_Plots):
            subs, names = subdivide_Yield(filteredYields[:,1:], threshold_size = subpopthreshold_size, disc_pos = 1,floor = subpopfloor)
            plot_Values(subs,names, work_name, fn,fovidx, intervall)
            plot_hists(subs)
            #plot_hists(filteredYields[:,1:])
        #For outside use, return our filtered yields
        outYields[f"{fn}"] = filteredYields
    return outYields
    
def reanalyze(yields, indexes):
    #If you want to reanalyze only specific numbered cells.
    manFilteredYields = [part for part in yields if part[0] in indexes]
    return manFilteredYields
    
def plot_Values(yields, names, jobname, filename, subjob, intervall = 5, rows = -1, columns = -1, mode = "Lines", floor = 0.2):
    #Assumes that yields is formatted as yields.shape = [n(subplots),n(samples),n(values)]
    color = None
    output_dir = f"Output/{jobname}"
    if(subjob%3==0):
        color = [0,0,np.random.rand()]
    elif(subjob%2 == 0):
        color = [0,np.random.rand(),0]
    else:
        color = [np.random.rand(),0,0]
    tot = len(names)
    if(tot == 1):
        #Only one dataset.
        columns = 1
        rows = 1
    if(rows == -1 or columns == -1):
        columns = 1
        if(tot % 2 == 0):
            columns = int(tot/2)
            
            rows = 2
        elif(tot % 2 == 1):
            columns = int(tot/2 + 1)
            rows = 2
    
    avg_lines = []
    avg_errors = []
    avg_sizes = []
    xlim = [0,0]
    try:
        xlim =[0,len(yields[0][0])*intervall]
    except:
        print("Not enough data points. Shutting down")
        return     
    ylim = [floor,0.7] 
    Position = range(1,tot + 1)
    #fig = plt.figure(figsize=(5*rows, 3*columns))
    plt.close(f"{jobname}: Subpopulations")
    fig = plt.figure(f"{jobname}: Subpopulations",figsize = (6*columns,rows*5))
    fig.suptitle("Subpopulations")
    for k in range(tot):
        # add every single subplot to the figure with a for loop
        ax = fig.add_subplot(rows,columns,Position[k])
        ax.set_title(names[k])
        ax.set_ylabel("Yield")
        ax.set_xlabel("Minutes")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        avg_line = (np.mean(yields[k][:],axis=0))
        avg_error = (np.std(yields[k][:],axis=0))
        pop_size = len(yields[k][:])
        avg_lines.append(avg_line)
        avg_errors.append(avg_error)
        avg_sizes.append(pop_size)
        for part in yields[k]:
            ax.plot(range(xlim[0],xlim[1],intervall),part, marker='o', markersize = 3, linewidth = 0.5)
    fig.tight_layout(pad = 3.0)
    fig.savefig(fname =f"{output_dir}/{jobname}_{subjob}_total_yields")
    
    #plt.close(f"{jobname}: Average_Yield")
    fig2 = plt.figure(f"{jobname}: Average_Yield")
    fig2.suptitle(f"{jobname}: Average Yield")              
    for idx, avgs in enumerate(avg_lines):
        avg_color = color + avgs[idx]
        plt.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], label=f"Sample size: {avg_sizes[idx]}", linewidth = 3, linestyle = 'dashed', c=avg_color, capsize = 5, elinewidth = 1, errorevery =(1,10))
        #plt.plot(range(xlim[0],xlim[1],intervall),avgs, label=f"{filename}: S:{avg_sizes[idx]}", c=avg_color)
        plt.ylabel("Yield")
        plt.xlabel("Minutes")
        plt.xlim(xlim)
        plt.ylim(ylim)
        #gca = plt.gca()
        #yax = gca.yaxis()
        #yax.set_major_locator(plt.MultipleLocator(0.1))
        #yax.set_major_formatter(np.arange(ylim[0], ylim[1], step=(ylim[1]-ylim[0])/10)
        #yax.set
        plt.yticks(np.arange(ylim[0], ylim[1], step=(ylim[1]-ylim[0])/10), labels=None)
    fig2.legend(bbox_to_anchor=(0.95,0.85), ncol = 2)
    fig2.tight_layout()
    fig2.savefig(fname = f"{output_dir}/{jobname}_Average_Yields")
    #1 big plot of just means
    #figure of subplots with subpopulation datapoints compared to all means 


def plot_hists(yields):
        fig3 = plt.hist(yields)
        
def filter_Yields(cellyields, meths):
    SYD = False
    sydthres = 1
    #For uniqueness, use set
    remidxs = set()
    outputmsg = ""
    if('SYD' in meths.keys()):
        #Sudden Yield Drop. Filter based on corresponding threshold
        sydthres = meths['SYD']
        SYD = True
    for idx,cell in enumerate(cellyields):
        if(SYD):
          for time in range(1,len(cell)-1):
                  if(cell[time]-cell[time+1] > sydthres):
                     #Remove, likely fell away
                     #Tail function so we don't remove on single frame brightness
                      if(time+2 <= len(cell)-2):
                          if(cell[time+3] <= cell[time]):
                              remidxs.add(idx)
    outputmsg += f"Filtered: {len(remidxs)} based on threshold: {sydthres}. {list(remidxs)}"
    print(outputmsg)
    return np.delete(cellyields,list(remidxs),axis=0)
              
def subdivide_Yield(cellyields, method = "Static Bins",floor = 0, threshold_size = 0.1, disc_pos = 1):
    #Several modes for finding subpopulations
    #Static Value = Divides up in subpopulations based on the yield value in disc pos 
    #using static bins of threshold_size
    #Distribution = Subpopulates based on distributions based on threshold size so that each
    #subpopulation contains sample_size*threshold_size samples
    #disc_pos selects what time-position to use when comparing values
    subpops = []  
    names = []
    if(method == "Static Bins"):
        subs = int((1-floor)/threshold_size)
        #We can't know the sizes of the populations from the start
        subpops = []
        for idx in range(0,subs):
            temp = []
            for cell in cellyields:
                if(floor + (idx*threshold_size)) < cell[disc_pos] and cell[disc_pos] <= (floor+((idx+1)*threshold_size)):          
                    temp.append(cell)
            #Else it is empty
            if(len(temp) > 0):
                subpops.append(temp)
                name = f"Subpopulation: {idx}, size: {len(temp)}. Static threshold: {(floor + (idx*threshold_size)):.3f}-{(floor+((idx+1)*threshold_size)):.3f}"
                print(name + f" {len(temp)}")
                names.append(name)          
    if(method == "Distribution"):
        #Create equally sized populations
        subs = len(cellyields)*threshold_size
        #We can know the sizes of the populations from the start.
        subpops = np.nan(shape = int(1/threshold_size,subs,cellyields.shape[2]))
        #Sort based on disc_pos value
        ntile_size = len(cellyields*threshold_size)
        cellyields[cellyields[:,disc_pos].argsort()]
        for idx in range(0,int(1/threshold_size)):
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
        Mask = np.where(Fo[i] > 7,1,0)
        #Emulate remove outliers from imageJ (Which is just a median filter)
        #Mask = Mask.astype(np.uint8)        
        #Mask = cv2.bilateralFilter(Mask,2,3,3)
        Fv = np.subtract(Fm[i],Fo[i],dtype = np.float32)
        #Floor to zero
        Fv = np.clip(Fv,0,255)*Mask
        cYield = np.divide(Fv,Fm[i],out=np.zeros_like(Fv),where=Fm[i]!=0)
        Yield.append(cYield)
    return np.asarray(Yield)
    
def create_Masks(imgstack, maskthres = 20):
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
    #Only keep pixels as part of a cell if they have an average of threshold pixel
    #intensity over the entire stack, i.e. the length of the stack.
    summed = summed*255
    print(np.median(summed))
    threshold = imgstack.shape[0]*maskthres
    #Simple threshold
    summed[summed < threshold] = 0
    summed = summed.astype(np.uint8)
    #summed = cv2.medianBlur(summed,3)    
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

global data
if __name__ == '__main__':
    global data
    batch_flag = False
    fp = ""
    job_name = ""
    args = sys.argv[1:]
    print("Args:" + str(args))
    if("/help" in args or "/h" in args or len(args)== 0 in args):
        print("This is the PAMalysis software which computes quantum yields of tiffstacks output"
              + " from ImagingWinGigE v2.51d.\n")
        print("Following options are available:\n")
        print("[/batch /b] Sets the analysis to be in batch mode(analyses all files in a directory). Directory must then be provided with directory flag argument")
        print("[/dir /d ] followed by DIRECTORY. string argument specificying relative path.\n ")
        print("[/file /f] Only useable if batch mode is not on. Analyses a single stack. \n")
        print("[/job /j] Sets the output folder name. \n")
        print("Made by Olle Pontén, Uppsala University, 2021")
        input("Press key to exit")
        sys.exit()
    if ("/Batch" in args or "/batch" in args or "/b" in args):
        batch_flag = True
        if("/dir" in args or "/d" in args):
            try:
                findex = args.index("/dir") + 1
            except ValueError:
                try:
                    findex = args.index("/d") + 1     
                except:
                    input("No directory flag found. Enter key to exit")
                    sys.exit()
            fp = args[findex]
        else:
            #Assume we are just grabbing all tif stacks in current folder
            fp = str(pathlib.Path().absolute())
    elif("/file" in args or "/f" in args):
        try:
            findex = args.index("/file") + 1
        except ValueError:
            try:
                findex = args.index("/f") + 1     
            except:
                input("No filename flag found. Enter key to exit")
                sys.exit()
        fp = args[findex]
    try:
        jindex = args.index("/job") + 1
    except ValueError:
        try:
            jindex = args.index("/j") + 1
        except:
            input("No jobname found. Enter key to exit")
            sys.exit()
        job_name = args[jindex]
    data = perform_Analysis(fp,job_name, batch = batch_flag)
    
def cleanup():
    import cv2
    import matplotlib.pyplot as plt
    cv2.destroyAllWindows()
    plt.close('all')
    