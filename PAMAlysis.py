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

def perform_Analysis(fp,work_name, batch = False,debug = True):
    
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
    AOI_mode = "Projection"
    tifs = []
    border = 100
    intervall = 5
    minsize = 5
    maxsize = 60
    subpopthreshold_size = 0.3
    subpopfloor = 0.1
    threshold = 20
    create_Hists = False
    create_Plots = True
    Debug=False
    globalcoordinates = False
    start_point=0
    legends = True
    end_point=-1
    filterMethods = {"SYD":0.3}
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
                create_Hists = bool(int(settings['Histogram']))
            except:
                print("Histogram badly formatted")
                pass        
        if('Plots' in keys):
            try:
                create_Plots = bool(int(settings['Plots']))
            except:
                print("Plots badly formatted")
                pass
        if('start_point' in keys):
            try:
                start_point = int(settings['start_point'])
            except:
                print("Startpoint badly formatted")
        if('end_point' in keys):
            try:
                end_point = int(settings['end_point'])
            except:
                print("Endpoint badly formatted")
        if('AOI_Mode' in keys):
            try:
                AOI_mode = str(settings['AOI_Mode'])
            except:
                print("AOI Mode unknown")
        if('legends' in keys):
            try:
                legends = bool(int(settings['legends']))
            except:
                print("AOI Mode unknown")
        #if('Sorting_Pos' in keys):
    outYields = dict()
    if(batch):      
        fl = os.listdir(fp)
        tifs = [fp+ "/" + file for file in fl if ".tif" in file]
        print("Running batch mode")
        print(f"Analysing following files: {str(tifs)}")
    else:
        tifs = [fp]
    #Clean up output folder
    output_folder = f'Output/{work_name}'
    try:
        os.makedirs(output_folder)
    except OSError:
        print("Could not create output folder")
    try:
        os.remove(f'{output_folder}/' + work_name+'AllYields.csv')
    except:
        pass
    #Start of actual analysis: Read files
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
            #Analyse everythin
            if(end_point == -1):
                end_point = len(tif)
            print(f"Analysing time points: {start_point}:{end_point}")
            #Analyse everything           
            tif = tif[4+(start_point*2):4+(end_point*2)]
        else:
            tif = tif[0:2]
        tif_tags = {}
        #Get all tags
        globalx = 0
        globaly = 0
        with tifffile.TiffFile(i) as tif_file:
            for tag in tif_file.pages[0].tags.values():
                name, value = tag.name,tag.value
                tif_tags[name] = value
            if('ImageDescription' in tif_tags.keys() and globalcoordinates):
                desc = tif_tags['ImageDescription']
                try:
                    globalx = tif_tags["xposition"]
                except:
                    xposind = desc.find("xposition")+len("xposition")+1
                    globalx = int(desc[xposind:xposind+desc[xposind:].find("\n")])
                try:
                    globaly = tif_tags["yposition"]
                except:
                    yposind = desc.find("yposition")+len("yposition")+1  
                    try:
                        globaly = int(desc[yposind:yposind+desc[yposind:].find("\n")])    
                    except:
                        globaly = int(desc[yposind:])   
        imgwidth = 640
        imgheight = 480
        if(border > 0):
            yields = [frame[border:imgheight-border,border*2:imgwidth-border*2]for frame in tif]
            yields = make_Yield_Images(yields)
        else:
            yields = make_Yield_Images(tif)
        if(Debug):
            cv2.imshow("First yield image", np.asarray(yields[0]*255,dtype=np.uint8))
        yields_for_img = (yields*255).astype(dtype = np.uint8)
        tifffile.imwrite(f'{output_folder}/' + f"Yields_{work_name}_{fn}.tif",data = yields_for_img)
        mask = 0
        if("Projection" in AOI_mode):
            mask = create_Masks(yields,threshold)  
        elif("Ft_Masks" in AOI_mode):
            mask = create_Masks_Ft([frame[border:imgheight-border,border*2:imgwidth-border*2]for frame in tif],threshold)
            #yields = yields*mask
            #cnts = []
            #for(yieldimg in yields):
            #    newcnts, hrs = cv2.findContours(yieldimg,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            #    cnts.append(newcnts)
        disp_mask = np.asarray(mask,dtype=np.uint8)
        if(Debug):
            cv2.imshow("Mask",disp_mask)
        cnts,hrs = cv2.findContours((mask*255).astype(np.uint8),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        contimg = cv2.cvtColor(np.zeros_like(mask),cv2.COLOR_GRAY2BGR)
        drawn = cv2.drawContours(contimg,cnts,-1,(0,0,255),-1)
        if(Debug):
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
        if(Debug):
            cv2.imshow("Filtered conts",filteredMask)  
        #Output table
        meanYields = np.zeros(shape=(len(size_filt_cnts),len(yields)+3))
        for cellidx, cnt in enumerate(size_filt_cnts):
            #Get minimum bounding rect
            rect = cv2.boundingRect(cnt)
            if(globalcoordinates):
                cellcenter = [globalx+(rect[1]+rect[3])/2,globaly+(rect[0]+rect[2])/2]
            else:
                cellcenter = [(rect[1]+rect[3])/2,(rect[0]+rect[2])/2]
            cellMinis = yields[:,rect[1]-1:rect[1]+rect[3]+1,rect[0]-1:rect[0]+rect[2]+1]
            filteredMask = cv2.putText(filteredMask,str(cellidx), (rect[0],rect[1]+int(rect[3]/2)), cv2.FONT_HERSHEY_PLAIN, 0.5,(0,255,0),thickness = 1)
            zeroedidx = []
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
                        if(cellidx not in zeroedidx):
                            print(f"Nan/Zero yield enc. Timeindex: {timeidx}. Cellindex: {cellidx}. Setting zero")
                            zeroedidx.append(cellidx)
                        meanYield = 0
                    if(timeidx == 0):
                        meanYields[cellidx,0] = cellidx
                        meanYields[cellidx,1] = cellcenter[0]
                        meanYields[cellidx,2] = cellcenter[1]
                meanYields[cellidx,timeidx+3] = meanYield      
        #THRESHOLD FILTER
        filteredYields = filter_Yields(meanYields[:,:], filterMethods)
        print(f"Yields remaining after filter: {len(filteredYields)}")
        cv2.imshow("Numbered masks",filteredMask)
        cv2.imwrite(f'{output_folder}/' + work_name + '_' + fn + "numbered masks.tif", filteredMask)       
        with open(f'{output_folder}/' + work_name+'AllYields.csv', mode = 'a', newline="") as tot_file:
            tot_yield_writer = csv.writer(tot_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
            with open(f'{output_folder}/' + work_name +'_'+ fn + '.csv', mode = 'w',newline="") as pos_file:
                yield_writer = csv.writer(pos_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
                yield_writer.writerow(fn)
                if(isinstance(settings,dict)):
                    yield_writer.writerow(settings.items() )
                for idx,part in enumerate(filteredYields):
                    yield_writer.writerow(part)    
                    tot_yield_writer.writerow(part)
                #filteredYields.shape[1]-2
                subs, names = subdivide_Yield(filteredYields[:,3:], threshold_size = subpopthreshold_size, disc_pos = 0,floor = subpopfloor)
        if(create_Plots):
            plot_Values(subs,names, work_name, fn,fovidx, intervall,floor=subpopfloor,legends = legends)      
        #For outside use, return our filtered yields
        outYields[f"{fn}"] = filteredYields
    if(create_Hists):
        unsorted_yields_start = []
        unsorted_yields_end = []
        for x in list(outYields.values()):
            unsorted_yields_start = np.concatenate((unsorted_yields_start, x[:,3+start_point]))
            if(start_point != end_point):                
                unsorted_yields_end = np.concatenate((unsorted_yields_end, x[:,int(end_point/2)]))
        #print(len(list(outYields.values())))
        #unsorted_Yields = np.concatenate((list(outYields.values())),axis=1)       
        plot_hists(unsorted_yields_start,work_name,subpopfloor,start_point,"red")
        if(start_point != end_point):
            plot_hists(unsorted_yields_end,work_name,subpopfloor,int((end_point/2)-1)*intervall,"green")
            
    return outYields

def reanalyze(yields, indexes):
    #If you want to reanalyze only specific numbered cells.
    manFilteredYields = [part for part in yields if part[0] in indexes]
    return manFilteredYields
    

def plot_Values(yields, names, jobname, filename, subjob, intervall = 5, rows = -1, columns = -1, mode = "Lines", floor = 0.2,legends=True):
    #Assumes that yields is formatted as yields.shape = [n(subplots),n(samples),n(values)]
    rows = -1
    columns = -1
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
    else:
    #if(rows == -1 or columns == -1):
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
        print(yields)
        print("Population empty, could not plot values. ")
        return     
    ylim = [0,1] 
    Position = range(1,tot + 1)
    #fig = plt.figure(figsize=(5*rows, 3*columns))
    plt.close(f"{jobname}: Subpopulations")
    fig = plt.figure(f"{jobname}: Subpopulations",figsize = (6*columns,rows*5))
    fig.suptitle("Subpopulations")
    for k in range(tot):
        # add every single subplot to the figure with a for loop
        ax = fig.add_subplot(rows,columns,Position[k])
        ax.set_title(names[k])
        ax.set_ylabel("Fv/Fm")
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
    fig2.suptitle(f"{jobname}: Average Fv/Fm")       
    avg_of_all = np.mean(avg_line,axis=0)       
    for idx, avgs in enumerate(avg_lines):
        if(legends):
            plt.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], label=f"{names[idx]}", markersize = 3, marker='o',linewidth = 1.5, capsize = 2, elinewidth = 1, errorevery =(1,3))
        else:
            plt.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], markersize = 3, marker='o',linewidth = 1.5, capsize = 2, elinewidth = 1, errorevery =(1,3))           
        #plt.plot(range(xlim[0],xlim[1],intervall),avgs, label=f"Sample size: {avg_sizes[idx]}", linewidth = 3, linestyle = 'dashed')
        #plt.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], label=f"Sample size: {avg_sizes[idx]}", linewidth = 3, linestyle = 'dashed', capsize = 5, elinewidth = 1, errorevery =(1,10))
        plt.plot(avg_of_all)
        plt.ylabel("Yield")
        plt.xlabel("Minutes")
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.yticks(np.arange(ylim[0], ylim[1], step=(ylim[1]-ylim[0])/10), labels=None)
        xstep = 0
        if(xlim[1] <= 100):
            xstep = 25
        else:
            xstep = round(int((xlim[1]+1)/10)/100)*100
        plt.xticks(np.arange(0,xlim[1],step=xstep))
        plt.minorticks_on()
    fig2.legend(bbox_to_anchor=(0.95,0.85), ncol = 2)
    #fig2.legend(loc="upper left", ncol = 2)
    fig2.tight_layout()
    fig2.savefig(fname = f"{output_dir}/{jobname}_Average_Yields")

def plot_hists(yields,jobname, floor=0.2,time_point=0, i_color = "red"):
    print(f"Creating histograms for {jobname}")
    output_dir = f"Output/{jobname}"
    col_fig = plt.figure(f"{jobname}")
    plt.xlabel("Fv/Fm")
    plt.ylabel("Count")
    yield_bins=np.linspace(floor,0.8,num=(round((0.8-floor)/0.05)+1))
    yields = yields[(yields>0)]
    below = len([i for i in yields if i < 0.1])
    arr = plt.hist(yields, bins=yield_bins, alpha=0.7, label = f"n: {len(yields)}. <0.1: {below}. T: {time_point} mins", color = i_color, edgecolor="black")
    col_fig.suptitle(f"{jobname}: Histogram of Fv/Fm.")   
    plt.xticks(yield_bins)
    avg=np.mean(yields)
    roof = round(max(arr[0])/100+1,0)*100
    for i in range(len((yield_bins))-1):
        if(i_color=="red"):
            plt.text(arr[1][i]+0.01,arr[0][i]+0.2,str(int(arr[0][i])),color="red")
            plt.axvline(avg,linestyle='dashed',color="red")
            plt.text(avg-0.13,roof*1.01,f"Average: {avg:.3f}",color="red")
        else:
            plt.text(arr[1][i]+0.03,arr[0][i]+0.2,str(int(arr[0][i])),color="green")
            plt.axvline(avg,linestyle='dashed',color="green")
            plt.text(avg+0.02,roof*0.96,f"{avg:.3f}",color="green")
    plt.ylim(0,roof)
    plt.legend(loc="best")
    col_fig.savefig(f"{output_dir}/{jobname}")
        
def filter_Yields(cellyields, meths):
    SYD = False
    sydthres = 1
    #For uniqueness, use set
    remidxs = set()
    outputmsg = ""
    if('SYD' in meths.keys()):
        #Sudden Yield Drop. Filter based on corresponding threshold
        print("Filtering using sudden yield drop")
        sydthres = meths['SYD']
        SYD = True
    for idx,cell in enumerate(cellyields):
        if(SYD):
          for time in range(3,len(cell)-1):
              if(cell[time]-cell[time+1] > sydthres):
                  #Remove, likely fell away
                  #Tail function so we don't remove on single frame brightness
                  if(time+2 <= len(cell)-2):
                      if(cell[time+3] <= cell[time]):
                          remidxs.add(idx)
                  else:
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
                if(floor + (idx*threshold_size)) <= cell[disc_pos] and cell[disc_pos] < (floor+((idx+1)*threshold_size)):          
                    temp.append(cell)
            #Else it is empty
            if(len(temp) > 0):
                subpops.append(temp)
                name = f"{idx}. n: {len(temp)}. Threshold: {(floor + (idx*threshold_size)):.3f}-{(floor+((idx+1)*threshold_size)):.3f}"
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
        Mask = Mask.astype(np.uint8)        
        Mask = cv2.medianBlur(Mask,3)
        Mask = np.where(Mask>0,1,0)
        Fv = np.subtract(Fm[i],Fo[i],dtype = np.float32)
        #Floor to zero
        Fv = np.clip(Fv,0,255)*Mask
        cYield = np.divide(Fv,Fm[i],out=np.zeros_like(Fv),where=Fm[i]!=0)
        Yield.append(cYield)
    return np.asarray(Yield)
    
def create_Masks(imgstack, maskthres = 0.1):
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
    #summed = (summed*255)/imgstack.shape[0]
    summed = (summed)/imgstack.shape[0]
    threshold = maskthres
    #Simple threshold
    summed[summed < threshold] = 0
    #summed = cv2.medianBlur(summed,3)    
    #np.clip(summed,0,255,out=summed)
    summed[summed > 0] = 1
    summed = summed.astype(np.uint8)
    return summed    

def create_Masks_Ft(imgstack,maskthres=0.1):
    """
    Creates masks based on method from yield macro. Based on Ft values

    Parameters
    ----------
    imgstack : TYPE
        DESCRIPTION.
    maskthres : TYPE, optional
        DESCRIPTION. The default is 10.

    Returns
    -------
    None.

    """
    Fo = imgstack[::2]
    th = []
    threshold = maskthres*255
    for i in range(len(Fo)):
        src = np.array(Fo[i])
        ret,img=cv2.threshold(src,int(threshold),1, cv2.THRESH_BINARY)
        th.append(cv2.medianBlur(img,3))
    mask = np.sum(th, axis=0)
    mask[mask>0]=1
    mask = mask.astype(np.uint8)
    return mask
    
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
    