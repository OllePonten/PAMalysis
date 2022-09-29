# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Fri May 7 13:36:59 2021
Software for analysing tiff image stacks created by pacman
@author: Olle
"""
import sys
import datetime
import tifffile
import os
import csv   
import warnings
import argparse
import pprint
import random

import numpy as np
import matplotlib.pyplot as plt
import cv2

#CALL SIG:
#runfile('C:/Users/ollpo511/Documents/GitHub/PAMalysis/PAMAlysis.py', wdir='C:Users/ollpo511/Documents', args = '{Project_Name} {PAMSet filename}')


global filenames, DEBUG
DEBUG = False

def load_PAM_Params(fp = "PAMset.txt"):
    try:
        with open(fp, mode='r') as file:
            lines = file.readlines()
            lines = [line for line in lines if(line[0] != '[' and line[-1] != ']')]
            lines = [line for line in lines if(line[0] != "#")]
            lines = [line.rstrip('\n') for line in lines]
            lines = [line.split('=') for line in lines]
            params = dict(lines)
            return params
    except FileNotFoundError:
        print("No pamset text file found, proceeding with defaults")
        return []

def perform_Analysis(fp,work_name, job_folder, batch = False, pamset=None):
    
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
    global PlotAvg
    cv2.destroyAllWindows()
    plt.close('all')
    #### DEFAULT SETTINGS ####
    AOI_mode = "Projection"
    tifs = []
    border = 75
    intervall = 5
    minsize = 5
    maxsize = 60
    subpopthreshold_size = 0.9
    subpopfloor = 0.1
    threshold = 0.048
    create_Hists = False
    create_Plots = True
    globalcoordinates = False
    hist_start=1
    start_point=0
    legends = True
    errorbars = True
    end_point=-1
    hist_end=-1
    filterMethods = {"SYD":0.3}
    cell_mask_fp=""
    sorting_meth = "Static_Bins"
    sorting_pos=1
    #########################
    settings = load_PAM_Params(pamset)
    output_folder = ""
    filenames = []
    PlotAvg = False
    cell_limit=2000
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
                if(subpopfloor + subpopthreshold_size > 1):
                    subpopthreshold_size = float(1 - subpopfloor)
            except:
                print("subpopthreshold badly formatted")
        if('sorting_method' in keys):
            try:
                sorting_meth = str(settings['sorting_method'])
            except:
                print("Unknown sorting method")
        if('Sorting_Pos' in keys):
            try:
                sorting_pos = int(settings['Sorting_Pos'])
            except:
                print("Unknown sorting ")
        if('border' in keys):
            try:
                border = int(settings['border'])
            except:
                print("border badly formatted")
        if('threshold' in keys):
            try:
                threshold = float(settings['threshold'])
            except:
                print("Threshold badly formatted")
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
        if('Hist_end' in keys):
            try:
                hist_end = int(settings['Hist_end'])
            except:
                print("Hist end point badly formatted")
                pass
        if('Hist_start' in keys):
            try:
                hist_start = int(settings['Hist_start'])
            except:
                print("Hist end point badly formatted")
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
        if('errorbars' in keys):
            try:
                errorbars = bool(int(settings['errorbars']))
            except:
                print("Bad format")
        if(AOI_mode=="Cell_mask"):
            try:                
                cell_mask_fp=job_folder + "/" + settings['Cell_mask_fp']
                print(f"Cell_mask_fp: {cell_mask_fp}")
            except:
                print("Cell map fp badly formatted")
        if('global_coordinates' in keys):
            try:
                globalcoordinates=bool(int(settings['global_coordinates']))
            except:
                print("Global coordinates badly formatted")
        if('Max_Cells' in keys):
            try:
                cell_limit=int(settings['Max_Cells'])
            except:
                print("Max cells badly formatted")
        pprint.pprint(settings)
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
    for fovidx, current_tif in enumerate(tifs):
        if('/' in current_tif):
            fn = current_tif.split('/')[-1][:-4]
        else:
            fn = current_tif
        try:
            tif = tifffile.imread(current_tif)
            print(f"Analysing {fn}.tif")
        except:
            input(f"Could not load image: {fn} Please check filename")
            sys.exit()
        print(f"Read in file: {current_tif} with {int((len(tif)-6)/2)} time points")
        if(len(tif) > 6):
            #Analyse everything
            if(end_point == -1):
                end_point = len(tif)-4
            print(f"Analysing slices: {4+start_point*2}:{4+end_point*2}")
            #Assume first are Black/NIR/Start images.     
            tif = tif[4+(start_point*2):4+(end_point*2)]
        else:
            tif = tif[0:2]
        tif_tags = {}
        #Get all tags
        globalx = 0
        globaly = 0
        with tifffile.TiffFile(current_tif) as tif_file:
            for tag in tif_file.pages[0].tags.values():
                name, value = tag.name,tag.value
                tif_tags[name] = value
            if('ImageDescription' in tif_tags.keys() and globalcoordinates):
                desc = tif_tags['ImageDescription']
                try:
                    globalx = tif_tags["xposition"]
                except:
                    xposind = desc.find("xposition")+len("xposition")+1
                    globalx = desc[xposind:xposind+desc[xposind:].find("\n")]
                    globalx = int(''.join([i for i in globalx if i.isdigit()]))
                try:
                    globaly = tif_tags["yposition"]
                except:
                    yposind = desc.find("yposition")+len("yposition")+1  
                    try:
                        globaly = desc[yposind:yposind+desc[yposind:].find("\n")]   
                        globaly = int(''.join([i for i in globaly if i.isdigit()]))
                    except:
                        globaly = int(desc[yposind:])   
        imgwidth = 640
        imgheight = 480
        if(border > 0):
            yields = [frame[border:imgheight-border,border*2:imgwidth-border*2]for frame in tif]
            yields = make_Yield_Images(yields)
        else:
            yields = make_Yield_Images(tif)
        if(DEBUG):
            cv2.imshow("First yield image", np.asarray(yields[0],dtype=np.uint8))
        yields_for_img = (yields*255).astype(dtype = np.uint8)
        tifffile.imwrite(f'{output_folder}/' + f"Yields_{work_name}_{fn}.tif",data = yields_for_img)
        mask = 0
        if("Projection" in AOI_mode):
            mask = create_Masks(yields,threshold)  
        elif("Ft_Masks" in AOI_mode):
            mask = create_Masks_Ft([frame[border:imgheight-border,border*2:imgwidth-border*2]for frame in tif],threshold)
            yields = np.multiply(yields,mask)
        elif("Cell_mask" in AOI_mode):
            print(f"Reading: {cell_mask_fp} as cell mask image")
            cell_mask = tifffile.imread(cell_mask_fp,)
            mask = cell_mask
        mask = mask.astype(dtype=np.uint8)
        if(DEBUG):
            cv2.imshow("Mask",mask)
        cnts,hrs = cv2.findContours(mask,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        contimg = cv2.cvtColor(np.zeros_like(mask),cv2.COLOR_GRAY2BGR)
        drawn = cv2.drawContours(contimg,cnts,-1,(0,0,255),-1)
        if(DEBUG):
            cv2.imshow("Conts",drawn)
        print(f"Raw contours found: {len(cnts)}")
        #Remove all contours below a given size
        size_filt_cnts = []
        for cnt in cnts:
            mu = cv2.moments(cnt)
            if(mu['m00'] > minsize and mu['m00'] < maxsize):
                x,y,w,h = cv2.boundingRect(cnt)
                maxY = imgheight-border
                maxX = imgwidth-border*2
                minX = border*2
                minY = border
                if(x+border*2 > minX and x+w+border*2 < maxX and y+border > minY and y+h+border < maxY):  
                    if(len(size_filt_cnts) >= cell_limit):
                        size_filt_cnts.pop(random.randint(0,len(size_filt_cnts)-1))
                    size_filt_cnts.append(cnt)
        print(f"Contours remaining filtering with size thresholds: {minsize}-{maxsize} pixels are {len(size_filt_cnts)}")
        #Draw a filtered mask
        filteredMask = cv2.drawContours(cv2.cvtColor(np.zeros_like(mask,dtype=np.uint8),cv2.COLOR_GRAY2BGR),size_filt_cnts,-1,(255,255,255),-1)
        cell_mask_out = cv2.cvtColor(filteredMask,cv2.COLOR_BGR2GRAY)
        tifffile.imwrite(f'{output_folder}/' + work_name + '_' + fn + "cell_mask.tif", cell_mask_out)
        if(DEBUG):
            cv2.imshow("Filtered conts",filteredMask)  
        #Output table
        meanYields = np.zeros(shape=(len(size_filt_cnts),len(yields)+3))
        numberedMask = cv2.resize(filteredMask,(filteredMask.shape[1] *2, filteredMask.shape[0] * 2),interpolation = cv2.INTER_CUBIC)
        for cellidx, cnt in enumerate(size_filt_cnts):

            #Get minimum bounding rect
            rect = cv2.boundingRect(cnt)
            if(globalcoordinates):
                cellcenter = [border + globalx+int((rect[0]+rect[2])/2)*4,border*2 + globaly+int((rect[1]+rect[3])/2)*4]
            else:
                cellcenter = [border + int((rect[0]+rect[2])/2)*4,border*2 + int((rect[1]+rect[3])/2)*4]
            cellMinis = yields[:,rect[1]-1:rect[1]+rect[3]+1,rect[0]-1:rect[0]+rect[2]+1]         
            numberedMask = cv2.putText(numberedMask,str(cellidx), (rect[0]*2,(rect[1]+int(rect[3]/2))*2), cv2.FONT_HERSHEY_PLAIN, 1,(0,255,0),thickness = 1)
            zeroedidx = []
            for timeidx,img in enumerate(cellMinis):  
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        #First we make all zeros into np.nan Np.nanmean then allows us to sum
                        #only over nonnan numbers, giving an accurate average yield
                        #We do this within warnings to catch completely nan-filled areas
                        #(Likely cells who have wandered off)
                        #We round to 3 digits because the base data is 8-bit. We should not create data where there could be none.
                        meanYield = round(np.nanmean(np.where(img!=0,img,np.nan)),3)
                    except Warning:
                        if(cellidx not in zeroedidx):
                            if(DEBUG):
                                cv2.imshow("Cell_Mini",np.asarray(img*255,np.uint8))
                            print(f"Nan/Zero yield enc. Timeindex: {timeidx+start_point}. Cellindex: {cellidx}. Setting zero")
                            zeroedidx.append(cellidx)
                        meanYield = 0
                    if(timeidx == 0):
                        #The first digit tells what field of view a cell is from
                        meanYields[cellidx,0] = (fovidx*1000)+cellidx
                        meanYields[cellidx,1] = cellcenter[0]
                        meanYields[cellidx,2] = cellcenter[1]
                meanYields[cellidx,timeidx+3] = meanYield     
                meanYields[cellidx,timeidx+3] = round(meanYields[cellidx,timeidx+3],3)
        #THRESHOLD FILTER
        filteredYields = filter_Yields(meanYields[:,:], filterMethods)
        print(f"Yields remaining after filter: {len(filteredYields)}")
        if(len(filteredYields)==0):
            print("No yields remaining. Ending analysis of current file")
            continue
        if(DEBUG):
            cv2.imshow("Numbered masks",numberedMask)
        cv2.imwrite(f'{output_folder}/' + work_name + '_' + fn + "numbered_masks.tif", numberedMask) 
        subs, names,sortedYields = subdivide_Yield(filteredYields, method = sorting_meth, threshold_size = subpopthreshold_size, disc_pos = sorting_pos,floor = subpopfloor)
        with open(f'{output_folder}/' + work_name+'AllYields.csv', mode = 'a', newline="") as tot_file:
            tot_yield_writer = csv.writer(tot_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
            times = ["Index","XPosition","YPosition"] + list(range(0,len(filteredYields[0]-3)*(intervall),intervall))
            if(isinstance(settings,dict) and fovidx==0):
                tot_yield_writer.writerow(list(settings.items()) + [str(datetime.date.today())])
                tot_yield_writer.writerow(times)
            with open(f'{output_folder}/' + work_name +'_'+ fn + '.csv', mode = 'w',newline="") as pos_file:
                yield_writer = csv.writer(pos_file, delimiter = ",",quotechar = '"', quoting = csv.QUOTE_MINIMAL)
                yield_writer.writerow(fn)
                if(isinstance(settings,dict)):
                    yield_writer.writerow(settings.items() )
                    yield_writer.writerow(times)
                for idx,part in enumerate(sortedYields):
                    yield_writer.writerow(part)    
                    tot_yield_writer.writerow(part)    
        if(create_Plots):
            plot_Values(subs,names, work_name, fn,fovidx, intervall,floor=subpopfloor,legends = legends,errorbars=errorbars)      
        #For outside use, return our filtered yields
        outYields[f"{fn}"] = filteredYields
    if(create_Hists):
        unsorted_yields_start = []
        unsorted_yields_end = []
        for x in list(outYields.values()):
            unsorted_yields_start = np.concatenate((unsorted_yields_start, x[:,2+hist_start]))
            if(hist_end != -1):            
                unsorted_yields_end = np.concatenate((unsorted_yields_end, x[:,2+hist_end-start_point]))
        #print(len(list(outYields.values())))
        #unsorted_Yields = np.concatenate((list(outYields.values())),axis=1)       
        plot_hists(unsorted_yields_start,work_name,subpopfloor,(hist_start-1)*intervall,"blue")
        if(hist_end != -1):
            plot_hists(unsorted_yields_end,work_name,subpopfloor,(hist_end-start_point)*intervall,"green")
            
    PlotAvg = False
    return outYields

def reanalyze(yields, indexes):
    #If you want to reanalyze only specific numbered cells.
    manFilteredYields = [part for part in yields if part[0] in indexes]
    return manFilteredYields
    

def plot_Values(yields, names, jobname, filename, subjob, intervall = 5, rows = -1, columns = -1, mode = "Lines", floor = 0.2,legends=True, errorbars=True):
    #Assumes that yields is formatted as yields.shape = [n(subplots),n(samples),n(values)]
    global PlotAvg
    rows = -1
    columns = -1
    color = None
    base = 3
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
    xstep = 0
    
    try:
        xlim =[0,(len(yields[0][0])-base)*intervall]
    except:
        print(yields)
        print("Population empty, could not plot values. ")
        return     
    if(xlim[1] <= 100):
        xstep = 25
    else:
        xstep = round((xlim[1]+1)/100)*10
    ylim = [0,0.7] 
    Position = range(1,tot + 1)
    #fig = plt.figure(figsize=(5*rows, 3*columns))
    plt.close(f"{jobname}: Subpopulations")
    fig = plt.figure(f"{filename}: Subpopulations",figsize = (6*columns,5*rows))
    fig.suptitle("Subpopulations")
    for k,subyields in enumerate(yields):
        # add every single subplot to the figure with a for loop
        ax = fig.add_subplot(rows,columns,Position[k])
        ax.set_title(names[k])
        ax.set_ylabel("$F_{V}$/$F_{m}$")
        ax.set_xlabel("Minutes")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        subyields = [trim_yield[3:] for trim_yield in subyields]
        avg_line = (np.mean(subyields,axis=0))
        avg_error = (np.std(subyields,axis=0))
        pop_size = len(subyields)
        avg_lines.append(avg_line)
        avg_errors.append(avg_error)
        avg_sizes.append(pop_size)
        for part in subyields:
            ax.plot(range(xlim[0],xlim[1],intervall),part, marker='o', markersize = 3, linewidth = 0.5)
        plt.yticks(np.arange(ylim[0], ylim[1], step=(ylim[1]-ylim[0])/7), labels=None)
        plt.xticks(np.arange(0,xlim[1],step=xstep))
        ax.plot(range(xlim[0],xlim[1],intervall),avg_line, marker='o', markersize = 3, linewidth = 6,linestyle='dashed',color="black")
        plt.minorticks_on()
        plt.grid(axis="y")

    #fig.tight_layout(pad = 3.0)
    fig.savefig(fname =f"{output_dir}/{jobname}_{subjob}_total_yields")
    
    
    fig2 = plt.figure(f"{jobname}: Average_Yield",figsize = [6,5.5])
    fig2.suptitle(f"{jobname}: Average "+"$F_{V}$/$F_{m}$.")
    ax = plt.subplot(111)
    for idx, avgs in enumerate(avg_lines):
        if(legends and errorbars):
            ax.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], label=f"{filename}, " + f"n(cells): {str().join([s for s in names[idx][3:12] if s.isdigit()])}", markersize = 3, marker='o',linewidth = 2, capsize = 2, elinewidth = 1, errorevery =(1,5))
        elif(errorbars):
            ax.errorbar(range(xlim[0],xlim[1],intervall),avgs, yerr = avg_errors[idx], markersize = 3, marker='o',linewidth = 2, capsize = 2, elinewidth = 1, errorevery =(1,5))
        else:
            ax.plot(range(xlim[0],xlim[1],intervall),avgs, markersize = 3, marker='o',linewidth = 4, color="black")
            #plt.errorbar(range(xlim[0],xlim[1],intervall),avgs, markersize = 3, marker='o',linewidth = 2, capsize = 2, elinewidth = 1, errorevery =(1,3),color="black")
        #plt.plot(range(xlim[0],xlim[1],intervall),avgs, label=f"Sample size: {avg_sizes[idx]}", linewidth = 3, linestyle = 'dashed')
        
        plt.ylabel("$F_{V}$/$F_{m}$")
        plt.xlabel("Minutes")
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.yticks(np.arange(ylim[0], ylim[1], step=(ylim[1]-ylim[0])/7), labels=None)
        plt.xticks(np.arange(0,xlim[1],step=xstep))
        plt.minorticks_on()
        plt.grid(True, axis="y")
    
    lgd = fig2.legend(loc = 'upper center', ncol = 2)

    if(not PlotAvg):
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height * 0.9])
        PlotAvg = True
    fig2.savefig(f"{output_dir}/{jobname}_Average_Yields", bbox_extra_artists=(lgd,), bbox_inches='tight')
    #fig2.legend(loc="upper left", ncol = 2)
    #fig2.tight_layout()
    
    
def plot_hists(yields,jobname, floor=0.2,time_point=0, i_color = "red"):
    print(f"Creating histograms for {jobname}")
    output_dir = f"Output/{jobname}"
    col_fig = plt.figure(f"{jobname}")
    plt.xlabel("$F_{V}$/$F_{m}$")
    plt.ylabel("Count")
    yield_bins=np.linspace(floor,0.8,num=(round((0.8-floor)/0.05)+1))
    below = len([i for i in yields if i <= floor])
    yields = yields[(yields>=floor)]
    avg=np.mean(yields)
    arr = plt.hist(yields, bins=yield_bins, alpha=0.7, label = f"n: {len(yields)}. {below} <= {floor}. T: {time_point} mins. Mean: {avg:.3f}", color = i_color, edgecolor="black")
    plt.xticks(yield_bins)
    roof = round(max(arr[0])/100+1,0)*100
    #Make sure roof stays the same to keep scaling.
    if(roof < plt.ylim()[1]):
        roof = plt.ylim()[1]
    plt.ylim(0,roof)
    for i in range(len((yield_bins))-1):
        if(i_color=="blue"):
            plt.text(arr[1][i]+0.01,arr[0][i]+0.2,str(int(arr[0][i])),color="blue")
            plt.axvline(avg,linestyle='dashed',color="blue")
            #plt.text(max(yield_bins)-0.12,plt.ylim()[1]*1.02,f"Mean: {avg:.3f}",color="blue")
        else:
            plt.text(arr[1][i]+0.03,arr[0][i]+0.2,str(int(arr[0][i])),color=i_color)
            plt.axvline(avg,linestyle='dashed',color=i_color)
            #plt.text(max(yield_bins)-0.12,plt.ylim()[1]*0.96,f"Mean: {avg:.3f}",color=i_color)
    

    if(len(col_fig.axes[0].get_lines())>20):
        lgd = col_fig.legend(bbox_to_anchor=(0.45,1))
    col_fig.savefig(f"{output_dir}/{jobname}", bbox_inches='tight')
        
def filter_conts(cnts,distance):    
    """
    Filters contours so that they are at least distance pixels away from each other  

    Parameters
    ----------
    cnts : TYPE
        DESCRIPTION.
    distance : TYPE
        Distance in pixel from edge to edge of minimum enclosing circles.

    Returns
    -------
    None.

    """    
    cps = []
    radii = []
    discardlist = []
    for idx,cnt in enumerate(cnts):
        (x,y),radius = cv2.minEnclosingCirlce(cnt)
        cps.append[int(x),int(y)]
        radii.append(radius)
    for idx, circle in enumerate(radii):
        for jdx, compCircle in enumerate(radii):
            #Check for self-intersection
            if(cps[jdx] != cps[idx]):
                cpdist = np.sqrt((cps[idx][0]-cps[idx][1])**2 + (cps[jdx][0]-cps[jdx][1])**2)
                if (cpdist+radii[idx]+radii[jdx] <= distance):
                    #Discard
                    discardlist.append(jdx)
    for cnt in discardlist:
        cnts.remove(cnts)
    return cnts


    
def filter_Yields(cellyields, meths):
    SYD = False
    sydthres = 1
    #For uniqueness, use set
    remidxs = set()
    tail = 5
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
                  if(time+tail <= len(cell)-1):
                      if(cell[time+tail] == 0):
                          remidxs.add(idx)
                  else:
                    remidxs.add(idx)
    if(len(list(remidxs))>0):
        outputmsg += f"Filtered: {len(remidxs)} based on an sudden decrease of more than: {sydthres}. {list(remidxs)}"
        print(outputmsg)
    return np.delete(cellyields,list(remidxs),axis=0)
              
def subdivide_Yield(cellyields, method = "Static_Bins",floor = 0, threshold_size = 0.25, disc_pos = 1):
    #Several modes for finding subpopulations
    #Static Value = Divides up in subpopulations based on the yield value in disc pos 
    #using static bins of threshold_size
    #Distribution = Subpopulates based on distributions based on threshold size so that each
    #subpopulation contains sample_size*threshold_size samples
    #disc_pos selects what time-position to use when comparing values
    sortedYields = []  
    names = []
    base=2
    if(method == "Static_Bins"):
        subs = int((1-floor)/threshold_size)
        #We can't know the sizes of the populations from the start
        subpops = []
        for idx in range(0,subs):
            temp = []
            for cell in cellyields:
                if(floor + (idx*threshold_size)) <= cell[disc_pos+base] and cell[disc_pos+base] < (floor+((idx+1)*threshold_size)):          
                    temp.append(cell)
                    sortedYields.append(cell)
            #Else it is empty
            if(len(temp) > 0):
                subpops.append(temp)
                name = f"{idx}. n: {len(temp)}. Threshold: {(floor + (idx*threshold_size)):.3f}-{(floor+((idx+1)*threshold_size)):.3f}"
                names.append(name)          
    elif(method == "Distribution"):
        print(f"Sorting based on percentile. Sorting position: {disc_pos}. Percentile size: {threshold_size}")
        #Create equally sized populations
        #We can know the sizes of the populations from the start.
        subpops = []
        ntile_size = int(cellyields.shape[0]*threshold_size)-1
        #Theoretically: ntile_size = int(cellyields.shape[0]*(1-floor)*threshold_size)-1
        #Sort based on disc_pos value
        cellyields = cellyields[cellyields[:,disc_pos+base].argsort()]
        #Theoretically: cellyields = cellyields[:cellyields[floor*cellyields.shape[0]:-1,disc_pos+base].argsort()]
        sortedYields = cellyields
        for idx in range(0,int(1/threshold_size)):
            subpops.append(cellyields[0:ntile_size])
            if(cellyields.shape[0] > ntile_size):
                cellyields = np.delete(cellyields,np.s_[0,ntile_size],0)                
            #Grab percentile between idx to idx + threshold            
            subpops[idx] = cellyields[idx*ntile_size:(idx+1)*ntile_size]
            names.append(f"n:{len(subpops[idx])} Quantile: {idx*threshold_size:.2f}-{(idx+1)*threshold_size:.2f}")     
        subpops=np.asarray(subpops)
    else:
        raise NameError("Sorting method is not valid.")
    return subpops, names, sortedYields
     
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
        #Remove s&p noise
        Mask = np.where(Fo[i] > int(0.048*255),1,0)
        Mask = Mask.astype(np.uint8)        
        Mask = cv2.medianBlur(Mask,3)
        Mask = np.where(Mask>0,1,0)
        #ipdb.set_trace()
        Fv = np.subtract(Fm[i],Fo[i],dtype = np.int8)
        #Floor to zero
        Fv = np.multiply(np.clip(Fv,0,255),Mask)
        #Fv = np.clip(Fv,0,255)
        #ipdb.set_trace()
        #cYield = np.divide(Fv,Fm[i],out=np.zeros_like(Fv),where=Fm[i]!=0)
        cYield = np.divide(Fv.astype(np.float16),Fm[i].astype(np.float16),out=np.zeros_like(Fv, dtype=np.float16),where=Fm[i]!=0)
        cYield = np.round(cYield,decimals=3)
        Yield.append(cYield)
    return np.asarray(Yield)
    
def create_enchancement_mask(mask,FM,threshold):
    """
    Applies enchancement mask based on 10.1.1.6 of ImagingPAM manual. To be used in conjunction with
    cell_mask AOI mode. Only valid if dark-adapted before so that Fm truly gives Fm values. 
    Fv/Fm = 0 if Fm < 0.048 in that pixel.

    Parameters
    ----------
    imgstack : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    E_mask = np.where(FM>int(threshold*255),1,0)
    res = np.multiply(mask,E_mask)
    res[res>0]=1
    res = res.astype(np.uint8)
    return res
    

def create_Masks(imgstack, maskthres = 0.048):
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
    threshold = maskthres*255
    summed[summed < threshold] = 0
    
    #np.clip(summed,0,255,out=summed)
    summed[summed > 0] = 1
    summed = summed.astype(np.uint8)
    return summed    

def create_Masks_Ft(imgstack,maskthres=0.048):
    """
    Creates masks. Based on Ft values

    Parameters
    ----------
    imgstack : numpy array
        should be a MxNxP stack of yield images.
    maskthres : float, optional
        DESCRIPTION. The default is 0.048.

    Returns
    -------
    Mask.

    """
    Fo = imgstack[::2]
    th = []
    threshold = maskthres*256
    for i in range(len(Fo)):
        src = np.array(Fo[i])
        ret,img=cv2.threshold(src,int(threshold),1, cv2.THRESH_BINARY)
        th.append(cv2.medianBlur(img,3))
    mask = np.sum(th,axis=0)
    mask[mask>0]=1
    mask = mask.astype(np.uint8)
    return mask
    

def cleanup():
    import cv2
    import matplotlib.pyplot as plt
    cv2.destroyAllWindows()
    plt.close('all')

global outputdata
parser = argparse.ArgumentParser(description ="PAMalysis: A Python analysis script for analysing Microscopy-IPAM tif images. Made by Olle Pontén, Uppsala University 2021. GPL license v3.")
parser.add_argument('ProjectName',type=str, help="Project/output name")
parser.add_argument('-PAMSet','--PS',dest='PAMSet',type=str,default="PAMset.txt", help="Path/Name of PAMset file")
parser.add_argument('-FilePath','--FP', dest='proj_fp',help="Name of data file or data folder(if batch mode enabled). If batch mode is on and this argument is empty PAMalysis will analyse current folder.")
parser.add_argument('-b','--batch', dest='batch_flag', action='store_true', help="Enable batch mode")
parser.add_argument('-Debug',action='store_true', help="Switches on verbose output and debug code.")

args = parser.parse_args()
cur_dir = os.getcwd()
if(args.proj_fp is None):
    if(args.batch_flag):
        outputdata = perform_Analysis(cur_dir,args.ProjectName,cur_dir,args.batch_flag,args.PAMSet)
    else:
        print("Specify file path to analyze if batch mode not enabled. Exiting")
        sys.exit()
else:
    outputdata = perform_Analysis(args.proj_fp,args.ProjectName,cur_dir,args.batch_flag,args.PAMSet)



    