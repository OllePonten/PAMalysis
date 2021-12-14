# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 13:50:34 2021

@author: ollep
"""

def read_all_files(fp):
    import os, tifffile
    imgs = []
    fl = os.listdir(fp)
    tifs = [fp+ "/" + file for file in fl if ".tif" in file]
    for fovidx, i in enumerate(tifs):        
        tif = tifffile.imread(i)
        imgs.append(tif)
        
    return imgs

def create_stitched_image(fp, borderx,bordery,timeframe):
    import os, tifffile, cv2, numpy as np

    imgs = []
    fl = os.listdir(fp)
    tifs = [fp+ "/" + file for file in fl if ".tif" in file]
    for fovidx, i in enumerate(tifs):        
        tif = tifffile.imread(i)
        imgs.append(tif[timeframe])
    i=0
    j=0
    vconcat=[]
    #cv2.imshow("Durr",np.asarray(imgs[0],np.uint8))
    seq=[imgs[0:4],imgs[4:8:]]
    patch_list=[]
    while j<2:
        patch_list.append([])
        while i<4:
            patch_list[j].append(seq[j][i][50:-50:,:])
            #cv2.imshow(f"patch img{j}", seq[j][i])
            i+=1
        i=0
        #print(len(patch_list[0]))
        vconcat.append(np.vstack(patch_list[j]))
        j+=1
    
    vconcat = np.asarray(vconcat,np.uint8)
    #cv2.imshow("V img", vconcat)
    hconcat = np.hstack((vconcat[:]))
    cv2.imshow("Large img", hconcat)
    return hconcat
            
def create_ICF(images, kernel, only_Fo=True):
    #Based on Pipeline for illumination correction of images for high-throughput microscopy
    #S. SINGH, M.-A. BRAY, T.R. JONES* & A. E . CARPENTER, doi: 10.1111/jmi.12178
    import cv2, numpy as np, scipy as sp
    import scipy.ndimage as nimg
    #Assume images are in form (x,n*2,640x480)
    Firsts = []
    for stack in images:
        Firsts.append(stack[1])
    mean_img = np.sum(Firsts,axis=0,dtype=np.uint16)
    smoothed_img = nimg.median_filter(mean_img,size=(kernel,kernel),mode="reflect")
    smoothed_img = smoothed_img.astype(np.float64)
    smoothed_img *= (255/np.max(smoothed_img))
    smoothed_img = smoothed_img.astype(np.uint8)
    cv2.imshow("Smooth",smoothed_img)

        