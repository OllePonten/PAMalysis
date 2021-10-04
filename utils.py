# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 13:50:34 2021

@author: ollep
"""



def create_stitched_image(fp, borderx,bordery,timeframe):
    import os
    
    imgs = []
    fl = os.listdir(fp)
    tifs = [fp+ "/" + file for file in fl if ".tif" in file]
    for fovidx, i in enumerate(tifs):        
        tif = tifffile.imread(i)
        imgs.append(tif[timeframe])
    

i=0
j=0
vconcat=[]
seq=[imgs[0:4],imgs[4:8:],imgs[11:7:-1]]
patch_list=[]
while j<3:
    patch_list.append([])
    while i<4:
        patch_list[j].append(seq[j][i][0][50:-50:,:])
        i+=1
    i=0
    vconcat.append(np.vstack(patch_list[j]))
    j+=1


hconcat = np.hstack((vconcat[0], vconcat[1], vconcat[2]))
cv2.imshow("Large img", hconcat)