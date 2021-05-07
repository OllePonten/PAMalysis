# -*- coding: utf-8 -*-
"""
Created on Fri May  7 13:36:59 2021
Software for analysing tiff image stacks created by pacman
@author: Olle
"""

from dataclasses import dataclass
import numpy as np
import cv2
import math
import time
import IPAMRH
import pacmangui
import scipy
   
    
def make_Yield_Images(img_stack):
    """
    Parameters
    ----------
    img_stack : Numpy, shape = 480,640,2x
        Input image stack of only Fo/Fm.

    Returns
    -------
    Yield image stack.

    """
    #Assumes the initial 4 images have been removed
    Fo = img_stack[:,:,::2]
    Fm = img_stack[:,:,1::2]
    #Yield is defined as Fv/Fm or (Fm-Fo)/Fm
    Yield = np.divide((Fm-Fo),Fm)
    return Yield
    
def create_Masks(img):
    """
    

    Parameters
    ----------
    img : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
def extract_Grid():
    #Extracting grid
    #Return binary mask of interest
    return 0
    
def template_matching():
    #For using a template image of the microwell instead
    #Return binary mask of interest
    return 0