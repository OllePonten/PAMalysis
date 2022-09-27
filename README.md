# PAMalysis
PAMalysis: A Python analysis script for analysing Microscopy-IPAM tif images. Made by Olle Pontén, Uppsala University 2021. GPL license v3.
#Required packages
  Python >3.8
  numpy
  opencv >3.1
  Matplotlib
  tiffile
  

#Useage
PAMAlysis can either analyse a single stack of images or run in batch mode (-b) where it reads all tiff-files in the same folder and applies the same analysis.

#Run


#Output
  A CSV-formatted text file containing individual cells position and Fv/Fm value at every time frame. This file also lists with what settings the current analysis was done.
  An image stack containing all Fv/Fm images for manual analysis
  A binary cell-mask which shows the contours of all cells
  The same cell mask with cell indexes mapping all contours. The indexes correspond to the indexes in the CSV file


