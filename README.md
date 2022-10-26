# PAMalysis
PAMalysis: A Python analysis script for analysing Microscopy-IPAM tif images. Made by Olle Pontén, Uppsala University 2021. GPL license v3.
#Required packages
  Python >3.8
  numpy
  opencv >3.1
  Matplotlib
  tifffile
  

# Useage
PAMAlysis can either analyse a single stack of images or run in batch mode (-b) where it reads all tiff-files in the same folder and applies the same analysis.

# How to run
PAMalysis is run as directly from the shell or via IPython or the Jupyter interactive shell. A typical run signature from cmd could look as follows:
> python PAMAlysis.py Test_Project --FP C:/My_Experimental_Data/ -b
This would run an analysis on all .tif files (-b = Batch Mode) contained in the "C:/My_Experimental_Data/" folder and output the results into a folder named "Test_Project"
The equivalent IPython call:
runfile('C:/PAMalysis/PAMAlysis.py', wdir='C:/', args ="Test_Project --FP My_Experimental_Data -b")
The call
> python PAMAlysis.py Test_Project --FP C:/My_Experimental_Data/ -PAMset PAMset.txt -b
would perform the same analysis but applying the analysis parameters found in PAMset.txt. For more information on how to structure the settings file look under the heading PAMset.

# Output
  A CSV-formatted text file containing individual cells position and Fv/Fm value at every time frame. This file also lists with what settings the current analysis was done.
  An image stack containing all Fv/Fm images for manual analysis
  A binary cell-mask which shows the contours of all cells
  The same cell mask with cell indexes mapping all contours. The indexes correspond to the indexes in the CSV file

# PAMset
PAMset is a settings file that direct various parameters during the execution of PAMalysis. All options are specified in the key=value form, which is the read as a dictionary by PAMalysis upon execution.
The following key-value pairs are settable.
- minsize : minimum size (in pixels) of objects
- maxsize : maximum size (in pixels) of objects
- start_point : start frame to analyze
- end_point : end frame to analyze
- AOI_Mode : Cell segmentation method. Ft_Masks, cell_mask or Projection are valid arguments. See article: for description on impact.
- Cell_mask_fp : If AOI_Mode=cell_mask this should be the path to the cell mask to use.
- global_coordinates : Whether to correlate metadata from tif-files to map cells position using a shared coordinate system.
- max_cells : Maximum amount of cells to count. If more cells are detected they are discarded.
- floor : All cells with an Fv/Fm below this are discarded.
- subpop_size : How big are the subpopulations. Specified as 0.25 e.g. for 4 equivalently larger subpopulation of [0:0.25:0.5:0.75:1.0]
- sorting_method : How to sort cells. Valid arguments are Static_Bins and Distribution. Static_Bins creates static bins as: [floor+subpopsize, floor+subpopsize*2,....]. Distribution instead distributes all cells into percentiles based on subpop_size.
- sorting_position : What frame to use for sorting
- border : Removes outer borders of this size to avoid vignetting effects.
- intervall : Time between each frame in minutes.
- SYD : Sudden Yield Drop. Valid arguments are 0 or 1. Removes cells that rapidly fall in Fv/Fm (0.3) in a single frame and stay at this low level. Used to remove cells that move from their current resting place.
- histogram : Whether to create histograms showing the distribution of cells Fv/Fm
- hist_start : From what frame to draw data from.
- hist_end : From what frame to draw comparison data from.
- plots : Whether to output plots. Valid arguments are 0 or 1.
- legends : Whether to include legends on the plots.
- errorbars : Whether to include errorbars on the plots.
