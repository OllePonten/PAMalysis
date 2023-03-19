# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Fri May 7 13:36:59 2021
Software for analysing tiff image stacks created by pacman
@author: Olle
"""
import sys, datetime, os, configparser
import tifffile
import csv
import warnings
import argparse
import pprint

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
import cv2

# CALL SIG:
#runfile('C:/Users/ollpo511/Documents/GitHub/PAMalysis/PAMAlysis.py', wdir='C:Users/ollpo511/Documents', args = '{Project_Name} {PAMSet filename}')


global filenames, DEBUG

DEBUG = False

if (DEBUG):
    import ipdb


def load_PAM_Params(fp="PAMset.txt"):
    print(f"Looking for PAMset at {fp}")
    try:
        with open(fp, mode='r') as file:
            lines = file.readlines()
            lines = [line for line in lines if(
                line[0] != '[' and line[-1] != ']')]
            lines = [line for line in lines if(line[0] != "#")]
            lines = [line.rstrip('\n') for line in lines]
            lines = [line.split('=') for line in lines]
            params = dict(lines)
            return params
    except FileNotFoundError:
        print("No pamset text file found, proceeding with defaults")
        return []


class PAMalysis:

    def __init__(self):
        self.settings = dict()
        self.hist_fig = None
        font = {'family': 'arial',
                'size': 18}
        plt.rc('font', **font)

    def alt_PAMset(file_path):
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file {file_path} was not found.")
        config = configparser.ConfigParser()
        config.read(file_path)
        ini_dict = {}
        for section in config.sections():
            section_dict = {}
            for key, value in config.items(section):
                if value.lower() == 'true':
                    section_dict[key] = True
                elif value.lower() == 'false':
                    section_dict[key] = False
                else:
                    section_dict[key] = value
            ini_dict[section] = section_dict
        return ini_dict

        

    def map_PAMSet(self, file_settings):
        """
        Maps keywords and their values to parameter names in PAMalysis in th PAMalysis.setttings dict.

        Parameters
        ----------
        file_settings : Dictionary
            Dictionary originating from PAMset file.

        Returns
        -------
        None.

        """
        #So that people can use whatever capitalization in PAMset
        file_settings = {k.lower(): v for k, v in file_settings.items()}
        keys = file_settings.keys()
        print(f"Found settings for: {keys}. Mapping.")
        if('minsize' in keys):
            try:
                self.settings['Minsize'] = int(file_settings['minsize'])
            except:
                print("Minsize badly formatted")
        if('maxsize' in keys):
            try:
                self.settings['Maxsize'] = int(file_settings['maxsize'])
            except:
                print("Maxsize badly formatted")
        if('floor' in keys):
            try:
                self.settings['floor'] = float(file_settings['floor'])
            except:
                print("floor badly formatted")
        if('subpop_size' in keys):
            try:
                self.settings['subpopthreshold_size'] = float(
                    file_settings['subpop_size'])
            except:
                print("subpopthreshold badly formatted")
        if('sorting_method' in keys):
            try:
                self.settings['sorting_meth'] = str(
                    file_settings['sorting_method'])
            except:
                print("Unknown sorting method")
        if('sorting_position' in keys):
            try:
                self.settings['sorting_pos'] = int(
                    file_settings['sorting_position'])
            except:
                print("Unknown sorting ")
        if('border' in keys):
            try:
                self.settings['Border'] = int(file_settings['border'])
            except:
                print("border badly formatted")
        if('threshold' in keys):
            try:
                self.settings['Threshold'] = float(file_settings['threshold'])
            except:
                print("Threshold badly formatted")
        if('intervall' in keys):
            try:
                self.settings['Intervall'] = int(file_settings['intervall'])
            except:
                print("Intervall badly formatted")
        if('position_intervalll' in keys):
            try:
                self.settings['Position_Intervalll'] = int(file_settings['position_intervalll'])
            except:
                print("Position Intervall badly formatted")
        if('syd' in keys):
            try:
                self.settings['SYD'] = bool(int(file_settings['syd']))
                print(self.settings['SYD'])
                if(self.settings['SYD'] == False):
                    try:
                        self.settings['filter_methods'].pop("SYD")
                    except:
                        pass
            except:
                print("SYD badly formatted")
        if('outline' in keys):
            try:
                self.settings['outline'] = int(file_settings['outline'])
            except:
                print("Outline thickness badly formatted")
                pass
        if('histogram' in keys):
            try:
                self.settings['create_Hists'] = bool(int(file_settings['histogram']))
            except:
                print("Histogram badly formatted")
                pass
        if('histogram_end' in keys):
            try:
                self.settings['Histogram_end'] = int(file_settings['histogram_end'])
            except:
                print("hist_end point badly formatted")
                pass
        if('histogram_start' in keys):
            try:
                self.settings['Histogram_start'] = int(file_settings['histogram_start'])
            except:
                print("Hist start point badly formatted")
        if('plots' in keys):
            try:
                self.settings['create_Plots'] = bool(
                    int(file_settings['plots']))
            except:
                print("Plots badly formatted")
                pass
        if('start_point' in keys):
            try:
                self.settings['start_point'] = int(
                    file_settings['start_point'])
            except:
                print("Startpoint badly formatted")
        if('end_point' in keys):
            try:
                self.settings['orig_end_point'] = int(
                    file_settings['end_point'])
            except:
                print("Endpoint badly formatted")
        if('AOI_Mode' in keys):
            try:
                self.settings['AOI_mode'] = str(file_settings['AOI_Mode'])
            except:
                print("AOI Mode unknown")
        if('legends' in keys):
            try:
                self.settings['legends'] = bool(int(file_settings['legends']))
            except:
                print("AOI Mode unknown")
        if('errorbars' in keys):
            try:
                self.settings['errorbars'] = bool(
                    int(file_settings['errorbars']))
            except:
                print("Bad format")
        if(self.settings['AOI_mode'] == "cell_mask"):
            try:
                self.settings['cell_mask_fp'] = self.job_folder + \
                    "/" + file_settings['Cell_mask_fp']
                print(f"Cell_mask_fp: {self.settings['cell_mask_fp']}")
            except:
                print("Cell map fp badly formatted")
        if('global_coordinates' in keys):
            try:
                self.settings['globalcoordinates'] = bool(
                    int(file_settings['global_coordinates']))
            except:
                print("Global coordinates badly formatted")
        if('max_cells' in keys):
            try:
                self.settings['cell_limit'] = int(file_settings['max_cells'])
            except:
                print("Max cells badly formatted")
        if('cell_distance' in keys):
            try:
                self.settings['cell_dist'] = int(
                    file_settings['cell_distance'])
            except:
                print("Cell distance badly formatted")
        if('font_size' in keys):
            try:
                self.settings['font_size'] = int(file_settings['font_size'])
                font = {'family': 'Arial',
                        'size': self.settings['font_size']}
                plt.rc('font', **font)
            except:
                print("Font size badly formatted")
        if("Show_Cell_Mask" in keys):
            try:
                self.settings['show_cell_mask'] = bool(
                    file_settings["Show_Cell_Mask"])
            except:
                print("Show cell mask badly formatted.")
        pprint.pprint(self.settings)

    def perform_analysis(self, fp, work_name, job_folder, batch=False, pamset=None):
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
        tifs = []
        #### DEFAULT SETTINGS ####
        self.settings['AOI_mode'] = "Ft_Masks"
        self.settings['Border'] = 50
        self.settings['Intervall'] = 5
        self.settings['Minsize'] = 5
        self.settings['Maxsize'] = 60
        self.settings['subpopthreshold_size'] = 1
        self.settings['Floor'] = 0.05
        self.settings['Threshold'] = 0.048
        self.settings['create_Hists'] = False
        self.settings['create_Plots'] = True
        self.settings['globalcoordinates'] = False
        self.settings['start_point'] = 0
        self.settings['legends'] = True
        self.settings['errorbars'] = True
        self.settings['orig_end_point'] = -1
        self.settings['end_point'] = -1
        self.settings['Histogram_start'] = 1
        self.settings['Histogram_end'] = -1
        self.settings['filter_methods'] = {"SYD": 0.2}
        self.settings['cell_mask_fp'] = ""
        self.settings['sorting_meth'] = "Static Bins"
        self.settings['sorting_pos'] = 1
        self.settings['PlotAvg'] = False
        self.settings['cell_limit'] = 2000
        self.settings['cell_dist'] = -1
        self.settings['font_size'] = 16
        self.settings['show_cell_mask'] = False
        self.settings['verbosity'] = 1
        self.settings['outline']=1
        #########################
        fp_PAMset = None
        if(batch):
            fp_PAMset = fp + "/" + pamset
        else:
            fp_PAMset = fp[:fp.rindex("/")]+"/"+pamset

        PAMsettings = load_PAM_Params(fp_PAMset)
        if(len(PAMsettings) > 0):
            self.map_PAMSet(PAMsettings)

        font = {'family': 'Arial',
                'size': self.settings['font_size']}
        plt.rc('font', **font)

        self.outYields = dict()
        self.figures = dict()
        self.project_name = work_name
        self.output_folder = fp+"/Output/"+work_name
        if(batch):
            fl = os.listdir(fp)
            tifs = [fp + "/" + file for file in fl if ".tif" in file]
            print("Running batch mode")
            print(f"Analysing following files: {str(tifs)}")
            if(len(tifs) == 0):
                print("No tif files found at provided filepath. Aborting")
                sys.exit()
            self.no_files = len(tifs)
        else:
            tifs = [fp]
            self.output_folder = "/".join(fp.split("/")
                                          [:-1])+"/Output/"+work_name
            self.no_files = 1
        # Clean up output folder
        try:
            os.makedirs(self.output_folder)
        except OSError:
            print("Could not create output folder")
        try:
            os.remove(f'{self.output_folder}/' + work_name+'AllYields.csv')
        except:
            pass
        print(f"Saving output to {self.output_folder}")
        # Start of actual analysis: Read files
        for fovidx, current_tif in enumerate(tifs):
            if('/' in current_tif):
                fn = current_tif.split('/')[-1][:-4]
            else:
                fn = current_tif
            try:
                tif = tifffile.imread(current_tif)
                print(f"Analysing {fn}.tif")
            except:
                input(
                    f"Could not load image: {fn} Please check filename. Aborting")
                sys.exit()
            print(
                f"Read in file: {current_tif} with {int((len(tif)-4)/2)} time points")
            if(len(tif) > 6):
                # Analyse everything
                if(self.settings['orig_end_point'] == -1):
                    self.settings['end_point'] = len(tif)
                    tif = tif[4+(self.settings['start_point']*2)
                                 :4+(self.settings['end_point'])]
                    if(self.settings['verbosity'] >= 2):
                        print(
                            f"Analysing slices: {5+self.settings['start_point']*2}:{self.settings['end_point']}")
                else:
                    self.settings['end_point'] = self.settings['orig_end_point']
                    tif = tif[4+(self.settings['start_point']*2)
                                 :4+(self.settings['end_point']*2)]
                    if(self.settings['verbosity'] >= 2):
                        print(
                            f"Analysing slices: {5+self.settings['start_point']*2}:{self.settings['end_point']*2}")
                # Assume first are Black/NIR/Start images.
            else:
                tif = tif[0:2]
            tif_tags = {}
            globalx = 0
            globaly = 0
            # Get all tags
            with tifffile.TiffFile(current_tif) as tif_file:
                for tag in tif_file.pages[0].tags.values():
                    name, value = tag.name, tag.value
                    tif_tags[name] = value
                if('ImageDescription' in tif_tags.keys() and self.settings['globalcoordinates']):
                    desc = tif_tags['ImageDescription']
                    try:
                        globalx = tif_tags["xposition"]
                    except:
                        xposind = desc.find("xposition")+len("xposition")+1
                        globalx = desc[xposind:xposind +
                                       desc[xposind:].find("\n")]
                        globalx = int(
                            ''.join([i for i in globalx if i.isdigit()]))
                    try:
                        globaly = tif_tags["yposition"]
                    except:
                        yposind = desc.find("yposition")+len("yposition")+1
                        try:
                            globaly = desc[yposind:yposind +
                                           desc[yposind:].find("\n")]
                            globaly = int(
                                ''.join([i for i in globaly if i.isdigit()]))
                        except:
                            globaly = int(desc[yposind:])
            self.imgwidth = 640
            self.imgheight = 480
            if(self.settings['Border'] > 0):
                yields = [frame[self.settings['Border']:self.imgheight-self.settings['Border'],
                                self.settings['Border']*2:self.imgwidth-self.settings['Border']*2]for frame in tif]
                yields = make_yield_images(yields)
            else:
                yields = make_yield_images(tif)
            if(DEBUG):
                cv2.imshow("First yield image", np.asarray(
                    yields[0], dtype=np.uint8))
            yields_for_img = (yields*255).astype(dtype=np.uint8)
            tifffile.imwrite(f'{self.output_folder}/' +
                             f"Yields_{work_name}_{fn}.tif", data=yields_for_img)
            mask = 0
            if("Projection" in self.settings['AOI_mode']):
                mask = create_Masks(yields, self.settings['Threshold'])
            elif("Ft_Masks" in self.settings['AOI_mode']):
                mask = create_Masks_Ft([frame[self.settings['Border']:self.imgheight-self.settings['Border'], self.settings['Border']
                                       * 2:self.imgwidth-self.settings['Border']*2]for frame in tif], self.settings['Threshold'])
                #yields = np.multiply(yields,mask)
            elif("Cell_mask" in self.settings['AOI_mode']):
                print(
                    f"Reading: {self.settings['cell_mask_fp']} as cell mask image")
                try:
                    cell_mask = tifffile.imread(self.settings['cell_mask_fp'],)
                except ValueError:
                    print("No/Incorrect binary cell mask image provided")
                    exit()
                mask = cell_mask
            mask = mask.astype(dtype=np.uint8)
            if(DEBUG):
                cv2.imshow("Mask", mask)
            cnts, hrs = cv2.findContours(
                mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            contimg = cv2.cvtColor(np.zeros_like(mask), cv2.COLOR_GRAY2BGR)
            
            drawn = cv2.drawContours(contimg, cnts, -1, (0, 0, 255), -1)
            if(DEBUG):
                cv2.imshow("Conts", drawn)
            print(f"Raw contours found: {len(cnts)}")
            filt_cnts = self.filter_conts(cnts)
            print(f"Contour list pruned to {len(filt_cnts)}")
            
            
            # Draw a filtered mask
            filteredMask = cv2.drawContours(cv2.cvtColor(np.zeros_like(
                mask, dtype=np.uint8), cv2.COLOR_GRAY2BGR), filt_cnts, -1, (255, 255, 255), -1)
            #Draw thicker outlines to accomodate "wiggling" cells, controlled by PAMset: outline=
            if(self.settings['outline']>1):
                filteredMask = cv2.drawContours(filteredMask, filt_cnts, -1,(255,255,255),self.settings['outline'])
            #filteredMask = cv2.drawContours(filteredMask, filt_cnts, -1, (255, 255, 255), 2)
            cell_mask_out = cv2.cvtColor(filteredMask, cv2.COLOR_BGR2GRAY)
            tifffile.imwrite(f'{self.output_folder}/' + work_name +
                             '_' + fn + "cell_mask.tif", cell_mask_out)
            if(DEBUG):
                cv2.imshow("Filtered conts", filteredMask)
                
                
            # Output table
            meanYields = np.zeros(shape=(len(filt_cnts), len(yields)+3))
            numberedMask = cv2.resize(
                filteredMask, (filteredMask.shape[1] * 2, filteredMask.shape[0] * 2), interpolation=cv2.INTER_CUBIC)
            for cellidx, cnt in enumerate(filt_cnts):
                # Get minimum bounding rect
                rect = cv2.boundingRect(cnt)
                if(self.settings['globalcoordinates']):
                    cellcenter = [self.settings['Border'] + globalx+int(
                        (rect[0]+rect[2])/2)*4, self.settings['Border']*2 + globaly+int((rect[1]+rect[3])/2)*4]
                else:
                    cellcenter = [self.settings['Border'] + int(
                        (rect[0]+rect[2])/2)*4, self.settings['Border']*2 + int((rect[1]+rect[3])/2)*4]
                cellMinis = yields[:, rect[1]:rect[1]+rect[3], rect[0]:rect[0]+rect[2]]
                mini_mask = cell_mask_out[rect[1]:rect[1]+rect[3], rect[0]:rect[0]+rect[2]]
                cellMinis=cellMinis[:]*(mini_mask//255)
                numberedMask = cv2.putText(numberedMask, str(cellidx), (rect[0]*2, (rect[1]+int(
                    rect[3]/2))*2), cv2.FONT_HERSHEY_PLAIN, 1, (0, 255, 0), thickness=1)
                zeroedidx = []
                cellSize = np.count_nonzero(cellMinis[0])
                for timeidx, img in enumerate(cellMinis):
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        try:
                            meanYield = np.sum(img)/cellSize
                            #cv2.imshow("cell mask",cmo)
                            #print(np.count_nonzero(img))
                            #print(cv2.contourArea(cnt))
                        except Warning:
                            if(cellidx not in zeroedidx):
                                if(DEBUG):
                                    cv2.imshow("Cell_Mini", np.asarray(
                                        img*255, np.uint8))
                                if(self.settings['verbosity'] >= 2):
                                    print(
                                        f"Nan/Zero yield enc. Timeindex: {timeidx+self.settings['start_point']}. Cellindex: {cellidx}. Setting zero")
                                zeroedidx.append(cellidx)
                            meanYield = 0
                        if(timeidx == 0):
                            # The first digit tells what field of view a cell is from
                            meanYields[cellidx, 0] = (fovidx*1000)+cellidx
                            meanYields[cellidx, 1] = cellcenter[0]
                            meanYields[cellidx, 2] = cellcenter[1]
                    meanYields[cellidx, timeidx+3] = meanYield
                    meanYields[cellidx, timeidx+3] = round(meanYields[cellidx, timeidx+3], 3)
            
            # THRESHOLD FILTER
            filteredYields = self.filter_yields(
                meanYields[:, :], self.settings['filter_methods'], self.settings['Floor'])
            print(f"Yields remaining after filter: {len(filteredYields)}")

            if(self.settings['show_cell_mask'] or DEBUG):
                cv2.imshow("Numbered cell mask", numberedMask)
            cv2.imwrite(f'{self.output_folder}/' + work_name +
                        '_' + fn + "numbered_masks.tif", numberedMask)

            if(len(filteredYields)==0):
                raise ValueError("No contours remaining after filters applied. Exiting")
                
            if(DEBUG):
                cv2.imshow("Numbered masks", numberedMask)
            subs, names, sortedYields = self.subdivide_yield(
                filteredYields, threshold_size=self.settings['subpopthreshold_size'], floor=self.settings['Floor'])
            
            with open(f'{self.output_folder}/' + work_name+'AllYields.csv', mode='a', newline="") as tot_file:
                tot_yield_writer = csv.writer(
                    tot_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
                times = ["Index", "XPosition", "YPosition"] + list(range(0, len(
                    sortedYields[0]-3)*(self.settings['Intervall']), self.settings['Intervall']))
                if(fovidx == 0):
                    tot_yield_writer.writerow(
                        list(self.settings.items()) + [str(datetime.date.today())])
                    tot_yield_writer.writerow(times)
                with open(f'{self.output_folder}/' + work_name + '_' + fn + '.csv', mode='w', newline="") as pos_file:
                    print(f"Writing CSV to: {pos_file.name}")
                    yield_writer = csv.writer(
                        pos_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    yield_writer.writerow(fn)
                    yield_writer.writerow(self.settings.items())
                    yield_writer.writerow(times)
                    for idx, part in enumerate(sortedYields):
                        yield_writer.writerow(part)
                        tot_yield_writer.writerow(part)
                        
            if(self.settings['create_Plots']):
                try:
                    self.plot_values(subs, names, work_name, fn, fovidx)
                except Exception as e:
                    print(e)
                    print("Could not plot values.")
                    pass
            # For outside use, return our filtered yields
            self.outYields[f"{fn}"] = sortedYields

        if(self.settings['create_Hists']):
            unsorted_yields_start = []
            unsorted_yields_end = []
            for x in list(self.outYields.values()):
                #unsorted_yields_start = np.concatenate(unsorted_yields_start, x[:,2+hist_start])
                x = np.asarray(x)
                #unsorted_yields_start = np.concatenate((x[:,2+hist_start]))
                unsorted_yields_start.extend(
                    x[:, 2+self.settings['Histogram_start']])
                if(self.settings['Histogram_end'] != -1):
                    #unsorted_yields_end = np.concatenate((unsorted_yields_end, x[:,2+self.settings['hist_end']-self.settings['start_point']]))
                    unsorted_yields_end.extend(
                        x[:, 2+self.settings['Histogram_end']-self.settings['start_point']])
            # print(len(list(outYields.values())))
            #unsorted_Yields = np.concatenate((list(outYields.values())),axis=1)
            self.hist_fig = self.plot_histograms(
                unsorted_yields_start, (self.settings['Histogram_start']-1)*self.settings['Intervall'], "blue", "//")
            if(self.settings['Histogram_end'] != -1):
                self.plot_histograms(unsorted_yields_end, (
                    self.settings['Histogram_end']-self.settings['start_point'])*self.settings['Intervall'], "green", "\\\\")

            self.hist_fig.legend(loc="upper right", bbox_to_anchor=(0.9, 0.88), prop={'size': 18})
            self.hist_fig.savefig(f"{self.output_folder}/{work_name}_Histogram", bbox_inches='tight')
        if(self.settings['create_Plots'] and self.settings['legends']):
            self.figures[f"{work_name} Average FvFm"].legend(
                loc="lower left", bbox_to_anchor=(0.12, 0.1), ncol=2, prop={'size': 18})
            self.figures[f"{work_name} Average FvFm"].savefig(
                f"{self.output_folder}/{work_name}_Average FvFm", )

        PlotAvg = False
        return self.outYields, self.figures, self.hist_fig

    def reanalyze(self, yields, indexes):
        # If you want to reanalyze only specific numbered cells.
        manFilteredYields = [part for part in yields if part[0] in indexes]
        return manFilteredYields

    def plot_values(self, yields, names, jobname, filename, subjob, rows=-1, columns=-1, mode="Lines", floor=0.2, legends=True, errorbars=True):
        # Assumes that yields is formatted as yields.shape = [n(subplots),n(samples),n(values)]
        global PlotAvg, sub_legends
        sub_legends = []
        rows = -1
        columns = -1
        color = None
        base = 3
        tot = len(names)
        if(tot == 1):
            columns = 1
            rows = 1
        else:
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
        xlim = [0, 0]
        xstep = 0
       
        
        xlim = [0, round((len(yields[0][0])-base) *
                         self.settings['Intervall'])]

        if(xlim[1] <= 100):
            xstep = 10
        else:
            xstep = round((xlim[1]/100)+1)*10
            if(xstep < 120 and xstep > 80):
                xstep = 100
        ylim = [0, 0.7]
        Position = range(1, tot + 1)
        plt.close(f"{jobname}: Subpopulations")
        fig = plt.figure(f"{filename} Subpopulations",
                         figsize=(12*columns, 10*rows))
        cmap = plt.cm.get_cmap("viridis_r")
        fig.suptitle("Subpopulations")
        for k, subyields in enumerate(yields):
            #norm = mcolors.Normalize(vmin=0,vmax=len(subyields))
            fig.suptitle("Subpopulations")
            # add every single subplot to the figure with a for loop
            ax = fig.add_subplot(rows, columns, Position[k])
            ax.set_title(names[k])
            ax.set_ylabel("$F_{V}$/$F_{m}$")
            ax.set_xlabel("Time [minutes]")
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)
            subyields = [trim_yield[3:] for trim_yield in subyields]
            norm = mcolors.Normalize(vmin=ylim[0]+0.2, vmax=0.55)
            avg_line = (np.mean(subyields, axis=0))
            avg_error = (np.std(subyields, axis=0))
            pop_size = len(subyields)
            avg_lines.append(avg_line)
            avg_errors.append(avg_error)
            avg_sizes.append(pop_size)
            for index, part in enumerate(subyields):
                ax.plot(range(xlim[0], xlim[1], self.settings['Intervall']), part,
                        marker='o', markersize=3, linewidth=0.5, color=cmap(norm(part[-1])))
            plt.yticks(np.arange(ylim[0], ylim[1], 0.1), labels=None)
            plt.xticks(np.arange(0, xlim[1], step=xstep))
            ax.tick_params(which='major', width=1.2,length=6)
            ax.tick_params(which='minor', width=0.75,length=2.5)
            ax.xaxis.set_major_locator(MultipleLocator(xstep))
            ax.xaxis.set_minor_locator(MultipleLocator(xstep/5))
            plt.minorticks_on()
            plt.grid(axis="y")

        fig.savefig(
            fname=f"{self.output_folder}/{jobname}_{subjob}_total_yields")

        cmap = plt.cm.get_cmap("viridis")
        norm = mcolors.Normalize(vmin=0, vmax=self.no_files-1)
        fig2 = plt.figure(f"{jobname} Average FvFm", figsize=[12, 10])
        ax = plt.subplot(111)
        plt.ylabel("$F_{V}$/$F_{m}$")
        plt.xlabel("Time [minutes]")
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.yticks(np.arange(ylim[0], ylim[1], step=(
            ylim[1]-ylim[0])/7), labels=None)
        plt.xticks(np.arange(0, xlim[1], step=xstep))
        plt.minorticks_on()
        ax.tick_params(which='major', width=1.2,length=6)
        ax.tick_params(which='minor', width=0.75,length=2.5)
        ax.xaxis.set_major_locator(MultipleLocator(xstep))
        ax.xaxis.set_minor_locator(MultipleLocator(xstep/5))
        plt.grid(True, axis="y")
        lw = 3
        msize = 4
        elw = 2
        opac = 0.75
        for idx, avgs in enumerate(avg_lines):
            if(legends and errorbars):
                ax.errorbar(range(xlim[0], xlim[1], self.settings['Intervall']), avgs, yerr=avg_errors[idx], color=cmap(norm(
                    subjob)), label=f"{filename}, " + f"n(cells)= {str().join([s for s in names[idx][3:12] if s.isdigit()])}", markersize=msize, marker='o', linewidth=lw, capsize=2, elinewidth=elw, errorevery=(1+subjob, self.no_files), alpha=opac)
            elif(errorbars):
                ax.errorbar(range(xlim[0], xlim[1], self.settings['Intervall']), avgs, yerr=avg_errors[idx], color=cmap(norm(
                    subjob)), markersize=msize, marker='o', linewidth=lw, capsize=2, elinewidth=elw, errorevery=(1+subjob, self.no_files), alpha=opac)
            elif(legends):
                ax.plot(range(xlim[0], xlim[1], self.settings['Intervall']), avgs, color=cmap(norm(
                    avgs[-1])), label=f"{filename}, " + f"n(cells)= {str().join([s for s in names[idx][3:12] if s.isdigit()])}", markersize=msize, marker='o', linewidth=lw, alpha=opac)
            else:
                ax.plot(range(xlim[0], xlim[1], self.settings['Intervall']), avgs, color=cmap(
                    norm(avgs[-1])), markersize=msize, marker='o', linewidth=lw, alpha=opac)

        self.figures[f"{jobname}_{subjob}_total_yields"] = fig
        self.figures[f"{jobname} Average FvFm"] = fig2

    def plot_histograms(self, yields, time_point=0, i_color="red", hatch_char="/"):
        print(f"Creating histograms for {self.project_name}")
        col_fig = plt.figure(f"{self.project_name}", figsize=[12, 10])
        plt.xlabel("$F_{V}$/$F_{m}$")
        plt.ylabel("Count [cells]")
        plt.minorticks_on()
        plt.grid(visible=True, which='major', axis="y")
        yield_bins = np.linspace(0, 0.7, num=(round((0.7)/0.05)+1))
        yields = np.asarray(yields)
        avg = np.mean(yields)
        arr = plt.hist(yields, bins=yield_bins, alpha=0.7, label=f"n: {len(yields)}. T: {time_point} mins. Mean: {avg:.3f}",
                       edgecolor=i_color, align="mid", fill=False, orientation="vertical", hatch=hatch_char)
        plt.xticks(yield_bins, rotation=-45)
        roof = round(max(arr[0])/100+1, 0)*100
        # Make sure roof stays the same to keep scaling.
        if(roof < plt.ylim()[1]):
            roof = plt.ylim()[1]
        plt.ylim(0, roof)

        plt.axvline(avg, linestyle='dashed', color=i_color)

        return col_fig

    def filter_conts(self, cnts):
        """
        Filters contours so that they are at least {distance} pixels away from each other  

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

        # Filter based on proximity between cells
        dist_filtered_cnts = []
        cps = []
        radii = []
        discardlist = []
        import ipdb
        if(self.settings['cell_dist'] > 0):
            for idx, cnt in enumerate(cnts):
                (x, y), radius = cv2.minEnclosingCircle(cnt)
                cps.append([int(x), int(y),idx])
            for idx in range(len(cps)):
                for jdx in range(len(cps)):
                    # Check for self-intersection
                    if(cps[jdx] != cps[idx]):
                        cpdist = np.sqrt(
                            (cps[idx][0]-cps[jdx][0])**2 + (cps[idx][1]-cps[jdx][1])**2)
                        if (cpdist <= self.settings['cell_dist']):
                            # Discard
                            discardlist.append(idx)
            print(f"Removing {len(discardlist)} contours due to proximity.")
            dist_filtered_cnts = [i for j, i in enumerate(cnts) if j not in frozenset(discardlist)]
        else:
            dist_filtered_cnts = cnts
        
        
        size_filtered_cnts=[]
        # Filter based on size and being outside of border
        cent_dist = []
        for cnt in dist_filtered_cnts:
            mu = cv2.moments(cnt)
            if(mu['m00'] > self.settings['Minsize'] and mu['m00'] < self.settings['Maxsize']):
                x, y, w, h = cv2.boundingRect(cnt)
                maxY = self.imgheight-self.settings['Border']
                maxX = self.imgwidth-self.settings['Border']*2
                minX = self.settings['Border']*2
                minY = self.settings['Border']
                if(x+self.settings['Border']*2 > minX and x+w+self.settings['Border']*2 < maxX and 
                   y+self.settings['Border'] > minY and y+h+self.settings['Border'] < maxY):
                    size_filtered_cnts.append(cnt)
                    cntX = (x+w/2)-((maxX-minX)/2)
                    cntY = (y+h/2)-((maxY-minY)/2)
                    cent_dist.append(np.linalg.norm([cntX, cntY]))

        print(
            f"Contours remaining after filtering with size thresholds: {self.settings['Minsize']}-{self.settings['Maxsize']} pixels are {len(size_filtered_cnts)}")
        # Filter if there are too many contours
        while(len(size_filtered_cnts) > self.settings['cell_limit']):
            rmdist = max(cent_dist)
            ind = cent_dist.index(max(cent_dist))
            cent_dist.pop(ind)
            size_filtered_cnts.pop(ind)
            if(self.settings['verbosity'] >= 2):
                print(
                    f"Removed contours with {rmdist} from center. {len(size_filtered_cnts)} contours remaining. Limit: {self.settings['cell_limit']}")

        
        return size_filtered_cnts

    def filter_yields(self, cellyields, meths, floor=0):
        SYD = False
        sydthres = 1
        # For uniqueness, use set
        remidxs = set()
        tail = 5
        outputmsg = ""
        if('SYD' in meths.keys()):
            # Sudden Yield Drop. Filter based on corresponding threshold
            print("Filtering using sudden yield drop")
            sydthres = meths['SYD']
            SYD = True
        for idx, cell in enumerate(cellyields):
            if(SYD):
                for time in range(3, len(cell)-1):
                    if(cell[time]-cell[time+1] > sydthres):
                        # Remove, likely fell away
                        # Tail function so we don't remove on single frame brightness
                        if(time+tail <= len(cell)-1):
                            if(cell[time+tail] == 0):
                                remidxs.add(idx)
                        else:
                            remidxs.add(idx)
            if(floor > 0):
                if(cell[self.settings['sorting_pos']] < floor):
                    remidxs.add(idx)
        if(len(list(remidxs)) > 0):
            outputmsg += f"Filtered: {len(remidxs)} based on an sudden decrease of more than: {sydthres} or being lower than {floor} at t={self.settings['sorting_pos']}. {list(remidxs)}"
            print(outputmsg)
        if(len(remidxs) >= len(cellyields)):
            print("No cells remainig. Exiting.")
            sys.exit()
        return np.delete(cellyields, list(remidxs), axis=0)

    def subdivide_yield(self, cellyields, floor=0, threshold_size=0.25):
        # Several modes for finding subpopulations
        # Static Value = Divides up in subpopulations based on the yield value in disc pos
        # using static bins of threshold_size
        # Distribution = Subpopulates based on distributions based on threshold size so that each
        # subpopulation contains sample_size*threshold_size samples
        sortedYields = []
        names = []
        # First 3 positions are Index/X/Y
        base = 2
        if(self.settings['sorting_meth'] == "Static Bins"):
            subs = int(1/threshold_size)
            # We can't know the sizes of the populations from the start
            subpops = []
            for idx in range(0, subs):
                temp = []
                for cell in cellyields:
                    if((idx*threshold_size)) <= cell[self.settings['sorting_pos']+base] and cell[self.settings['sorting_pos']+base] < (((idx+1)*threshold_size)):
                        temp.append(cell)
                        sortedYields.append(cell)
                # Else it is empty
                if(len(temp) > 0):
                    subpops.append(temp)
                    name = f"{idx}. n: {len(temp)}. Threshold: {(idx*threshold_size):.3f}-{((idx+1)*threshold_size):.3f}"
                    names.append(name)
        elif(self.settings['sorting_meth'] == "Distribution"):
            print(
                f"Sorting based on percentile. Sorting position: {self.settings['sorting_pos']}. Percentile size: {threshold_size}")
            # Create equally sized populations
            # We can know the sizes of the populations from the start.
            subpops = []
            ntile_size = int(cellyields.shape[0]*threshold_size)-1
            #Theoretically: ntile_size = int(cellyields.shape[0]*(1-floor)*threshold_size)-1
            #Sort based on self.settings['sorting_pos'] value
            cellyields = cellyields[cellyields[:,self.settings['sorting_pos']+base].argsort()]
            sortedYields = cellyields
            for idx in range(0, int(1/threshold_size)):
                subpops.append(cellyields[0:ntile_size])
                if(cellyields.shape[0] > ntile_size):
                    cellyields = np.delete(cellyields, np.s_[0, ntile_size], 0)
                # Grab percentile between idx to idx + threshold
                subpops[idx] = cellyields[idx*ntile_size:(idx+1)*ntile_size]
                names.append(
                    f"n:{len(subpops[idx])} Quantile: {idx*threshold_size:.2f}-{(idx+1)*threshold_size:.2f}")
            subpops = np.asarray(subpops)
        elif(self.settings['sorting_meth'] == "None"):
            print("No sorting requested")
            names.append(f"{self.project_name}")
            subpops = cellyields
            sortedYields = cellyields
        else:
            print(self.settings['sorting_meth'])
            raise NameError("Sorting method is not valid.")
        if(len(sortedYields)==0):
            print("Could not sort cell values. Returning unsorted.")
            return "Unsorted",["Unsorted"],cellyields
        return subpops, names, sortedYields


def make_yield_images(img_stack):
    """
    Parameters
    ----------
    img_stack : Numpy, shape = 2x, 480, 640
        Input image stack of only Fo/Fm.

    Returns
    -------
    Yield image stack: x,480,640. Range of values are 0.0 to 1.0

    """
    # Assumes the initial 4 images have been removed
    Fo = img_stack[::2]
    Fm = img_stack[1::2]
    # Yield is defined as Fv/Fm or (Fm-Fo)/Fm
    Yield = []
    # Remove s&p noise
    Mask = np.where(Fm[0] >= int(0.048*255), 1, 0)
    Mask = Mask.astype(np.uint8)
    #Mask = cv2.medianBlur(Mask,3)
    #Mask = np.where(Mask>0,1,0)
    for i in range(len(Fo)):   
        Fv = np.subtract(Fm[i], Fo[i], dtype=np.float64)
        # Floor to zero
        Fv = np.clip(Fv, 0, 255)
        Fv = np.multiply(Fv, Mask)
        cYield = np.divide(Fv.astype(np.float32), Fm[i].astype(np.float32), out=np.zeros_like(
            Fv, dtype=np.float32), where=(np.logical_and(Fm[i] >= 3,Fo[i]>=1) ))

        Yield.append(cYield)
    return np.asarray(Yield)


def create_enchancement_mask(mask, FM, threshold):
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
    E_mask = np.where(FM > int(threshold*255), 1, 0)
    res = np.multiply(mask, E_mask)
    res[res > 0] = 1
    res = res.astype(np.uint8)
    return res


def create_Masks(imgstack, maskthres=0.048):
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
    summed = np.sum(imgstack, axis=0)
    # Only keep pixels as part of a cell if they have an average of threshold pixel
    # intensity over the entire stack, i.e. the length of the stack.
    #summed = (summed*255)/imgstack.shape[0]
    summed = (summed)/imgstack.shape[0]
    threshold = maskthres*255
    summed[summed < threshold] = 0

    # np.clip(summed,0,255,out=summed)
    summed[summed > 0] = 1
    summed = summed.astype(np.uint8)
    return summed


def create_Masks_Ft(imgstack, maskthres=0.048):
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
    threshold = maskthres*255
    for i in range(len(Fo)):
        src = np.array(Fo[i])
        ret, img = cv2.threshold(src, int(threshold), 1, cv2.THRESH_BINARY)
        th.append(cv2.medianBlur(img, 3))
        #th.append(img)
    mask = np.sum(th, axis=0)
    mask[mask > 0] = 1
    mask = mask.astype(np.uint8)
    return mask


def cleanup():
    import cv2
    import matplotlib.pyplot as plt
    cv2.destroyAllWindows()
    plt.close('all')


global outputdata
parser = argparse.ArgumentParser(
    description="PAMalysis: A Python analysis script for analysing Microscopy-IPAM tif images. Made by Olle Pontén, Uppsala University 2021. GPL license v3.")
parser.add_argument('ProjectName', type=str, help="Project/output name")
parser.add_argument('-PAMSet', '--P', dest='PAMSet', type=str,
                    default="PAMset.txt", help="Path/Name of PAMset file")
parser.add_argument('-Filepath', '--FP', dest='proj_fp',
                    help="Name of data file or data folder(if batch mode enabled). If batch mode is on and this argument is empty PAMalysis will analyse current folder.", default=None)
parser.add_argument('-batch', '--b', dest='batch_flag',
                    action='store_true', help="Enable batch mode")
parser.add_argument('-Debug', action='store_true',
                    help="Switches on verbose output and debug code.")
parser.add_argument('-Interactive', '--I', action='store_true',
                    dest='Interactive', help="Goes into interactive mode")

args = parser.parse_args()
cur_dir = os.getcwd()
PAMAlysis_Inst = PAMalysis()

# Interactive, gradual version
if(args.Debug):
    DEBUG = True
if(args.Interactive):
    print("Thank for you for using the interactive mode. Please note that this is still a beta feature, so I appreciate all feedback when using this analysis mode.\n")
    chosenPAMSet = False
    chosenFP = False
    chosenPN = False
    chosenBatch = False
    batch = False
    while(not chosenBatch):
        resp = input(
            "Do you want to analyze more than one IPAM tif file simultaneously?[Y/N]\n")
        if(resp == "Y"):
            batch = True
            chosenBatch = True
        elif(resp == "N"):
            batch = False
            chosenBatch = True
        else:
            resp = ""
            print("Aborted.")
    while(not chosenFP):
        FP = input(
            "Please provide a complete filepath to your singular datafile or to the folder where your data is stored.\n")
        confirm = input(f"Please confirm: {FP} [Y/N]\n")
        if(confirm == "Y"):
            chosenFP = True
        else:
            FP = ""
            confirm = ""
            print("Aborted.")
    while(not chosenPAMSet):
        resp = input(
            "Is your PAMSet file and your image data in the same folder? [Y/N]\n")
        if(resp == "Y"):
            PAMset = input(
                "Please input the file name for your PAMset file.\n")
            if(batch):
                PAMsetFP = FP+"/"+PAMset
            else:
                PAMsetFP = "/".join(FP.split("/")[:-1])+"/"+PAMset
        else:
            PAMsetFP = input("Please provide a complete filepath.\n")
        confirm = input(f"Please confirm: {PAMsetFP} [Y/N]\n")
        if(confirm == "Y"):
            chosenPAMSet = True
        else:
            FP = ""
            confirm = ""
            print("Aborted.")

    outputdata, plots, hists = PAMAlysis_Inst.perform_analysis(
        FP, args.ProjectName, cur_dir, batch, PAMsetFP)

elif(args.proj_fp is None):
    if(args.batch_flag):
        outputdata, plots, hists = PAMAlysis_Inst.perform_analysis(
            cur_dir, args.ProjectName, cur_dir, args.batch_flag, args.PAMSet)
    else:
        print("Specify file path to analyze if batch mode not enabled. Exiting")
        sys.exit()
else:
    outputdata, plots, hists = PAMAlysis_Inst.perform_analysis(
        args.proj_fp, args.ProjectName, args.proj_fp, args.batch_flag, args.PAMSet)
