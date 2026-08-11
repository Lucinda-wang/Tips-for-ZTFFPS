#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 6 2026
@author: lucindawang
"""

# Standard library imports
import re
import os
import astropy 
import openpyxl

# Basic package imports
import numpy as np 
import pandas as pd
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from matplotlib import lines
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
#np.trapz = trapezoid if you had the old version of numpy, use this 


#========================================================================
#========================================================================
class SN_reader_ztf():
    def __init__(self, sn_name):
        self.sn_name = sn_name
        self.sheet = pd.read_excel(self.path, engine = "openpyxl") 
        self.dpath = f'/the/way/to/your/file/{sn_name}_ztf_mag.txt'

    def read_line(self, num:int): # The script used to read the line on by one
        with open(self.dpath, 'r') as f:
            file = f.readlines()
            new_data = file[num].replace(" ", ",").split(",") # use comma to replace the space
            clean_data = [str(x) for x in new_data if x != '']
        return clean_data
    
    def read_index(self):
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions = []
        # Find lines containing target_word (skip other # lines)
        for line_number, line in enumerate(lines, start = 0):
            # Skip all lines containing #, but keep target_line
            if target_word in line:
                pass
            else:
                positions.append((line_number)) #, index
        return np.array(positions)
    
    def read_data_index(self):
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions = []
        # Find lines containing target_word (skip other # lines)
        for line_number, line in enumerate(lines, start = 0): #originally start = 0
            start = 0
            while True:
                index = line.find(target_word, start)
                if index == -1:
                    break
                positions.append((line_number, index+1))
                start = index + len(target_word)
        last_1 = int(positions[-1][0])
        last_2 = int(positions[-2][0])
        return last_2, last_1

class Lightcurve_ztf():
    def __init__(self, sn_name, sn_type):
        self.sn_name = sn_name
        self.sn_type = sn_type
        
    def earn_data(self):
        SNT_box = []
        MJD_box = []
        sci_box = []
        info_box= []
        seeing_box = []
        filter_box = []
        zpdiff_box = []
        forcediffimflux_box = []
        forcediffimfluxunc_box=[]
        index = SN_reader_ztf(self.sn_name, self.sn_type).read_index() #列出所有＃的位置
        data_index = SN_reader_ztf(self.sn_name, self.sn_type).read_data_index() #只列出最後兩個＃的位置
        index_array = SN_reader_ztf(self.sn_name, self.sn_type).read_line(index[0]) 
        end_num = data_index[1]
        start_num  = data_index[0]
        data_range  = int(end_num - start_num -1) #Remove the begining num and the last one

        """
        Here we comes to the main part, we will read the data from the txt file and skip the null value.
        """
        jd_index = index_array.index("jd")
        filter_index = index_array.index("filter")
        zpdiff_index = index_array.index("zpdiff")
        scisigpix_index = index_array.index("scisigpix") # Create new filtering criteria
        infobitssci_index  = index_array.index("infobitssci") #new filtering criteria
        sciinpseeing_index = index_array.index("sciinpseeing") #new filtering criteria
        forcediffimflux_index = index_array.index("forcediffimflux")
        forcediffimfluxunc_index = index_array.index("forcediffimfluxunc")
        
        for j in range(data_range): # Revise part
            data_array = SN_reader_ztf(self.sn_name, self.sn_type).read_line(start_num+1+j)
            if data_array[forcediffimflux_index] == 'null': # If flux is null, skip this line
                continue
            else:
                band_v = data_array[filter_index].split('_')
                filter_box.append('ztf_' + band_v[1])
                zpdiff_box.append(float(data_array[zpdiff_index]))
                sci_box.append(float(data_array[scisigpix_index])) #new
                info_box.append(float(data_array[infobitssci_index])) #new
                seeing_box.append(float(data_array[sciinpseeing_index])) #new
                MJD_box.append(float(data_array[jd_index]) - 2400000.5)
                forcediffimflux_box.append(float(data_array[forcediffimflux_index]))
                forcediffimfluxunc_box.append(float(data_array[forcediffimfluxunc_index]))
                SNT_box.append(float(data_array[forcediffimflux_index]) / float(data_array[forcediffimfluxunc_index]))
        SNT_box = np.array(SNT_box)
        MJD_box = np.array(MJD_box)
        sci_box = np.array(sci_box)
        info_box = np.array(info_box)
        seeing_box = np.array(seeing_box)
        filter_box = np.array(filter_box)
        zpdiff_box = np.array(zpdiff_box)
        forcediffimflux_box = np.array(forcediffimflux_box)
        forcediffimfluxunc_box = np.array(forcediffimfluxunc_box)
        sp = np.where((forcediffimflux_box != 'null') & (info_box < 33554432) & 
                      (sci_box < 25) & (seeing_box < 4) & (SNT_box > 3)) #此為最終篩選
        
        return  MJD_box[sp], filter_box[sp], forcediffimflux_box[sp], forcediffimfluxunc_box[sp], zpdiff_box[sp]