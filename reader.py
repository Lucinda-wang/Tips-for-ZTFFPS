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
#np.trapz = trapezoid if you use the old version of numpy, use this 

"""
Now, we have to read the data from ZTFFPS file, which is a txt file that you recieved from ZTF.
Please insert the path to your ZTFFPS file in the following line:
"""

class SN_reader_ztf():
    """
    Parameters
    ----------
    dpath : str
        The path to the ZTF data file (e.g., self.dpath = f'/the/way/to/your/ztffps_file.txt')
    n : int
        The line number to read from the ZTF data file (0-indexed).
    """
    def __init__(self, dpath):
        self.dpath = dpath
    """
    1. The script used to read the line on by one. Firstly, we replace the space with comma.
    2. Find lines containing target_word (skip other # lines). 
    3. Skip all lines containing #, but keep target_line.
    """
    def read_line(self, n:int): 
        with open(self.dpath, 'r') as f:
            file = f.readlines()
            new_data = file[n].replace(" ", ",").split(",") 
            clean_data = [str(x) for x in new_data if x != '']
        return clean_data
    
    def read_index(self): 
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions = []
        for line_number, line in enumerate(lines, start = 0): 
            if target_word in line: 
                pass
            else:
                positions.append((line_number))
        return np.array(positions)
    
    def read_data_index(self): 
        with open(self.dpath, "r", encoding = "utf-8") as file:
            lines = file.readlines()
            target_word = '#'
            positions = []

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

    def earn_data(self, key=None):
        """
        Reads the ZTF data file and extracts relevant information based on specified criteria.
        Parameters
        ----------
        key: sequence of str, only access with the following keys (e.g. 'mjd', 'filter', 'flux', 'flux_unc', 'zpdiff').
        
        Returns
        -------
        data_dict: dict
        """
        with open(self.dpath, 'r') as f:
            lines = f.readlines()

        # Read the index
        index = SN_reader_ztf(self.dpath).read_index()
        data_index = SN_reader_ztf(self.dpath).read_data_index()
        
        # Catch the header line and clean it
        index_array = [item.strip(',#') for item in lines[index[0]].split()]

        # Define the data range (start_num + 1 to end_num - 1)
        start_num = data_index[0]
        end_num = data_index[1]
        target_lines = lines[start_num + 1 : end_num]

        # Capture the indices of required columns
        jd_index                 = index_array.index("jd")
        filter_index             = index_array.index("filter")
        zpdiff_index             = index_array.index("zpdiff")
        scisigpix_index          = index_array.index("scisigpix")
        infobitssci_index        = index_array.index("infobitssci")
        sciinpseeing_index       = index_array.index("sciinpseeing")
        forcediffimflux_index    = index_array.index("forcediffimflux")
        forcediffimfluxunc_index = index_array.index("forcediffimfluxunc")

        # Delete 'null' and trailing commas while reading
        parsed_rows = []
        for line in target_lines:
            if not line.strip():
                continue
            # Delete trailing commas after each field value
            row = [item.strip(',') for item in line.split()]
            if row[forcediffimflux_index] != 'null':
                parsed_rows.append(row)

        # 若無有效資料直接回傳空陣列
        if not parsed_rows:
            return np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

        # 6. 轉為 NumPy 矩陣，進行高效的向量化提取與轉型
        data_matrix = np.array(parsed_rows)

        jd_box                 = data_matrix[:, jd_index].astype(float)
        MJD_box                = jd_box - 2400000.5
        zpdiff_box             = data_matrix[:, zpdiff_index].astype(float)
        sci_box                = data_matrix[:, scisigpix_index].astype(float)
        info_box               = data_matrix[:, infobitssci_index].astype(float)
        seeing_box             = data_matrix[:, sciinpseeing_index].astype(float)
        forcediffimflux_box    = data_matrix[:, forcediffimflux_index].astype(float)
        forcediffimfluxunc_box = data_matrix[:, forcediffimfluxunc_index].astype(float)

        # Convert to lowercase and standardize the format (Ex. 'ZTF_g' -> 'ztf_g') 
        raw_filters            = data_matrix[:, filter_index]
        filter_box             = np.array(['ztf_' + f.split('_')[1] for f in raw_filters])
        # Calculate the SNT
        SNT_box                = forcediffimflux_box / forcediffimfluxunc_box

        # Boolean Masking
        sp = (info_box < 33554432) & (sci_box < 25) & (seeing_box < 4) & (SNT_box > 3)
        data_dict = {
            'mjd': MJD_box[sp],
            'filter': filter_box[sp],
            'flux': forcediffimflux_box[sp],
            'flux_unc': forcediffimfluxunc_box[sp],
            'zpdiff': zpdiff_box[sp]
        }
        if key:
            return data_dict[key.lower()]
        
        return data_dict