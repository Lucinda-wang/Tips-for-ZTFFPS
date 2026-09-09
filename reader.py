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
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
np.trapz = trapezoid 

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
        SNT = 3 #signal-to-noise ratio
        reader      = SN_reader_ztf(self.dpath)
        index       = reader.read_index() #列出所有＃的位置
        data_index  = reader.read_data_index() #只列出最後兩個＃的位置
        index_array = reader.read_line(index[0]) 
        end_num     = data_index[1]
        start_num   = data_index[0]

        col = {name: i for i, name in enumerate(index_array)}

        #Read the row data once, avoid repeat calling
        row_line = [reader.read_line(i) for i in range(start_num + 1, end_num)]
        sn_mjd, sn_filter, sn_zp, sn_flux, sn_fluxerr = [], [], [], [], []
        """
        Here is the criteria:
        sp = np.where((forcediffimflux_box != 'null') & (info_box < 33554432) & 
                                (sci_box < 25) & (seeing_box < 4) & (SNT_box > 3)) #此為最終篩選
        """
        for data in row_line:

            if data[col["forcediffimflux"]] == 'null' or data[col["forcediffimfluxunc"]] == 'null':
                continue
            f   = float(data[col["forcediffimflux"]])
            ferr= float(data[col["forcediffimfluxunc"]])
            
            if ferr <= 0 or (f/ferr) <= 3:
                continue
            
            info   = float(data[col["infobitssci"]])
            sci    = float(data[col["scisigpix"]])
            seeing = float(data[col["sciinpseeing"]])
            if info >= 33554432 and sci >= 25 and seeing >= 4:
                continue
            
            band = data[col["filter"]].split('_')[1]
            sn_filter.append(f'ztf_{band}')
            sn_mjd.append(float(data[col["jd"]]) - 2400000.5)
            sn_zp.append(float(data[col["zpdiff"]]))
            sn_flux.append(f)
            sn_fluxerr.append(ferr)

        data_dict = {
            'mjd': np.array(sn_mjd),
            'filter': np.array(sn_filter),
            'flux': np.array(sn_flux),
            'flux_unc': np.array(sn_fluxerr),
            'zp': np.array(sn_zp)
        }
        if key:
            return data_dict[key.lower()]
        return data_dict
    
    def earn_mag(self, key=None):
        SNT = 3 #signal-to-noise ratio
        lightcurve = SN_reader_ztf(self.dpath)
        flux    = lightcurve.earn_data()['flux']
        fluxerr = lightcurve.earn_data()['flux_unc']
        zp      = lightcurve.earn_data()['zp']
        mask  = ((flux / fluxerr) > SNT )
        mag   = zp[mask] - 2.5*np.log10(flux[mask])
        sigma = 1.0857*(fluxerr[mask]/flux[mask]) #1.0857 = 2.5/ln(10)
        data_dict = {
            'mag':mag,
            'sigma':sigma
        }
        if key:
            return data_dict[key.lower()]
        return data_dict
    
    def run_figure(self, sn_name:str):
        reader = SN_reader_ztf(self.dpath)
        time = reader.earn_data()['mjd']
        band = reader.earn_data()['filter']
        mag     = reader.earn_data()['mag']
        err     = reader.earn_data()['sigma']

        fmt_list = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'H', 'X']
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
        mpl.rcParams['font.family'] = "monospace" #"Times New Roman"
        mpl.rcParams['lines.linewidth'] = 10.0
        for i, b in enumerate(np.unique(band)):
            plt.errorbar(time[band == b], mag[band == b], label = f'{b}-band', yerr = err[band == b], fmt = fmt_list[i % len(fmt_list)], ecolor = 'grey', color = colors[i % len(colors)], elinewidth = 1, capsize = 3)
            plt.ylim(np.max(mag[band == b])+1, np.min(mag[band == b])-1)
            plt.legend()
            plt.ylabel('Apparent Magnitude')
            plt.xlabel('MJD')
            plt.title(f'griz Light curve of SN{sn_name}'.format(sn_name))
            plt.savefig(f'SN{sn_name}_ztffps.png', dpi=600, facecolor='w', edgecolor='w', orientation='portrait')
        return plt.show()

"""
Now, we have to read the data from ATLAS FPS file, which is a txt file that you recieved from ZTF.
Please insert the path to your ATLAS FPS file in the following line:
"""
class SN_reader_atlas():
    """
    Parameters
    ----------
    dpath : str
        The path to the ATLAS data file (e.g., self.dpath = f'/the/way/to/your/atlasfps_file.txt')
    n : int
        The line number to read from the ATLAS data file (0-indexed).
    """
    def __init__(self, dpath):
        self.dpath = dpath

    def read_line(self, n:int):
        f = open(self.dpath, 'r')
        file = f.readlines()
        if n == 0:
            data = file[:][n][1:]
            data2 = data.replace(" ", ",").split(",")
            data2 = [item for item in data2 if item.strip() != '']
        else:
            data = file[:][n][:]
            data2 = data.replace(" ", ",").split(",")
            data2 = [item for item in data2 if item.strip() != '']
        f.close()
        return data2
    
    def earn_data(self, key=None):
  
        index_value  = SN_reader_atlas(self.dpath).read_line(0)
        MJD_index    = index_value.index('##MJD') #set position for data
        uJy_index    = index_value.index('uJy')
        duJy_index   = index_value.index('duJy')
        ftype_index  = index_value.index('F')
        err_index    = index_value.index('err')
        x_index      = index_value.index('x')
        y_index      = index_value.index('y')
        maj_index    = index_value.index('maj')
        min_index    = index_value.index('min')
        apfit_index  = index_value.index('apfit')
        mag5sig_index= index_value.index('mag5sig')
        Sky_index    = index_value.index('Sky')

        # 1. Read a file once, deal with the whole row data
        with open(self.dpath, 'r') as f:
            lines = f.readlines()

        # Deal with the first row and split the following strip
        header = lines[0].lstrip('#').split()
        
        # 直接拿整張資料表做分割，不用每次呼叫 SN_reader_atlas
        data_rows = [line.split() for line in lines[1:] if line.strip()]

        # 2. 轉成 NumPy 二維陣列，進行一次性的「直行（Column）」提取與型態轉換
        data_matrix = np.array(data_rows)

        ftype_box   = data_matrix[:, ftype_index]
        err_box     = data_matrix[:, err_index].astype(float).astype(int)
        x_box       = data_matrix[:, x_index].astype(float)
        y_box       = data_matrix[:, y_index].astype(float)
        mag5sig_box = data_matrix[:, mag5sig_index].astype(float)
        MJD_box     = data_matrix[:, MJD_index].astype(float) #Firstly, make the float data
        maj_box     = data_matrix[:, maj_index].astype(float)
        min_box     = data_matrix[:, min_index].astype(float)
        uJy_box     = data_matrix[:, uJy_index].astype(float)
        duJy_box    = data_matrix[:, duJy_index].astype(float).astype(int)
        apfit_box   = data_matrix[:, apfit_index].astype(float)
        Sky_box     = data_matrix[:, Sky_index].astype(float)

        #new one
        mask  = ((duJy_box<10000)&(err_box==0)& 
                (x_box>100)&(x_box<10460)&
                (y_box>100)&(y_box<10460)&
                (maj_box<5)&(maj_box>1.6)&
                (min_box<5)&(min_box>1.6)&
                (apfit_box>-1)&(apfit_box<-0.1)&
                (mag5sig_box>17)&(Sky_box>17))

        date_array = np.array(MJD_box[mask])
        duJy_array  = np.array(duJy_box[mask])
        flux_array = np.array(uJy_box[mask])
        type_array = np.array(ftype_box[mask])
        unique_day = np.unique((date_array).astype(int)) #beginning with the average part
        avg_type_box = np.zeros(len(unique_day), dtype = str)
        avg_m_box = np.zeros(len(unique_day), dtype = float) #make average flux value for ezch day
        avg_sigma_box = np.zeros(len(unique_day), dtype = float)
        avg_err_box = np.zeros(len(unique_day), dtype = float) 
        avg_day_box = np.zeros(len(unique_day), dtype = float) 
        for i in range(len(unique_day)): 
            day_loc = np.where((date_array).astype(int) == unique_day[i]) #calculate how many times 
            n = len(day_loc[0])
            final_date = date_array[day_loc]
            final_type = type_array[day_loc]
            final_flux = flux_array[day_loc] #how many flux per day
            sigma = duJy_array[day_loc] #how many errors per day
            weight = 1/(sigma*sigma)
            #加權平均
            weight_avg = np.sum(weight*final_flux)/np.sum(weight)
            #加權標準差
            weight_std = 1/np.sqrt(np.sum(weight))
            #針對日期做平均
            weight_date = np.sum(final_date)/n

            avg_m_box[i] = weight_avg
            avg_err_box[i] = weight_std
            avg_day_box[i] = weight_date
            avg_sigma_box[i] = np.abs((2.5 / np.log(10)) * (weight_std / weight_avg)) #new method
            avg_type_box[i] = final_type[0] #it looks weird, but I don't give a fxxk

        sigma_3 = 23.9 - 2.5*np.log10(3*avg_sigma_box)
        o_loc = np.where((avg_type_box == 'o') & (avg_m_box > 3*sigma_3)) #make criteria for o-band data
        c_loc = np.where((avg_type_box == 'c') & (avg_m_box > 3*sigma_3))

        odate = avg_day_box[o_loc]
        cdate = avg_day_box[c_loc]
        omag  = 23.9 - 2.5*np.log10(avg_m_box[o_loc])
        cmag  = 23.9 - 2.5*np.log10(avg_m_box[c_loc])
        oflux = avg_m_box[o_loc]
        cfluc = avg_m_box[c_loc]
        oflux_err = avg_err_box[o_loc]
        cflux_err = avg_err_box[c_loc]
        omag_err  = avg_sigma_box[o_loc]
        cmag_err  = avg_sigma_box[c_loc]

        data_dict = {
            'odate': odate,
            'cdate':cdate,
            'omag':omag,
            'cmag':cmag,
            'omag_err':omag_err,
            'cmag_err':cmag_err
        }
        if key:
            return data_dict[key.lower()]

        return data_dict
    
    def run_figure(self, sn_name):
        reader= SN_reader_atlas(self.dpath)
        odate = reader.earn_data()['odate']
        cdate = reader.earn_data()['cdate']
        omag  = reader.earn_data()['omag']
        cmag  = reader.earn_data()['cmag']
        omag_err = reader.earn_data()['omag_err']
        cmag_err = reader.earn_data()['cmag_err']
        mpl.rcParams['font.family'] = "monospace" #"Times New Roman"
        mpl.rcParams['lines.linewidth'] = 10.0
        with plt.style.context(['science', 'high-vis']):
            plt.errorbar(odate, omag, label = 'o-band', yerr = omag_err, fmt = 'o', ecolor = 'grey', color = 'b', elinewidth = 2, capsize = 4)
            #The part is for c-band
            plt.errorbar(cdate, cmag, label = 'c-band', yerr = cmag_err, fmt = 'o', ecolor = 'grey', color = 'g', elinewidth = 2, capsize = 4)
            plt.ylim(np.max(omag)+2,np.min(omag)-2)
            plt.grid()
            plt.xlabel('MJD')
            plt.ylabel('Apparent Magnitude')
            plt.title('oc Light curve of '+ sn_name)
            plt.savefig(f'SN{sn_name}_atlasfps.png', dpi=600, facecolor='w', edgecolor='w', orientation='portrait')
            plt.legend(loc = 'upper right')
        return plt.show()