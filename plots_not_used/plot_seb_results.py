#!/usr/bin/env python3.8
# -*- coding: iso-8859-1 -*-

'''
plot_seb_results.py

Code for plotting seb results from summary files

To run the code:
   ./plot_seb_results.py

Author: gdaron
'''

import csv
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

'''
def import_summary_file(file_path):

    filename = '{0}/trees_removed_single_scen2_run_summary.txt'.format(file_path) 

    if not os.path.isfile(filename):
        return -1

    orig_horizon = pd.read_csv(filename)

    return(orig_horizon)
'''

def main():

    #1. Set input and output paths

    file_path = '/data/users/gdaron/cobbacombe/radar_seb/data'

    out_path = '/data/users/gdaron/cobbacombe/radar_seb/data/plots'

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    #2. Get summary file 

    filename = '{0}/trees_removed_single_scen2_run_summary.txt'.format(file_path) 

    #file = open(filename).readlines()

    #for line in file:
    	#if 'pixel method - pluvial:' in line:
    		#test = next(file)
    		#print(test)

    with open(filename, 'r') as infile:
    	for line in infile:
    		if 'pixel method - pluvial:' in line: 
    			T0_pix_pluvial = infile.readline()
    			print(T0_pix_pluvial, end="")
    			test = T0_pix_pluvial.partition(', ')[0]  
    			print(test)	
    			#print(next_line)
    			#value = next_line.partition('Properties =')
    			#print(value)
    		if 'pixel method - fluvial:' in line:
    			T0_pix_fluvial = infile.readline()
    			print(T0_pix_fluvial, end="")  	
    			T1_pix_fluvial = infile.readline()
    			print(T1_pix_fluvial, end="")  	
    			T2_pix_fluvial = infile.readline()
    			print(T2_pix_fluvial, end="")  	
    			T3_pix_fluvial = infile.readline()
    			print(T3_pix_fluvial, end="")  	
    			T4_pix_fluvial = infile.readline()
    			print(T4_pix_fluvial, end="")  	
    			T5_pix_fluvial = infile.readline()
    			print(T5_pix_fluvial, end="")  	
    			T6_pix_fluvial = infile.readline()
    			print(T6_pix_fluvial, end="")
    		if 'catch method - fluvial:' in line:
    			T0_catch_fluvial = infile.readline()
    			print(T0_catch_fluvial, end="")  	
    			T1_catch_fluvial = infile.readline()
    			print(T1_catch_fluvial, end="")  	
    			T2_catch_fluvial = infile.readline()
    			print(T2_catch_fluvial, end="")  	
    			T3_catch_fluvial = infile.readline()
    			print(T3_catch_fluvial, end="")  	
    			T4_catch_fluvial = infile.readline()
    			print(T4_catch_fluvial, end="")  	
    			T5_catch_fluvial = infile.readline()
    			print(T5_catch_fluvial, end="")  	
    			T6_catch_fluvial = infile.readline()
    			print(T6_catch_fluvial, end="")   	      
    
'''
    data_df = import_orig_horizon(file_path)

    #3. Edit and export horizon file for scenarios

    #SCENARIO 1 - prefer 1deg scan for 286-345 and 16-30 sectors
    data_df.loc[285:345, '16'] = 0.51
    data_df.loc[15:30, '16'] = 0.51

    out_path_scen1 = out_path + '/scen1'

    if not os.path.exists(out_path_scen1):
        os.makedirs(out_path_scen1)

    data_df.to_csv(os.path.join(out_path_scen1, 'constant_nonstd_radar16_horizon_obstacles'), index=False)

    #SCENARIO 2 - prefer 2deg scan for 286-345 and 16-30 sectors
    data_df.loc[285:345, '16'] = 1.01
    data_df.loc[15:30, '16'] = 1.01

    out_path_scen2 = out_path + '/scen2'

    if not os.path.exists(out_path_scen2):
        os.makedirs(out_path_scen2)

    data_df.to_csv(os.path.join(out_path_scen2, 'constant_nonstd_radar16_horizon_obstacles'), index=False)
'''

if __name__ == '__main__':
    main()


