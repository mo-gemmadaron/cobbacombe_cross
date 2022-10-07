#!/usr/bin/env python3.8
# -*- coding: iso-8859-1 -*-

'''
horizon_file_creation.py

Code for generating horizon files for Cobbacombe Cross for different scenarios

To run the code:
   ./horizon_file_creation.py

Author: gdaron
'''

import csv
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


def import_orig_horizon(file_path):

    filename = '{0}/constant_nonstd_radar16_horizon_obstacles'.format(file_path) 

    if not os.path.isfile(filename):
        return -1

    orig_horizon = pd.read_csv(filename)

    return(orig_horizon)


def main():

    #1. Set input and output paths

    file_path = '/data/users/gdaron/cobbacombe/horizon_files'

    out_path = '/data/users/gdaron/cobbacombe/horizon_files'

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    #2. Get original horizon file 

    data_df = import_orig_horizon(file_path)

    #3. Edit and export horizon file for scenarios
    '''
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
    #SCENARIO 3 - WORST CASE prefer 2deg scan for continuous 275-68 sector
    data_df.loc[275:359, '16'] = 1.01
    data_df.loc[0:68, '16'] = 1.01

    out_path_scen3 = out_path + '/scen3'

    if not os.path.exists(out_path_scen3):
        os.makedirs(out_path_scen3)

    data_df.to_csv(os.path.join(out_path_scen3, 'constant_nonstd_radar16_horizon_obstacles'), index=False)


if __name__ == '__main__':
    main()


