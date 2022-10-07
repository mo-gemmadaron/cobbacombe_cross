#!/usr/bin/env python3.8
# -*- coding: iso-8859-1 -*-

'''
horizon_file_from_tree_elevations.py

Code for sorting/grouping tree elevations by azimuth to generate test horizon file

To run the code:
   ./horizon_file_from_tree_elevations.py

Author: gdaron
'''

import csv
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


def import_tree_elevation(file_path):

    filename = '{0}/Cobbacombe_tree_evelations.csv'.format(file_path) 

    if not os.path.isfile(filename):
        return -1

    tree_elev_file = pd.read_csv(filename)

    return(tree_elev_file)


def main():

    #1. Set input and output paths

    file_path = '/data/users/gdaron/cobbacombe/tree_heights'

    out_path = '/data/users/gdaron/cobbacombe/tree_heights'

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    #2. Get tree elevation data 

    data_df = import_tree_elevation(file_path)
    print(data_df)

    #3. Round azimuth to nearest 1 degree (and force 360 data to 0)

    data_df.loc[data_df['azimuth'] >= 359.5, 'azimuth'] = 0
    	
    data_df['azimuth_round'] = data_df['azimuth'].round(decimals=0)

    #4. Group data by azimuth and find max elevation angle in that sector

    scenarios_list = ['Elevation_angle', 'Elevation_angle_1m', 'Elevation_angle_2m', 'Elevation_angle_3m', 'Elevation_angle_4m']

    compile_df = pd.DataFrame(index = pd.RangeIndex(start=0.0, stop=360.0, step=1.0))

    for scenario in scenarios_list:

    	horizon_file = data_df.groupby(['azimuth_round'])[scenario].max()

    	horizon_file[horizon_file < 0 ] = 0

    	horizon_file.to_csv(os.path.join(out_path, 'horizon_file_{0}.csv'.format(scenario)))

    	compile_df = pd.merge(compile_df, horizon_file, left_index = True, right_index = True, how = 'inner')

    compile_df.columns = scenarios_list

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    for col in compile_df.columns:
    	compile_df[col].plot(label=f'{col}', linewidth =1.0)
    ax1.set_xlabel('Azimuth / deg')
    ax1.set_xlim([0, 359])
    ax1.set_ylabel('Elevation / deg')
    plt.xticks(np.arange(0, 359, 10), rotation=90)
    plt.legend(ncol=2, loc='upper left', prop={'size':8}, title='Tree height scenario')
    ax1.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, 'horizon_file_from_tree elevations.png'))
    plt.show()


if __name__ == '__main__':
    main()


