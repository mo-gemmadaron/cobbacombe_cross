#!/usr/bin/env python3.8
"""
country_specific_stats.py

Reads in the benefits data for each radar site and flood condition and extracts
the benefits values grouped by country.
The results are written to a text file (one per site/flood type) in a format
that can be used to generate histograms of the data.
"""
import argparse
import h5py
import iris
import numpy as np

#------------------------------------------------------------------------------
# Read in what and where parameters
#------------------------------------------------------------------------------
def get_parameters(args):
    """
    Uses command line information to generate dictionaries containing
    filenames and area parameters for use during the modelling process.
    The parametersare returned in the "info! dictioanry and files in the
    "files" dictionary.
    """
    args = parse_args(args)

    info = {'obs_dir': args.obs_dir,
            'ancil_dir': args.ancil_dir,
            'wind': args.wind,
            'countries': args.countries,
            'radar_sites': args.radar_sites,
            'out_prefix': args.out_prefix,
           }

    return info

#------------------------------------------------------------------------------
# Parse command line arguments
#------------------------------------------------------------------------------
def parse_args(args=None):
    """
    :param args: List of string equivalent to command line. A value of None
                 reads the command line.
    :type args: list

    :rtype: :py:class:`argparse.ArgumentParser`
    """
    parser = argparse.ArgumentParser(
        description=('Run the multiple wind direction benefit code'))

    parser.add_argument('--obs_dir', action="store",
                        default='/data/users/gdaron/cobbacombe/radar_seb/data',
                        help='Location on disk of directory containing benefits data'\
                             '(use $RAD_COV_DATA)')

    parser.add_argument('--ancil_dir', action="store",
                        default='/data/users/gdaron/cobbacombe/radar_seb/ancillary',
                        help='Location on disk of directory containing country data'\
                             '(use $RAD_COV_ANCIL)')

    # Defaulted values
    #-----------------
    parser.add_argument('--wind', action="store",
                        default='single',
                        help='Name of wind data observation type')

    parser.add_argument('--radar_sites', nargs='*', action="store",
                        #default=['StBeesHead', 'HighPark', 'SilothAirfield', 'Loughermore',
                        #'ClarkyHill', 'Aberman', 'Anglesey', 'Gwaenysgar'],
                        #default=['Bride', 'LeakinHill', 'Criffel'],
                        default=['scen1', 'scen2'],
                        help='Names of potential radar sites to use')

    parser.add_argument('--countries', nargs='*', action="store",
                        default=['England', 'Scotland', 'Wales', 'NI'],
                        help='Space-seperated list of coutries to use'\
                             '(England, Scotland, Wales, NI)')

    # File names
    #-----------
    parser.add_argument('--out_prefix', action="store",
                        default="",
                        help='Prefix for output file')

    args = parser.parse_args(args)

    return args


#---------------------------------------------------------------------------
def get_area_masks(country_file):
    """
    Code to read in the coutry file which contains country identifiers for
    each (land-based) pixek.
    The data for each area is stored as an array that can be laid over the
    benefits data array to identify the local benefit values.
    """
    # Read in country h5 file
    file_in = h5py.File(country_file, "r")
    print('successfully loaded country file!')

    # Get data array and missing data value
    uk_areas = file_in['/dataset1/data1/data'][()]
    missing = file_in['/what'].attrs['integer_missing_data_value']

    # Set identifiers for each country (based on default country map)
    country_masks = {'England': {'id': 2},
                     'Guernsey': {'id': 3},
                     'Eire': {'id': 4},
                     'IOM': {'id': 5},
                     'Jersey': {'id': 6},
                     'NI': {'id': 7},
                     'Scotland': {'id': 8},
                     'Wales': {'id': 9},
                     'missing': {'id': missing},}

    # Generate mask for each country in the list
    for country in country_masks:

        region = np.where(uk_areas == country_masks[country]['id'], 1, 0)
        country_masks[country]['mask'] = np.asarray(region, dtype=float)

    return country_masks

#---------------------------------------------------------------------------
#def get_current_benefits_data(info, site, hour, flood_type):
def get_current_benefits_data(info):
    """
    Read in the netcdf file for the current time and site and and extract
    the required benefits data as an iris cube and re-shape to match the
    orientation of the country masks.
    """

    # Read in the catchment fluvial data as a netcdf file
#    file_to_use = iris.sample_data_path(info['in_file'])
    file_to_use = info['in_file']
    cubes = iris.load(file_to_use, info['cube_name'])

    data_cube = cubes[0]

    # Extract benefits data for whole country and reshape
    ben_data = np.asarray(data_cube.data)
    ben_data = np.flipud(ben_data)

    return ben_data

#---------------------------------------------------------------------------
def generate_area_breakdown(ben_data, info, country_masks, area_data):
    """
    Calculates benefits for each area of the country
    """

    # Set cumulative total to zero
    total = 0

    # Loop over each region
    for country in info['countries']:#['England', 'Scotland', 'Wales', 'NI']:

        area_value = country_masks[country]['mask'] * ben_data
        area_sum = np.nansum(area_value)
        area_data[country][info['hour']] = area_sum
        total += area_sum

    area_data['Total'][info['hour']] = total


#---------------------------------------------------------------------------
def write_out_results(area_data, info):
    """
    Write the data to screen
    """

    print('==========================================================')
    print(info['in_file'])
    print(info['out_file'])

    # Timestamps for data header
    timestamps = ', '.join('     T{0:1}'.format(n) for n in range(info['hours']))

    with open(info['out_file'], 'w') as s_out:

        # Name of site tested and timestamps
        s_out.write("{name}, {site}:\n".format(**info))
        s_out.write("{0:8s}, {1}\n".format('Country', timestamps))

        # Benefits data for each area - constructed line-by-line
        for country in ['England', 'Scotland', 'Wales', 'NI', 'Total']:

            # Create data string for each hour
            data = ','.join('{0:8.0f}'.format(area_data[country][n])
                            for n in range(info['hours']))
            s_out.write('{0:8s},'.format(country) + data + '\n')


#---------------------------------------------------------------------------
# Main Program
def main(args=None):
    """ Main code """

    # Read in command line parameters
    info = get_parameters(args)

    # Set file names
    '''
    pluvial_only = {'name': 'pluvial_only', 'hours':1,
                    'cube_id': 'aux_pluvial_array',
                    'in_prefix': 'EWSN'}

    pluvial_all = {'name': 'pluvial_all', 'hours':1,
                   'cube_id': 'aux_pluvial_array',
                   'in_prefix': 'pluv_and_fluv_EWSN'}

    fluvial_catch = {'name': 'catchment_fluvial', 'hours':7,
                     'cube_id': 'catch_fluvial_array',
                     'in_prefix': 'EWSN'}

    fluvial_pixel = {'name': 'pixel_fluvial', 'hours':7,
                     'cube_id': 'pixel_fluvial_array',
                     'in_prefix': 'EWSN'}
    '''
    pluvial_only = {'name': 'pluvial_only', 'hours':1,
                    'cube_id': 'pluvial_array',
                    'in_prefix': 'trees_removed'}

    pluvial_all = {'name': 'pluvial_all', 'hours':1,
                   'cube_id': 'aux_pluvial_array',
                   'in_prefix': 'trees_removed'}

    fluvial_catch = {'name': 'catchment_fluvial', 'hours':7,
                     'cube_id': 'catch_fluvial_array',
                     'in_prefix': 'trees_removed'}

    fluvial_pixel = {'name': 'pixel_fluvial', 'hours':7,
                     'cube_id': 'pixel_fluvial_array',
                     'in_prefix': 'trees_removed'}

    country_file = "{ancil_dir}/countries_british_isles_1km.h5".format(**info)

    # (1) Get countries mask
    country_masks = get_area_masks(country_file)

    # (2) Get current site and hour data
    for site in info['radar_sites']:

        info['site'] = site
        #info['in_dir'] = '{obs_dir}/{site}_{wind}/'.format(**info)
        info['in_dir'] = '{obs_dir}/{site}/'.format(**info)

        # Loop over each type of flooding
        #for flood_type in [pluvial_only, fluvial_catch, fluvial_pixel, pluvial_all]:
        for flood_type in [fluvial_catch]:

#        for flood_type in [pluvial_all]:

            # Add flood type information to dictionary
            info.update(flood_type)

            # To enable Loughermore to be included, change cube id
            #if (site == 'Loughermore') and (flood_type == fluvial_catch):
                #info['cube_id'] = 'pixel_fluvial_array'

            area_data = {'England': {}, 'Wales': {}, 'Scotland': {},
                         'NI': {}, 'Total':{}}

            info['out_file'] = \
              "{in_dir}/{out_prefix}{name}_{site}_{wind}_benefits.txt".format(**info)

            # Loop over each hour
            for hour in range(flood_type['hours']):
                print(hour)

                info['hour'] = hour

                # Get current file name
                info['in_file'] = \
                   "{in_dir}/{in_prefix}_T{hour:1d}_{site}_run_data_pixel_fluvial_array.nc".\
                     format(**info)
                   #"{in_dir}/{site}_T{hour:1d}_{in_prefix}_{wind}_data.nc".\
                     #format(**info)

                # Construct current cube name
                #info['cube_name'] = "{site}_{cube_id}_T{hour}".format(**info)
                info['cube_name'] = "{in_prefix}_catch_fluvial_array_T{hour}".format(**info)

                # Extract benefits data for current hour
                ben_data = get_current_benefits_data(info)
'''
                # Get breakdown of benefits per region
                generate_area_breakdown(ben_data, info, country_masks,
                                        area_data)
            # (3) Write the data to the screen
            write_out_results(area_data, info)
'''

if __name__ == '__main__':
    main()
