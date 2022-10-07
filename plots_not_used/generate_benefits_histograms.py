#!/usr/bin/env python2.7
"""
generate_benefits_histograms.py

Code to generate histogram to demonstrate pluvial and fluvial benefits for
each potential radar site based on QPE and QPF data for each scenario.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

#---------------------------------------------------------------------------
# Set access from the command line
#---------------------------------------------------------------------------
def parse_args(args=None):
    """
    :param args: List of string equivalent to command line. A value of None
                 reads the command line.
    :type args: list

    :rtype: :py:class:`argparse.ArgumentParser`
    """
    parser = argparse.ArgumentParser(
        description=('Run the catchment-averaged opportunity benefit code'))
 
    parser.add_argument('--model', action="store", default='catch',
        help='Choice of catchment or pixel-based model (catch or pixel)')

    parser.add_argument('--winds', action="store", default='single',
        help='Choice of multi or single wind direction')
    
    parser.add_argument('--obs_dir', action="store", 
          default="/data/users/gdaron/cobbacombe/radar_seb/data/scen1",
          help='Location on disk of directory containing data')

    parser.add_argument('--ymax', action="store", 
          default=850000,
          help='Max benefit on y-axis')

    parser.add_argument('--yticks', action="store", 
          default=100000,
          help='Tick spacing on y-axis')

    parser.add_argument('--countries', nargs='*', action="store",
        default=['England', 'Wales', 'Scotland','NI'],
        help='Space-seperated list of countries to use'\
             '(Options are England, Scotland, Wales, NI)')

    args = parser.parse_args(args)

    # Main information pulled from teh command line
    info = {'obs_dir': args.obs_dir,
            'mode': args.model,
            'winds': args.winds,
            'ymax': args.ymax,
            'yticks': args.yticks,
            'countries': args.countries}

    return info

#---------------------------------------------------------------------------
# Extract data from datafiles
#---------------------------------------------------------------------------
def extract_data():
    """
    Extracts data from formatted data files

    The data is stored in a dictionary called "benefits"
    The dictionary then has subdictionaties with identical structures:

     pixel:        Pixel-wise fluvial data for all times with weighting applied
     catchment:    Catchment-based fluvial data for all times with weighting applied
     pluvial_only: Benefits for properties ONLY at pluvial risk
     pluvial_all:  Benefits for ALL properties at pluvial risk (includes those
                   also at risk of fluvial)
     catch_fluv_T0: Unscaled T0 fluvial data using catchment based calculation
     pixel_fluv_T0: Unscaled T0 fluvial data using catchment based calculation
     max:          Maximum benefit values (uses QPE at T0 if higher than total
                   QPF for forecast period)
    
     """
    # Create dictionary to hold data at each site and define column names
    benefits = {}

    # Extract data from summary text file
    # Generate sub-directory benefits sictionary for the chosen flood type

    benefits['pixel'] = extract_benefits(pixel)
    benefits['catchment'] = extract_benefits(catch)

    benefits['pluvial_only'] = extract_benefits(pluvial_only)
    benefits['pluvial_all'] = extract_benefits(pluvial_all)

    benefits['catch_fluv_T0'] = extract_benefits(catch_fluv_T0)
    benefits['pixel_fluv_T0'] = extract_benefits(pixel_fluv_T0)

    benefits['max'] = extract_benefits(maximum)

    # Generate dictionaries for data for combinaton plots
    for item in ['QPE_and_pluv', 'fluv_and_pluv', 'additional',
                 'add_by_country', 'maxfluv_and_pluv', 'maxadditional',
                 'max_add_by_country']:

        benefits[item] = {}

    # Manipulate fluvial data to use QPE instead of QPE+QPF
    #print "Adjusting fluvial values to use a maximum"
    #use_maximum_fluvial(['Aberman'], benefits)

    return benefits

#---------------------------------------------------------------------------
# Update fluvial data to use QPE if specified.
#---------------------------------------------------------------------------
def use_maximum_fluvial(max_sites, benefits):
    """
    For sites requiring the use of the QPE data (rather than the QPE+QPF
    data) modify the benefits array to use the unweighted T0 value. (This
    is equal to the original QPE values).
    """
    for site in max_sites:
        benefits['max'][site][:, 0] *= 2.0
        benefits['max'][site][:, 1:7] = 0
        benefits['max'][site][:, 7] = benefits['max'][site][:, 0]
    return

#---------------------------------------------------------------------------
# Main process for extracting benefits array
#---------------------------------------------------------------------------
def extract_benefits(scenario):
    """
    Extracts all benefits data for each site and benefits scenario into a
    single dictionary. Each site is stored as a sub-dictionary of the scenario
    dictionary.
    It is assumed that the data is stored in the data file in csv format as:

    <scenatio>, <site>
    Country,    T0,    T1,    T2,    T3,    T4,     T5,    T6
    England,    ?,     ?,     ?,     ?,     ?,      ?,     ?  
    Scotland,   ?,     ?,     ?,     ?,     ?,      ?,     ?  
    Wales,      ?,     ?,     ?,     ?,     ?,      ?,     ?  
    NI,         ?,     ?,     ?,     ?,     ?,      ?,     ?  
    Total,      ?,     ?,     ?,     ?,     ?,      ?,     ?  
    """

    print "Extracting benefits for: {0}".format(scenario['name'])
    print "-----------------------------------------------------"

    # Create dictionary to hold data at each site for chosen model
    scenario_benefits = {}

    # Loop over each radar site and extract data
    for site in site_list:

        # Generate datafile name
        filename = "{0}{1}_{2}_{3}_benefits.txt".\
                   format(scenario['mode'], scenario['file_id'], site, winds)

        # Generate location on disk of data file
        data_loc = "{0}/{1}_{2}/{3}".format(data_dir, site, winds, filename)

        print "File: {0}".format(data_loc)

        # Extract data if file exists
        if os.path.isfile(data_loc):

            # Extract data as a csv array. This assumes that the file is
            # formatted correctly
            data = np.genfromtxt(data_loc, dtype=int, skip_header=2,
                              usecols=scenario['data_columns'], delimiter=',')

            # Read in the row headers (usually countries & total)
            keys = np.genfromtxt(data_loc, dtype=str, skip_header=2,
                                 usecols=0, delimiter=',')

            # Convert data extracted to 2D array
            data = np.array(data, ndmin=2)

            # If data is a list (i.e. T0 only) convert to column vector
            if (np.shape(data)[0] == 1) and (np.shape(data)[1] > 1):
                data = data.transpose()

            # Add extra columns ro contain totals
            data = np.append(data, np.zeros_like(data), axis=1)

            for row in range(len(keys)):
                data[row, scenario['hours']] = np.sum(data[row,
                                                       0:scenario['hours']])

            # Multiply catch_fluv_T0 or pixel_fluv_T0 by 2 to get QPE
            # (i.e. unweighted benefits at T0)
            if (scenario == catch_fluv_T0) or (scenario == pixel_fluv_T0):
                data *= 2.0

            # Save to array
            scenario_benefits[site] = data[:, 0:scenario['hours']+1]

        else:
            # Return zeros if no data available
            scenario_benefits[site] = np.zeros((scenario['hours']+1, 7))

    return scenario_benefits


#---------------------------------------------------------------------------
# Generate colours for histogram
#---------------------------------------------------------------------------
def get_discrete_colourscheme(colourscheme, flood):
    """ Selects discrete range of colours from specified coulourscheme """

    multi_colour = []
    cmap = plt.cm.get_cmap(colourscheme)

    if colourscheme == 'BrBG':
        for t in range(flood['hours']):
            multi_colour.append(cmap(1.0 - (float(t+1)/16.0)))
    else:
        for t in range(flood['hours']):
            multi_colour.append(cmap(1.0 - (float(t+1)/8.0)))

    return multi_colour

#---------------------------------------------------------------------------
# Set time-varying colours for each region
#---------------------------------------------------------------------------
def set_area_colours(multi_colour, hour):
    """
    Set colour ranges for each area
    """
    colours = []

    for site in site_list:
        # Set colour key
        if site in Cumbria:
            colours.append(multi_colour['A'][hour])
        elif site in SE_Wales:
            colours.append(multi_colour['B'][hour])
        elif site in NE_Wales:
            colours.append(multi_colour['C'][hour])
        elif site in NI:
            colours.append(multi_colour['D'][hour])
        else:
            colours.append(multi_colour['E'][hour])
    return colours

#---------------------------------------------------------------------------
# Set uniform time-varying colours
#---------------------------------------------------------------------------
def set_uniform_colours(multi_colour, hour):
    """
    Set colour intensity ranges for each hour
    """
    colours = []

    for site in site_list:

        colours.append(multi_colour['U'][hour])

    return colours
#---------------------------------------------------------------------------
# Set country-varying colours for each region
#---------------------------------------------------------------------------
def set_country_colours(region, multi_colour, hour):
    """
    Set colour ranges for each area
    """
    colours = []

    # Set colour key
    if region == 'England':
        colours.append(multi_colour['A'][hour])
    elif region == 'Scotland':
        colours.append(multi_colour['E'][hour])
    elif region == 'NI':
        colours.append(multi_colour['B'][hour])
    elif region == 'Wales':
        colours.append(multi_colour['C'][hour])

    return colours

#---------------------------------------------------------------------------
# Generates histogram from extracted datasets
#---------------------------------------------------------------------------
def create_histogram(axis, flood, benefits, region, params):
    """
    Extracts the 'value' data for each site at a given timestep and adds it
    to the bar-chart. The data for each successive hour is 'stacked' on top
    of the previous hour's data by offsetting the lower limit of the bar by
    the sum of the previous nours' data (for the site).
    """
    # Initialise histogram settings
    site_benefits = {}
    sectors = []
    labels = []
    countries = []

    baseline = params['baseline']
    # Generate discrete colourset
    multi_colour = {'A': get_discrete_colourscheme('Reds', flood),
                    'B': get_discrete_colourscheme('Greens', flood),
                    'C': get_discrete_colourscheme('Purples', flood),
                    'D': get_discrete_colourscheme('Greys', flood),
                    'E': get_discrete_colourscheme('Blues', flood),
                    'U': get_discrete_colourscheme('BrBG', flood),}

    # Generate the data-sets containing all site data for each hour

    for level in range(params['max_levels']):

        site_benefits[level] = []

        if params['type'] == 'pluv and fluv lt':
            colours = set_area_colours(multi_colour, level)

        elif params['type'] == 'fluv lt':
            colours = set_area_colours(multi_colour, level)

        elif params['type'] == 'fluv T0':
            colours = set_area_colours(multi_colour, level)

        elif params['type'] == 'pluv and fluv lt exec':
            colours = set_uniform_colours(multi_colour, level)

        else:
            colours = set_country_colours(region, multi_colour, level)

        # Loop over each site (for the current level)
        for site in site_list:

            temp_value = benefits[flood['name']][site][area[region], level]
            site_benefits[level].append(temp_value)
               
        # Increase baseline for next collection of site data by adding the
        # current collection of data to it
        baseline[level+1] = np.add(baseline[level], site_benefits[level])

        # Add current time-data to barchart
        sector = axis.bar(params['xpos'],          site_benefits[level],
                          width=params['width'],   bottom=baseline[level],
                          color=colours,           align='edge',
                          linewidth=0,             zorder=10,)

        # Update list for legend data
        sectors.append(sector)

        # Special case for country plot
        if level == 0:
            countries.append(sector[0])

        # Lead-time labels: Check special cases
        if params['max_levels'] == 1:
            labels.append('QPE')

        elif flood in [pluvial_only, pluvial_all]:
            labels = [flood['name']]

        else:
            labels.append('T{0:1d}'.format(level))

   
    return sectors, labels, countries


#---------------------------------------------------------------------------
# Add yellow additional benefit box to plot
#---------------------------------------------------------------------------
def add_additional_benefit_data(additional_type, params, benefits, axis):
    """
    Adds the additional beenfits percentage on the bar chart as a yellow box
    """

    # Generate additional benefit value:
    #-----------------------------------
    # Get total benefit (combined pluv and fluv) as a list....
    if additional_type == 'additional':
        uplift = total_scenario(fluv_and_pluv, benefits)

    elif additional_type == 'maxadditional':
        uplift = total_scenario(maxfluv_and_pluv, benefits)

    # .... then calculate uplift for each site
    for i in range(len(uplift)):
            uplift[i] *= 0.282

    benefits[additional_type] = uplift

    # Add this additional data as first layer on the bar chart
    #----------------------------------------------------------
    sector = axis.bar(params['xpos'],        uplift,
                      width=params['width'], bottom=params['baseline'][0],
                      color='y',             linewidth=0,
                      align='edge',          zorder=10,)

    # Increase baseline for next collection of site data by adding the
    # current collection of data to it
    params['baseline'][0] = np.add(params['baseline'][0], uplift)
    

    # Update list for legend data
    all_sectors = [sector[0]]
    all_labels = ['Additional']

    return all_sectors, all_labels

#---------------------------------------------------------------------------
# Generate lead-time plots
#---------------------------------------------------------------------------
def plot_lead_time_benefits(scenario, benefits, title, outname, param_type,
                            components):
    """
    Creates a histogram for each site with the total benefit for each hour a
    different colour.
    """
    print "Generating benefits by lead-time for {name}".format(**scenario)

    # Generate axes
#    fig = plt.figure(figsize=(11.5, 8.0), facecolor="w")
    fig = plt.figure(figsize=(8.0, 6.5), facecolor="w")
    axis = fig.add_subplot(111, axisbg='none')

    # Extract histogram datasets and add to axes
    params = {'width': 0.65,
              'max_levels': scenario['hours'],
              'xpos': np.arange(len(site_list)),
              'baseline': {0: [0]*len(site_list)},
              'type': param_type,
              }

    # Initialise histogram settings
    all_layers = []
    all_labels = []
    baseline = params['baseline'][0]

    # (1) Plot additional benefits if needed and add labels to legend list
    #=====================================================================
    # This is a single layer of data plotted above the existing baseline
    if 'add' in components:
        all_layers, all_labels = \
           add_additional_benefit_data('additional', params, benefits, axis)

    if 'maxadd' in components:
        all_layers, all_labels = \
        add_additional_benefit_data('maxadditional', params, benefits, axis)

    # (2) Plot pluvial-only data if required
    #========================================
    # This is a single layer of data plotted above the existing baseline
    if 'pluv' in components:

        # Get pluvial benefit for the sites in a list
        pluvial_only_ben = total_scenario(pluvial_only, benefits)

        # Add to the barchart as the next layer of data
        sector = axis.bar(params['xpos'],    pluvial_only_ben,
                      width=params['width'], bottom=params['baseline'][0],
                      color='gray',          linewidth=0,
                      align='edge',          zorder=10,)

        # Increase baseline for next layer of data
        params['baseline'][0] = np.add(params['baseline'][0], pluvial_only_ben)

        # Update list for legend data
        all_layers.append(sector[0])
        all_labels.append('Pluvial')

    # (3) Plot fluvial (lead-time) data or other specialised flood type
    #==================================================================
    # This is up to 7 layers of data plotted above the existing baseline
    layers, labels, countries = create_histogram(axis, scenario, benefits,
                                                 'Total', params)

    # Apply formatting to current histogram
    format_histogram(all_layers+layers, all_labels+labels[0:7], title)

    # Save ro file
    plt.savefig("{0}/histograms/{1}".format(data_dir, outname))

    return


#---------------------------------------------------------------------------
# Generate total model benefit
#---------------------------------------------------------------------------
def total_scenario(flood, benefits):
    """
    Returns a list of the TOTAL scenario benefit oredered by site.
    This is used either to get the grey background in the country data OR to
    get the pluvial benefit data for each site
    """
    # Add total as a grey background to each site
    total_scenario_data = []

    # Get total benefits for each site
    for site in site_list:
        temp_value = benefits[flood['name']][site][4, flood['hours']]
        total_scenario_data.append(temp_value)

    return total_scenario_data

#---------------------------------------------------------------------------
# Generate lead-time plots by country
#---------------------------------------------------------------------------
def plot_lead_time_by_country(scenario, benefits, title, outname, hours):
    """
    Plot individual country benefits broken down by lead-time in front of a
    grey box showing the total value
    """
    country_labels = ['England', 'Scotland', 'Wales', 'NI', 'Total'] 
    country_key = []

    params = {'width': 0.2,
              'max_levels': hours,
              'baseline': {0: [0]*len(site_list)},
              'type': 'by country',
             }
    print "Generating benefits by lead-time and country"

    # Generate canvas
#    fig = plt.figure(figsize=(11.5, 8.0), facecolor="w")
    fig = plt.figure(figsize=(8.0, 6.5), facecolor="w")
    axis = fig.add_subplot(111, axisbg='none')

    # (1) Plot each country as a time-resolved bar chart
    #---------------------------------------------------
    for region in country_labels[0:4]:

        # Set x-position for current country
        params['xpos'] = \
                  np.arange(len(site_list))+(area[region]*(params['width']))

        # Extract histogram datasets and add to axes
        layers, labels, countries = create_histogram(axis, scenario, benefits,
                                                      region, params)

        # Update country colours list
        country_key.append(countries[0])

    # (2) Plot site total as a full-width bar behind the country-specific plot
    #-------------------------------------------------------------------------
    # Get total benefit for the current scenario
    background_total = total_scenario(scenario, benefits)

    # Plot grey boxes
    total = axis.bar(np.arange(len(site_list)),   background_total,
                      width=params['width']*4,    color='LightGray',
                      linewidth=0,                zorder=9,)

    # Update labels with "total" label
    country_key.append(total[0])

    # (3) Apply formatting to whole plot
    #-----------------------------------
    format_histogram(country_key, country_labels, title)

    # Save to file
    plt.savefig("{0}/histograms/{1}".format(data_dir, outname))
    print "{0}/histograms/{1}".format(data_dir, outname)

    return


#---------------------------------------------------------------------------
# Applies formatting to cyrrent histogram
#---------------------------------------------------------------------------
def format_histogram(sectors, labels, title):
    """
    Apply formating to current histogram. This includes:
    - Labelling x-axis with site names
    - Labelling y-axis with financial units
    - Setting appropriately sized legends
    - Adding plot titles
    """ 
    # (1) Title
    #-----------
    plt.title(title, fontsize=18)

    # (2) x-axis
    #-----------
    x_labels = ['HighPark', 'StBeesHead', 'SilothAirfield', 'Loughermore',
                'ClarkyHill', 'Coed Cae\nAberman', 'Anglesey', 'Gwaenysgor']

    plt.xticks(np.arange(len(x_labels)), x_labels, rotation=60,
               fontsize=14, zorder=50)

    # (3) y-axis
    #-----------
    ytick_locs = np.arange(0, ymax, yticks)
    ytick_names = ytick_locs/1000
    plt.yticks(ytick_locs, ytick_names, fontsize=14)    

    plt.ylabel('Value ('u"\xA3"'k p.a.)', fontsize=16)

    # (4) Legend (Set legend special cases)
    #---------------------------------------
    if len(labels) > 8:        # For Pluv + fluv + additional
        labels = labels[0:2] + ['Fluvial']
        sectors = sectors[0:3]
        plt.legend(sectors, labels, loc=2, ncol=1, fontsize=18)

    elif len(labels) > 7:      # For Pluv + fluv only
        labels = labels[0:1] + ['Fluvial']
        sectors = sectors[0:3]
        plt.legend(sectors, labels, loc=2, ncol=1, fontsize=18)

    elif len(labels) > 1:      # For fluvial lead time only
        plt.legend(sectors, labels, loc=2, ncol=1, fontsize=18)

    # (5) Set frame parameters
    #--------------------------
    plt.grid(color='grey', which='major', axis='y', linestyle='--')
    plt.xlim(xmin=-0.05)
    plt.tight_layout()

    return
#---------------------------------------------------------------------------
# Data dictionaries
#------------------
additional = {'name': 'additional','hours':1}

add_by_country = {'name': 'add_by_country','hours':7}

area = {'England': 0, 'Scotland': 1, 'Wales': 2, 'NI': 3, 'Total': 4}
col = {'hour': 0, 'props': 1, 'value': 2}  # Data array column headers

catch = {'name': 'catchment', 'grid':'catch', 'file_id':'fluvial', 'hours':7,
         'data_columns': range(1,8,1), 'scale': 1.0, 'mode': 'catchment_'}

catch_fluv_T0 = {'name': 'catch_fluv_T0', 'grid':'catch',
                 'file_id':'fluvial', 'hours':1, 'data_columns': range(1,8,1),
                 'scale': 2.0, 'mode': 'catchment_'}

fluv_and_pluv = {'name': 'fluv_and_pluv','hours':7}

pixel = {'name': 'pixel', 'grid':'pixel', 'file_id':'fluvial', 'hours':7,
         'data_columns': range(1,8,1), 'scale': 1.0, 'mode': 'pixel_'}

pixel_fluv_T0 = {'name': 'pixel_fluv_T0', 'grid':'catch',
                 'file_id':'fluvial', 'hours':1, 'data_columns': range(1,8,1),
                 'scale': 2.0, 'mode': 'pixel_'}

pluvial_only = {'name': 'pluvial_only', 'grid':'pixel',
                'file_id':'pluvial_only', 'hours':1, 'data_columns': 1,
                'scale': 1.0, 'mode': ''}

pluvial_all = {'name': 'pluvial_all', 'grid':'pixel',
               'file_id':'pluvial_all', 'hours':1, 'data_columns': 1,
               'scale': 1.0, 'mode': ''}

maximum = {'name': 'max', 'grid':'catch', 'file_id':'fluvial', 'hours':7,
           'data_columns': range(1,8,1), 'scale': 1.0, 'mode': 'catchment_'}

maxfluv_and_pluv = {'name': 'maxfluv_and_pluv','hours':7}
maxadditional = {'name': 'maxadditional','hours':1}
max_add_by_country = {'name': 'max_add_by_country','hours':7}
QPE_and_pluv = {'name': 'QPE_and_pluv','hours':1}

#site_list = ['HighPark', 'StBeesHead', 'SilothAirfield', 'Loughermore',
#             'ClarkyHill', 'Aberman', 'Anglesey', 'Gwaenysgar']
#site_list = ['HighPark', 'StBeesHead', 'SilothAirfield', 'Bride',
#             'LeakinHill', 'Criffel']
site_list = ['Scen1', 'Scen2']


#Cumbria = ['HighPark', 'StBeesHead', 'SilothAirfield', 'Criffel',
#           'LeakinHill']
#IOM = ['Bride']
#NI = ['Loughermore']
#Scotland = ['ClarkyHill']
#SE_Wales = ['Aberman', 'Merthyr']
#NE_Wales = ['Anglesey', 'MynyddYCwm', 'Gwaenysgar']


#############################################################################

# Extract data from the command line and default arguments
info = parse_args(args=None)

# Specify single or multi (averaged) winds
winds = info['winds']

# Select model type (catchment based or pixel-based)
if info['mode'] == 'catch':
    mode = catch
    mode_T0 = catch_fluv_T0

elif info['mode'] == 'pixel':
    mode = pixel
    mode_T0 = pixel_fluv_T0

else:
    print "Error! No model set (catch or pixel)"

# Specify main directory containing site sub-directories
data_dir = info['obs_dir']

# Set axis formatting
ymax = int(info['ymax'])
yticks = int(info['yticks'])

#---------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------
def main(args=None):
    """ Main code """

    # Extract benefits summaries from data files
    benefits = extract_data()

    # Construct bespoke benefit lists for each site of interest
    #===========================================================
    for site in site_list:

        # Generate QPE + pluvial only
        #----------------------------
        benefits['QPE_and_pluv'][site] = \
            (benefits['pluvial_only'][site] + benefits[mode_T0['name']][site])

        # Generate fluvial + pluvial only  (do for both standard and max scenario)
        #-------------------------------------------------------------------------
        # First add in fluvial data
        benefits['fluv_and_pluv'][site] = 1.0 * benefits[mode['name']][site]
        benefits['maxfluv_and_pluv'][site] =  1.0 * benefits['max'][site]

        # Then add pluvial to T0 data and then update total column
        for item in ['fluv_and_pluv', 'maxfluv_and_pluv']:

            benefits[item][site][:, 0] += benefits['pluvial_only'][site][:, 0]
            benefits[item][site][:,-1] += benefits['pluvial_only'][site][:, -1]

        # Add to additional benefits list (this is total value of fluv + pluv + add)
        benefits['add_by_country'][site] = benefits['fluv_and_pluv'][site]*1.282
        benefits['max_add_by_country'][site] = benefits['maxfluv_and_pluv'][site]*1.282


    # PLOTTING COMMANDS
    #======================

    # (1) Get fluvial time-specific totals
    #====================================
    plot_lead_time_benefits(mode, benefits,
         "Fluvial flooding ({0} based - {1}-wind)".format(mode['name'],winds),
         "Fluvial_flooding_by_lead_time_{0}.png".format(winds), 'fluv lt', [])

    # (2) Get pluvial-only time-specific (T0 only) totals
    #=====================================================
    plot_lead_time_benefits(pluvial_only, benefits,
         "Pluvial-only flooding (Pixel based calculation",
         "Pluvial_only_flooding_{0}.png".format(winds), 'fluv lt', [])

    # (3) Get pluvial-all time-specific (T0 only) totals
    #=====================================================
    plot_lead_time_benefits(pluvial_all, benefits,
         "All pluvial flooding (Pixel based calculation",
         "Pluvial_all_flooding_{0}.png".format(winds), 'fluv lt', [])

    # (4) Get country-specific breakdowns by lead-time
    #=================================================
    plot_lead_time_by_country(mode, benefits,
        "Fluvial flooding by lead-time and country "\
        "({0} based - {1}-wind)".format(mode['name'], winds),
        "Fluvial_flooding_by_time_and_country_{0}.png".format(winds), 7)

    # (5) Get pluvial plus fluvial
    #==================================
    plot_lead_time_benefits(mode, benefits,
         "Fluvial + pluvial-only",
         "Fig30a_Fluvial_plus_Pluvial_only_{0}.png".format(winds),
         'pluv and fluv lt', ['pluv'])

    # (6) Get fluvial T0 unscaled
    #====================================
    plot_lead_time_benefits(mode_T0, benefits,
         "Fluvial flooding ({0} based)".format(mode['name']),
         "Fluvial_flooding_T0_unscaled_{0}.png".format(winds), 'fluv lt', [])

    # (8) Get country-specific breakdowns by fluvial T0
    #=================================================
    plot_lead_time_by_country(mode_T0, benefits,
        "Fluvial flooding ({0} based) by country".format(mode['name']),
        "Fluvial_flooding_T0_by_country.png", 1)

    # (9) Get country-specific totals for pluvial-only
    #==================================================
    plot_lead_time_by_country(pluvial_only, benefits,
        "Pluvial-only flooding (Pixel based calculation",
        "Pluvial_only_flooding_by_country.png", 1)

    # (10) Get country-specific totals
    #==================================
    plot_lead_time_by_country(pluvial_all, benefits,
        "All pluvial flooding (Pixel based calculation",
        "Pluvial_all_flooding_by_country.png", 1)

    # (11) Get raw fluvial plus pluvial-only
    #========================================
    plot_lead_time_benefits(mode_T0, benefits,
         "QPE + pluvial-only",
         "Fig29a_QPE_and_Pluvial_only.png", 'pluv and fluv lt', ['pluv'])

    # (12) Get country-specific totals for QPE + pluv
    #===============================================
    plot_lead_time_by_country(QPE_and_pluv, benefits,
        "QPE + pluvial-only (by country)",
        "Fig29b_QPE_and_Pluvial_only_by_country.png", 1)

    # (13) Get country-specific totals for fluv + pluv
    #===============================================
    plot_lead_time_by_country(fluv_and_pluv, benefits,
        "Fluvial + pluvial-only (by country)",
        "Fig30b_Fluvial_and_Pluvial_only_by_country.png", 6)

    # (15) Get pluvial plus fluvial + additional
    #==========================================
    plot_lead_time_benefits(mode, benefits,
         "Fluvial + pluvial-only + additional benefits",
         "Fig31a_Fluvial_plus_Pluvial_plus_additional_{0}.png".format(winds),
          'pluv and fluv lt', ['add', 'pluv'])

    # (16) Get country-specific totals for fluv + pluv + additional
    #===============================================
    plot_lead_time_by_country(add_by_country, benefits,
        "Fluvial + pluvial-only + additional benefits (by country)",
        "Fig31b_Fluvial_and_Pluvial_plus_add_by_country.png", 6)

    # (6max) Get pluvial plus MAXIMUM fluvial
    #==================================
    plot_lead_time_benefits(maximum, benefits,
         "Maximum fluvial + pluvial-only",
         "Fig50a_Max_Fluvial_plus_Pluvial_only_{0}.png".format(winds),
         'pluv and fluv lt', ['pluv'])

    # (13max) Get country-specific totals for fluv + pluv
    #===============================================
    plot_lead_time_by_country(maxfluv_and_pluv, benefits,
        "Maximum fluvial + pluvial-only (by country)",
        "Fig50b_MaxFluvial_and_Pluvial_only_by_country.png", 6)

   # (17) Get pluvial plus fluvial + additional
    #==========================================
    plot_lead_time_benefits(maximum, benefits,
         "Fluvial + pluvial + additional benefits",
         "MaxFig61a_Fluvial_plus_Pluvial_plus_additional_{0}.png".format(winds),
         'pluv and fluv lt exec', ['maxadd', 'pluv'])

    # (18) Get country-specific totals for fluv + pluv + additional
    #===============================================
    plot_lead_time_by_country(max_add_by_country, benefits,
        "Fluvial + pluvial + additional benefits (by country)",
        "MaxFig61b_Fluvial_and_Pluvial_plus_add_by_country.png", 6)

    plt.show()

if __name__ == '__main__':
    main()

