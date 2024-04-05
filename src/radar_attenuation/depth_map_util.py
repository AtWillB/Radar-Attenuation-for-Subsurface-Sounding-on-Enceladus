# depth_map_util.py was written by William Byrne in the Fall of 2023
# To be used by DLR - Institute of Planetary Research
# Script used to run attenuation calculations on different ice shell depth maps of Enceladus

import radar_attenuation.ConductiveProfile_equations as cpe
import radar_attenuation.attenuation_calculator as ac

import pandas as pd
import numpy as np
import struct
import time

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker




def hemingway_to_numpy(filepath):
    """
    hemingway_to_numpy takes in a filepath to one of Dr. Doug Hemingways Ice Shell depth maps(.tab)
    and converts the binary file into a numpy array of the shape (361, 181)
    The numpy shape downscaling is done by mean kernel averaging
   
    param filepath: file to one of Dr. Doug Hemingways Ice Shell depth maps. 
    return: a ice shell depth map as a numpy array of shape (361, 181)
    """

    resolution = 4.0
    num_rows = 720
    num_columns = 1440
    pixel_type = 'd'

    total_pixels = num_rows * num_columns

    with open(filepath, 'rb') as file:
        binary_data = file.read()

    pixel_format = f'>{total_pixels}{pixel_type}' # Read binary data in big endian format(>), read total_pixels number of values, and all values are floating point(pixel_type)
    pixels = struct.unpack(pixel_format, binary_data)

    depth_map = np.array(pixels) # flat
    depth_map = depth_map.reshape(1440, 720)
    depth_map = np.flip(depth_map).T # transpose to be correct orientation

    depth_map = downscale_hemingway(depth_map)


    return  depth_map


def downscale_hemingway(hemingway_depth_map):
    """
    downscale_hemingway takes heingway_depth_map(numpy array) and downscales it from (1440, 720) to (361, 181)
   
    param hemingway_depth_map: numpy array of ice shell depths in the shape (1440, 720) 
    return: a ice shell depth map as a numpy array of shape (361, 181)
    """

    # downsize to (360, 180)
    scaling_factors = (4, 4)
    new_shape = (hemingway_depth_map.shape[0] // scaling_factors[0], hemingway_depth_map.shape[1] // scaling_factors[1])
    downscaled_hem_depth_map = hemingway_depth_map[:new_shape[0]*scaling_factors[0], :new_shape[1]*scaling_factors[1]].reshape(new_shape[0], scaling_factors[0], new_shape[1], scaling_factors[1]).mean((1, 3))

    # Add a extra column and row to be (361, 181)
    downscaled_hem_depth_map = np.vstack((downscaled_hem_depth_map, downscaled_hem_depth_map[-2:-1]))
    rightmost_column = [[x] for x in downscaled_hem_depth_map[:, -1]]
    downscaled_hem_depth_map = np.hstack((downscaled_hem_depth_map, np.array(rightmost_column)))
    downscaled_hem_depth_map = np.flip(downscaled_hem_depth_map, axis = 1)
    downscaled_hem_depth_map = np.roll(downscaled_hem_depth_map, 180, axis = 1)
    downscaled_hem_depth_map = downscaled_hem_depth_map*1000

    return downscaled_hem_depth_map


def create_coords(rows, columns):
    """
    create_coords creates x and y 2-D arrays of shape (361, 181)
    Each cell (for array x and array y respectively) represent the coordinate location of values inputed in plt.pcolormesh
   
    param rows: number of rows for x and y coordinatea
    param columns: number of colums for x and y coordinatea
    return: a touple containing coordinate arrays
    """

    big_array = np.linspace(-1*rows+1, rows-1, columns)

    for i in range(0, rows-1):
        new_array = np.linspace(-1*rows+1, rows-1, columns)
        big_array = np.vstack((big_array, new_array))

    x = big_array.reshape(rows, columns)

    
    big_array = np.linspace((rows-1)/2, (rows-1)/2, columns)

    for i in range(1, rows):
        new_array = np.linspace((rows-1)/2-i, (rows-1)/2-i, columns)
        big_array = np.vstack((big_array, new_array))

    y = big_array.reshape(rows, columns)

    return (x, y)


def mk_depth_map(depth_map, x, y, cmap_name, cbar_label, levels, fig, tick_fontsize, ax = None, atten_limit = None, depth_map_limit = None):
    """
    plot_depth_map pltots depth maps(both attenuation and ice shell) using plt.pcolormesh, and formats it given the specified parameters

    param depth_map: ice shell depth/ attenuation map as a 2-d Numpy array
    param x: a array of same shape as depth_map containing all x coordinates for depth_map values
    param y: a array of same shape as depth_map containing all y coordinates for depth_map values
    param cmap_name: name of colormap pallete to be used
    param cbar_label: Label on colorbar
    param levels: levels used for contour plots(if plotting with contourf)
    param fig: figure object generated from plt.subplots() call
    param tick_fontsize: fontsize of all ticks on figure
    param ax: axes object generated from plt.subplots() call. Defaults to None
    param limit: value limit on depth map. Anything above this value is set to this value. Defaults to None
    return: colorbar object from 
    """

    temp_depth_map = np.copy(depth_map)

    if ax == None:
        ax = plt.gca()
    # if a limit is set on depth map, use pcolormesh and get exact color values back. If not, use contourf
    if atten_limit is not None:
        temp_depth_map[temp_depth_map >= atten_limit] = atten_limit
        norm = Normalize(vmin = 0, vmax = atten_limit)
        ax.pcolormesh(x, y, temp_depth_map, cmap = cmap_name, vmin = 0, vmax = atten_limit)#, levels = levels)
        
    if depth_map_limit is not None:
        norm = Normalize(vmin = 0, vmax = depth_map_limit)
        ax.contourf(x, y, temp_depth_map, cmap = cmap_name, vmin = 0, levels = levels)
    if cmap_name != "cividis":
        cmap_name  = "rocket_r"


    ax.grid()
    cmap = plt.get_cmap(cmap_name, levels)
    cbar = fig.colorbar(ScalarMappable(norm= norm, cmap=cmap), ax = ax)
    cbar.set_label(cbar_label, fontsize = tick_fontsize)
    cbar.ax.tick_params(labelsize=tick_fontsize-2)
    plt.xlabel("Longitude [deg]")
    x_ticks = ["120 °W", "60 °W", "0°","60 °E","120 °E"]
    ax.set_xticks([-120, -60, 0, 60, 120])
    ax.set_xticklabels(x_ticks) 


    ax.set_ylabel("Latitude [deg]")
    y_ticks = ["60 °S", "30 °S", "0°","30 °N","60 °N"]
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_yticklabels(y_ticks)

    return cbar
    

def write_map(filepath, attenuation_map, depth_map_title, regolith_map_title, attenuation_model, unit = 'attenuation'):
    """
    write_map writes a calculated attenuation map as a csv to a txt file. 
    Parameters provide comment headers and filename structure

    param attenuation_map: attenuation map to be saved
    param depth_map_title: title of depth map used for attenuation caluclation(ex - hemingway_map_a)
    param regolith_map_title: title of regolith distribution used for attenuation caluclation(ex - exp(exponential))
    param attenuation_model: attenuation model used for attenuation caluclations
    param unit: datatype of value saved. Either in dB, in meters as Depth to 100 dB, or in % as relative depth penetration
    return: None
    """

    filename = f'{depth_map_title}_{regolith_map_title}_{unit}_{attenuation_model}.csv'
    
    if unit == 'attenuation':
        unit = '2-way Attenuation [dB]'
    elif unit == 'depth':
        unit = 'Depth to 100 dB [m]'
    elif unit == "relative":
        unit = "Relative Depth Penetration [%]"


    header = '# This file contains the attenuation results for a ice shell depth map of Enceladus\n'
    depth_map_title = f'# The tile of the attenuation map is: {depth_map_title}\n'
    regolith_map_title = f'# The title of the regolith map used was: {regolith_map_title}\n'
    unit = f'# The unit of measurement captured in each cell of this map was: {unit}\n'
    shape = f'# The shape of this file is {str(attenuation_map.shape)}\n'
    attenuation_model = f'# The loss type of this attenuation map is : {attenuation_model}\n'
    

    comments = header + depth_map_title + regolith_map_title + unit + shape + attenuation_model
    df = pd.DataFrame(attenuation_map)

    with open(filepath+filename, 'a') as f:
        f.write(comments)
        df.to_csv(f, index = False, header = False)
    
    return


def create_regolith_maps():
    """
    create_regolith_maps creates a dictionary of each regolith distribution map as a numpy array.
    Will each be a numpy array of shape (361, 181)
    Three kinds of maps are no, exp(exponential), and bc(best-case)
    All regolith map deposition values are in meters
   
    return: a dictionary of all regolith distributions
    """
    no_regolith_map = np.zeros(shape=(181, 361))

    exp_lats = [x for x in range(-90, 91)]
    exp_lats.reverse()
    exp_reg_depths = [pow(1.0073, x)*250 - 230 for x in range(0, 181)]
    exp_regolith_depth_map = [np.linspace(x, x, 361) for x in exp_reg_depths]
    exp_regolith_depth_map = np.array(exp_regolith_depth_map)

    lats = [x for x in range(-90, 91)]
    reg_depths = [0 for x in range(0, 6)]+ [pow(250, x/85) for x in range(86, 0, -1)]+[0 for x in range(0, 89)]
    reg_depths[84] = 1.6
    reg_depths[85] = 1.4
    reg_depths[86] = 1.2
    reg_depths[87] = 1.0
    reg_depths[88] = .8
    reg_depths[89] = .6 
    reg_depths[90] = .4
    reg_depths[91] = .2
    reg_depths[92] = 0
    reg_depths.reverse()
    bc_regolith_map = [np.linspace(x, x, 361) for x in reg_depths]
    bc_regolith_map = np.array(bc_regolith_map)

    regolith_map_dict = {}
    regolith_map_dict['no'] = no_regolith_map
    regolith_map_dict['exp'] = exp_regolith_depth_map
    regolith_map_dict['bc'] = bc_regolith_map

    return regolith_map_dict


def create_attuenation_map_data(depth_map, regolith_depth_map, attenuation_model = 'low_loss', datatype = 'attenuation'):
    """
    create_attuenation_map_data 

    param depth_map: numpy array of ice shell depth values of shape (361, 181)
    param regolith_depth_map: numpy array of regolith depth values of shape (361, 181)
    param attenuation_model: attenuation model used for attenuation caluclations
    param datatype: whether to store the attenuation calculations in dBs, as depth to 100 dB in meters, as attenuation @ eutectic_NH4CL, or attenuation @ eutectic_NH4CL
    return: a map with all depth values replaced with the proper attenuation datatype, or None if incorrect datatype specified
    """
    
    temp_depth_map= np.copy(depth_map)

    Ro = 252.1e3 # radius of Enceladus in meters
    d_ocean = 37e3 # ocean thickness in m
    Ts = 60 # surface temperature in K
    resolution = 50 # model resolution in m
    density_ice = 925 # kg/m3
    density_water = 1050 # kg/m3
    density_silicate = 2400 # kg/m3

    klayer1 = 0.025
    k0 = 651

    last_layer = ''
    start = time.time()

    # for every ice shell depth value, calculate that cells attenuation value and replace the depth value with attenuation calculation
    for idx, d_ice in np.ndenumerate(temp_depth_map):

        reg_depth = regolith_depth_map[idx[0], idx[1]]
        Rinterface = Ro - reg_depth


        temperature_varK_regolith0, depth = cpe.calc_ConductiveTempProfileVarkT(klayer1, k0, Rinterface, density_ice, density_water, density_silicate, d_ice, d_ocean, Ro, Ts, resolution, 5)
        layers_dict = ac.calc_twoway_attenuation(depth, temperature_varK_regolith0, attenuation_model) #.calc_A2(depth, temperature_varK_regolith0, loss_type)
        last_layer = len(layers_dict)-1
        if datatype == 'attenuation':
            temp_depth_map[idx[0], idx[1]] = layers_dict['twoway_loss'][last_layer]
        elif datatype == 'depth':
            twoway_loss_list = [x for x in range(0, len(layers_dict['layerNum'])) if layers_dict['twoway_loss'][x] >= 100]
            if len(twoway_loss_list) == 0:
                temp_depth_map[idx[0], idx[1]] = layers_dict['depth'][last_layer]
            else:
                depth_to_100_index = twoway_loss_list[0]
                temp_depth_map[idx[0], idx[1]] = layers_dict['twoway_loss'][depth_to_100_index]
        elif datatype == 'eutectic_NH4CL_attenuation':
            twoway_loss_list = [x for x in range(0, len(layers_dict['layerNum'])) if layers_dict['temp'][x] >= 257.79]
            index_of_temp_at_NH4CL = twoway_loss_list[0]
            temp_depth_map[idx[0], idx[1]] = layers_dict['twoway_loss'][index_of_temp_at_NH4CL]
        elif datatype == 'eutectic_NH4CL_depth':
            twoway_loss_list = [x for x in range(0, len(layers_dict['layerNum'])) if layers_dict['temp'][x] >= 257.79]
            index_of_temp_at_NH4CL = twoway_loss_list[0]
            temp_depth_map[idx[0], idx[1]] = layers_dict['depth'][index_of_temp_at_NH4CL]
        elif datatype == 'eutectic_NH3_attenuation':
            twoway_loss_list = [x for x in range(0, len(layers_dict['layerNum'])) if layers_dict['temp'][x] >= 175.45]
            index_of_temp_at_NH3 = twoway_loss_list[0]
            temp_depth_map[idx[0], idx[1]] = layers_dict['twoway_loss'][index_of_temp_at_NH3]
        elif datatype == 'eutectic_NH3_depth':
            twoway_loss_list = [x for x in range(0, len(layers_dict['layerNum'])) if layers_dict['temp'][x] >= 175.45]
            index_of_temp_at_NH3 = twoway_loss_list[0]
            temp_depth_map[idx[0], idx[1]] = layers_dict['depth'][index_of_temp_at_NH3]
        else:
                print("The datatype you have entered is not supported. Please use a supported datatype: ['attenuation', 'depth', 'eutectic_NH4CL']")
                return None

  

        if idx[0] % 20 == 0 and (idx[0] != 0 and idx[0] != 1) and idx[1] == 0:
            print(f'Finished latidue {idx[0]} at minute {abs(start - time.time())/60}')
        
    return temp_depth_map


def depth_map_formatting(ax, cbar = False, y_ticks = False, x_ticks = False, label_fontsize = 14):
    """
    depth_map_formatting provides utility for turning on/off different parts of a depth maps

    param ax: axis variable from plt.subplots() call
    param cbar: Boolean, determins whether or not to display the colorbar. Defaults to False
    param y_ticks: Boolean, determins whether or not to display y_ticks and label. Defaults to False
    param x_ticks: Boolean, determins whether or not to display x_ticks and label. Defaults to False
    param label_fontsize: Font size of custom colorbar label. Defaults to 14
    return: None
    """

    x_label = 'Longitude [deg]'
    y_label = "Latitude [deg]"
    tick_fontsize  = label_fontsize - 3

    if not y_ticks:
        ax.set_yticklabels([])
        ax.set_ylabel('')
    else:
        ax.set_ylabel(y_label, fontsize = label_fontsize)
        ax.tick_params(axis="y", labelsize=tick_fontsize)

    
    if not x_ticks:
        ax.set_xticklabels([])
        ax.set_xlabel('')
    else:
        ax.set_xlabel(x_label, fontsize = label_fontsize)
        ax.tick_params(axis="x", labelsize=tick_fontsize)
    
    if cbar:
        cbar.remove()
    
    return

def mk_custom_cbar(cbar_ax, fig, bar_label = "Relative Pen. Depth [%]", label_fontsize = 14):
    """
    mk_custom_cbar creates colorbar with manual dimensions on a new axis based on cbar_ax

    param cbar_ax: dimensions of new axis for custom colorbar
    param fig: Figure variable generated from plt.subplots() call
    param cbar_label: Label name for custom colorbar
    param label_fontsize: Font size of custom colorbar label. Defaults to 14
    return: None
    """

    norm = plt.Normalize(0, 100)
    rocket_colors = sns.color_palette('rocket', 80) + [(0.9955451416666666, 0.97512059, 0.960499675)]
    cmap = ListedColormap(sns.color_palette(rocket_colors, as_cmap=True))

    sm = plt.cm.ScalarMappable(cmap= cmap, norm=norm)
    sm.set_array([])

    cbar_ax = fig.add_axes(cbar_ax)
    cb = fig.colorbar(sm, orientation="vertical", pad=0.02, shrink=0.8, cax = cbar_ax)
    cb.set_label(bar_label, fontsize = label_fontsize)

    tick_locator = ticker.MultipleLocator(base=10)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_ticklabels([0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize = label_fontsize - 4)

    return

def calc_depth_map_stats(np_heatmap, name_of_depthmap, name_of_reg_dist, attenuation_model):
    """
    calc_depth_map_stats takes DataFrame generated via mk_depth_map and generates statistics

    param heatmap_df: DataFrame of relative penetration percentages. Generated using mk_depth_map
    return: new Dataframe containing stats based on input DataFrame
    """


    global_min = round(np_heatmap.min(), 3)
    global_max = round(np_heatmap.max(), 3)
    global_mean  = round(np_heatmap.mean(), 3)

    mean_at_equator = round(np_heatmap[88:93, :].mean(), 3)
    mean_at_north_pole = round(np_heatmap[0:5, :].mean(), 3)
    mean_at_south_pole = round(np_heatmap[-6:-1, :].mean(), 3)
    percent_at_ocean = round((np.where(np_heatmap == 100.0)[0].shape[0]/np.prod(np_heatmap.shape))*100, 3)

    df_dict = {'icesehll_map': name_of_depthmap,
               'regolith_distribution': name_of_reg_dist,
               'attenuation_model': attenuation_model,
               'global_min': global_min,
               'global_max': global_max,
               'global_mean': global_mean,
               'mean_at_north_pole': mean_at_north_pole,
               'mean_at_equator': mean_at_equator,
               'mean_at_south_pole': mean_at_south_pole,
                'percent_at_ocean': percent_at_ocean
               }
    stats_df = pd.DataFrame(df_dict, index = [1])


    return pd.DataFrame(stats_df)

def create_depth_pen(depth_to_100_map, depth_map):
    """
    create_depth_pen takes in a depth map containing values in meters of depth to 100 dB, and converts the values to a percentage

    param depth_to_100_map: a depth map containing values in meters of depth to 100 dB
    param depth_map: ice shell depth map from Dr. Hemingway
    return: new numpy array relative depth percentages instead of meters to 100 dB
    """

    return (np.round(depth_to_100_map, 3)/np.round(depth_map, 3))*100