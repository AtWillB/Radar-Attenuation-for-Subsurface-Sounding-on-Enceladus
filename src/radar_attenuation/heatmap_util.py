# heatmap_util.py was written by William Byrne in the Fall of 2023
# To be used by DLR - Institute of Planetary Research
# Script used to create heatmaps as combinations of parametric variations of temperature profiles


import radar_attenuation.attenuation_calculator as ac
import radar_attenuation.temp_profile_util as tpu

import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt




 # This function contains lot of requirments, but there is was no research reason to make this modular
def mk_heatmap(ice_depth, attenuation_model, txt_folder_path, attenuation_limit = 100, save = False, savepath = False): 
    """
    mk_heatmap creates a DataFrame used to generate heatmaps via seaborn.heatmap()
    Heatmap values represent relative depth penetration as a percent based on the input ice shell depth
    Requires 15 1-D temperature profiles in the same folder, all with the same ice shell depth but varying regolith thermal conductivity and regolith depth
    Variations must be identical to variations described in cond_list and reg_depth_list
    Files of every 1-D temperature profile must be named using the pattern f"Regolith_{reg_depth}_icedepth_{ice_depth}_conduct_{cond}"

    param ice_depth: ice shell depth of all 15 1-D temperature profiles
    param attenuation_model: attenuation model used for calulations. Can be either ['low_loss', 'mid_loss', 'high_loss']
    param txt_folder_path: path to folder where the 15 1-D temperature profiles are stored
    param save: boolean detemining whether to save the DataFrame and figure or not. Defaults to False
    param savepath: path to save the figure and DataFrame. Defaults to None
    return: Dataframe to be used to generate heatmap. 
    """

    heatmap_dict = {}
    cond_list = ["01", "0025", "0001"]
    reg_depth_list = ["0", "100", "250", "400", "700"]

    print(f"Ice Shell Depth: {ice_depth} km, Loss Type: {attenuation_model}")

    for cond in cond_list:
        max_list = []
        for reg_depth in reg_depth_list:
            fname = f"Regolith_{reg_depth}_icedepth_{ice_depth}_conduct_{cond}"
            max = tpu.calc_attenuation_from_temp_profile(f'{txt_folder_path}{fname}.txt', attenuation_model, fname).query("twoway_loss >= @attenuation_limit")["depth"].min()
            max_list.append(max)

        max_list.reverse()
        max_list = [((x/1000)/ice_depth)*100 for x in max_list]
        heatmap_dict[str(cond)+" W/mK"] = max_list
    
    heat_df = pd.DataFrame(heatmap_dict)
    heat_df.rename(columns = {"01 W/mK": "0.1 W/mK", "0025 W/mK":"0.025 W/mK", "0001 W/mK":"0.001 W/mK"}, inplace = True)
    heat_df["Regolith Depth [m]"] = ['700', '400', '250', '100', '0']
    heat_df.set_index('Regolith Depth [m]', inplace=True)
    heat_df.fillna(0, inplace=True)
    heat_df["0.1 W/mK"].replace(0, 100, inplace = True)
    heat_df["0.025 W/mK"].replace(0, 100, inplace = True)
    heat_df["0.001 W/mK"].replace(0, 100, inplace = True)

    if save == True:
        heat_df.to_csv(savepath+f"/heatmap_{ice_depth}_{attenuation_model}.csv")
        
    return heat_df


def mk_custom_cbar(cbar_ax, fig, cbar_label, ticklabels, label_fontsize = 14):
    """
    mk_custom_cbar creates colorbar with manual dimensions on a new axis based on cbar_ax

    param cbar_ax: dimensions of new axis for custom colorbar
    param fig: Figure variable generated from plt.subplots() call
    param cbar_label: Label name for custom colorbar
    param ticklabels: list of tick labels for the colorbar
    param label_fontsize: Font size of custom colorbar label. Defaults to 14
    return: None
    """

    tick_fontsize  = label_fontsize - 3

    cbar_ax = fig.add_axes(cbar_ax)
    norm = plt.Normalize(0, 100)
    sm = plt.cm.ScalarMappable(cmap='rocket', norm=norm)
    sm.set_array([]) # For istannce, I don't remeber why this needs to be empty, but it does
    cb = fig.colorbar(sm, orientation="horizontal", pad=0.02, shrink=0.8, cax = cbar_ax)
    cb.set_label(cbar_label, fontsize = label_fontsize)
    cb.set_ticklabels(ticklabels, fontsize = tick_fontsize)

    return


def heatmap_formatting(ax, y_ticks = False, x_ticks = False, cbar = False, label_fontsize = 14):
    """
    heatmap_formatting provides utility for turning on/off different parts of a heatmap

    param ax: axis variable from plt.subplots() call
    param y_ticks: Boolean, determins whether or not to display y_ticks and label. Defaults to False
    param x_ticks: Boolean, determins whether or not to display x_ticks and label. Defaults to False
    param cbar: Boolean, determins whether or not to display the colorbar. Defaults to False
    param label_fontsize: Font size of custom colorbar label. Defaults to 14
    return: None
    """

    x_label = "Porous Layer Conductivity [W/mK]"
    y_label = "Porous Layer Thickness [m]"
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
    
    if not cbar:
        cbar = ax.collections[0].colorbar
        cbar.remove()
    
    return


def custom_cmap():
    """
    custom_cmap generates a custom colormap in format acceptable by Matplotlib methods
    Discrete colormap with 20 colors, but the last color is 5x as bright as bightest value in "rocket" colormap with 20 colors

    return: custom colormap
    """

    colors = sns.color_palette('rocket', 19) + [(0.9910902833333334, 0.95024118, 0.92099935)]
    cmap = sns.color_palette(colors)

    return cmap


def calc_heatmap_stats(heatmap_df, iceshell_depth, attenuation_model):
    """
    calc_heatmap_stats takes DataFrame generated via mk_heatmap and generates statistics

    param heatmap_df: DataFrame of relative penetration percentages. Generated using mk_heatmap
    return: new Dataframe containing stats based on input DataFrame
    """

    min = heatmap_df.to_numpy().min()
    max = heatmap_df.to_numpy().max()
    mean = heatmap_df.to_numpy().mean()
    num_to_ocean = heatmap_df.query("@heatmap_df == 100").count().sum()

    df_dict = {'iceshell_depth': iceshell_depth, 'attenuation_model': attenuation_model, 'min':min, 'max':max, 'mean':mean, 'num_to_ocean': num_to_ocean}

    heatmap_dict_stats = pd.DataFrame(df_dict, index = [1])
    heatmap_dict_stats.set_index(['iceshell_depth', 'attenuation_model'], inplace = True)




    return pd.DataFrame(heatmap_dict_stats)
