# temp_profile_util.py was written by William Byrne in the Fall of 2023
# To be used by DLR - Institute of Planetary Research
# Script used to create attenuation profile on top of 1-D temperature profile as a color gradient

import radar_attenuation.attenuation_calculator as ac

import numpy as np
import pandas as pd

from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns




def calc_attenuation_from_temp_profile(txt_path, attenuation_model, sim_name):
    """
    calc_attenuation_from_temp_profile takes the path to a txt file storing a 1-D temperature profile and creates a new DataFrame that includes attentuation calculations at every depth

    param txt_path: txt_path to where the temperature profile is stored
    param attenuation_model: attenuation model used for calulations. Can be either ['low_loss', 'mid_loss', 'high_loss']
    param sim_name: Name given the the temperature profile. Gathered from profle's filename
    return: Pandas DataFrame with the appropriate column names with layer number set to index
    """

    temp_profile = np.loadtxt(txt_path)

    depth_array = temp_profile[:, 1] 
    temp_array = temp_profile[: , 2]

    layers_dict = ac.calc_twoway_attenuation(depth_array, temp_array, attenuation_model)
    df = ac.dict_to_df(layers_dict)

    sim_name_column = [sim_name for i in range(df.shape[0])]
    df['sim_name'] = sim_name_column
    df['sim_name'] = df['sim_name'].astype(str)

    return df


def plot_line_segments(df, ax, linewidth = 5):
    """
    plot_line_segmetns plots a line segment between every point of temperature/depth. Creates a gradient change with attenuation. 
    This is done by descretizing all values below 1 and to 1 and above 100 to 100. 
    All values in between are sent to ints, and they correspond to a point on the provided colormap. 

    param df: pandas DataFrame of temperature profile that is used to generate line segments
    param ax: axes to plot the line segments. Provided from previous plt.subplots call
    param linewidth: the line width of all of the sengments for the provided profile. Defaults to 5 for single lineplots
    return: Pandas DataFrame with the appropriate column names with layer number set to index
    """

    df = df.reset_index()
    df['depth'] = df['depth']/1000
    df['twoway_loss'].where(df['twoway_loss'] > 1, 1.0, inplace=True)
    df['twoway_loss'].where(df['twoway_loss'] < 100.0, 100.0, inplace=True)
    df['twoway_loss'] = df['twoway_loss'].astype(int)
    cmap = sns.cubehelix_palette(100)
    
    for i in range(df['twoway_loss'].size - 1):
        x = [df['temp'][i], df['temp'][i + 1]]
        y = [df['depth'][i], df['depth'][i + 1]]
        color = cmap[df['twoway_loss'][i]-1]
        ax.plot(x, y, color=color, lw = linewidth)

    return df


def mk_single_attenuation_profile(txt_path, attenuation_model, ax, save = False, savepath = None):
    """
    mk_single_attenuation_profile is a wrapper function that takes a temperature profile and calculates attenuation for a given attenuation_model

    param txt_path: txt_path to where the temperature profile is stored
    param attenuation_model: attenuation model used for calulations. Can be either ['low_loss', 'mid_loss', 'high_loss']
    param ax: axes to plot the line segments. Provided from previous plt.subplots call
    param save: boolean detemining whether to save the DataFrame and figure or not. Defaults to False
    param savepath: path to save the figure and DataFrame. Defaults to None
    return: DataFrame generated from attenuation calculations of provided temperature profile
    """

    sim_name = txt_path.split("/")[-1].split('.')[0]
    df = calc_attenuation_from_temp_profile(txt_path, attenuation_model, sim_name)
    df_copy = df.copy(deep = True)
    df_copy = plot_line_segments(df_copy, ax)    
    lineplot_formatting(ax, df_copy['depth'].max(), multi = False)

    if save == True:
        df.rename(columns= {"depth": "depth_[m]", "thickness":"thickness_[km]", "temp" : "temp_[k]", "twoway_loss":"twoway_loss_[db]"}, inplace = True)
        df.to_csv(f"{savepath}/{sim_name}_{attenuation_model}.csv")
        plt.savefig(f"{savepath}/{sim_name}_{attenuation_model}.png")

    return df


def mk_multi_attenuation_profile(txt_list, attenuation_model, legend_title, linewidth_list, legend_labels, ax, save = False, savepath = None):
    """
    mk_multi_attenuation_profile is a wrapper function that multiple temperature profile and the attenuation at each depth for each profile

    param txt_path: txt_path to where the temperature profile is stored
    param attenuation_model: attenuation model used for calulations. Can be either ['low_loss', 'mid_loss', 'high_loss']
    param legend_title: describes the parameters being varied for multi-line plot
    param linewidth_list: list of linewideths assigned to each attenuation profile
    param legend_labels: list of labels for each attenuation profle. Denoted by linewidths
    param ax: axes to plot the line segments. Provided from previous plt.subplots call
    param save: boolean detemining whether to save the DataFrame and figure or not. Defaults to False
    param savepath: path to save the figure and DataFrame. Defaults to None
    return: DataFrame generated from attenuation calculations of provided temperature profile
    """

    df_list = []
    graph_list = []

    for txt_path, curr_linestyle in zip(txt_list, linewidth_list):
        sim_name = txt_path.split("/")[-1].split('.')[0]
        df = calc_attenuation_from_temp_profile(txt_path, attenuation_model, sim_name)
        df_list.append(df)

        graph_df = df.copy(deep = True)
        graph_df = plot_line_segments(graph_df, ax, curr_linestyle)
        graph_list.append(graph_df)

    big_df = pd.concat(df_list)
    big_graph_df = pd.concat(graph_list)
    lineplot_formatting(ax, big_graph_df['depth'].max(), legend_title, linewidth_list, legend_labels, multi = True)


    if save == True:
        legend_title = legend_title.lower().replace(" ", "_")
        big_df.rename(columns= {"depth": "depth_[m]", "thickness":"thickness_[km]", "temp" : "temp_[k]", "twoway_loss":"twoway_loss_[db]"}, inplace = True)
        big_df.to_csv(f"{savepath}/{legend_title}_{attenuation_model}.csv")
        plt.savefig(f"{savepath}/{legend_title}_{attenuation_model}.png")

    return big_df

    
def lineplot_formatting(ax, max_depth, legend_title = None, linewidths = [], legend_labels = [], multi = False):
    """
    lineplot_formatting provides consistent lineplot formatting to both single and multi lineplots greated with this script

    param ax: axes to plot the line segments. Provided from previous plt.subplots call
    param max_depth: maximum depth of the DataFrame used to make lineplot(s). Important for multi-lineplots when varying ice thickness
    param legend_title: title of the legend describing parameter variation. Defaults to None
    param linewidth: the line width of all of the sengments for the provided profile. Defaults to []
    param sim_name_list: list of legend labels used to show which parameters are varied. Defaults to []
    param multi: boolean describing whether the current plot is multi-line or single-line. Defaults to False
    return: None
    """

    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    ax.set_ylabel('Depth [km]', fontsize = 27) # Y label
    ax.set_xlabel('Temperature [K]', fontsize = 27) # X label
    ax.tick_params(axis="x", labelsize=24)
    ax.tick_params(axis="y", labelsize=24)
    ax.set_ylim(max_depth, 0)

    cmap = sns.cubehelix_palette(120)
    handles = []
    line1 = Line2D([0], [0], label='< 1 dB', color=cmap[0])
    line2 = Line2D([0], [0], label='40 dB', color=cmap[40])
    line3 = Line2D([0], [0], label='80 dB', color=cmap[80])
    line4 = Line2D([0], [0], label='> 100 dB', color=cmap[-1])
    handles = [line1, line2, line3, line4]


    if multi == True:
        style_handels = []
        for linestyle in linewidths:

            index = linewidths.index(linestyle)
            line = Line2D([0], [0], label=legend_labels[index], lw=linestyle, color=cmap[-1])
            style_handels.append(line)
        legend_1 = ax.legend(handles = style_handels,fontsize = 18, loc = 'center left', title = legend_title, title_fontsize = 18)
        ax.add_artist(legend_1)

    ax.legend(handles = handles,fontsize = 18, loc = "lower left", title = "2-way Attenuation", title_fontsize =  18)
    plt.grid()

    return


def calc_multiplot_stats(df, attentuation_model):
    """
    calc_multiplot_stats takes DataFrame generated via calc_attenuation_from_temp_profile and generates statistics showing how much ice shell was penetrated

    param df: Dataframe of 1-D temperature profile with attenuation calulation. Generated using calc_attenuation_from_temp_profile
    return: new Dataframe containing stats based on input DataFrame
    """

    max_depth_groupby = df.groupby('sim_name')[['twoway_loss', 'depth']].max()
    max_depth_groupby.rename(columns  = {'depth':'max_depth'}, inplace = True)

    groupby_df = df.query('`twoway_loss` <= 100').groupby('sim_name')[['twoway_loss', 'depth']].max().join(max_depth_groupby['max_depth'])
    groupby_df['pen_percent'] = groupby_df['depth']/groupby_df['max_depth']*100
    groupby_df.rename(columns  = {'depth':'pen_depth'}, inplace = True)
    groupby_df.reset_index(inplace = True)

    groupby_df['attenuation_model'] = [attentuation_model for x in range(len(groupby_df))]
    groupby_df.set_index(['sim_name', 'attenuation_model'])

    return groupby_df