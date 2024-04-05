# attenuation_calculator.py was written by William Byrne in the Fall of 2023
# To be used by DLR - Institute of Planetary Research
# Attenuaiton calculations based on Kalousova et al. 2017.
# Attenuation models based on Souckek et al. 2023

import math
import numpy as np
import pandas as pd


def calc_electrical_conductivity(temp, cl_conc):
    """
    calc_electrical_conductivity calculates the eletrical conductivity for ice at a given temeprature with a given concentration of cloride

    param temp: temperature at a given layer (in kelvin)
    param cl_conc: concentration of cloride at a given layer (in micromoles)
    """
    Kb = 8.617*(10**-5)     # boltzman constant(eV K^−1)
    Tr = 251                #refrence temperature

    𝜎0 = 7.2                # electrical conductivity of pure ice
    C0 = 1                  # Auxiliary constant
    E0 = 0.55               # Activation energy of pure ice

    𝜎1 = 0.43               # Molar electrical conductivity of Cloride
    E1 = 0.19               # Activation energy of Cloride


    i1 =(𝜎0*C0*math.exp((-(E0)/(Kb)*((1/temp) - (1/Tr))))) # pure ice
    i2 = (𝜎1*cl_conc*math.exp((-(E1)/(Kb)*((1/temp) - (1/Tr))))) # cloride
    return i1+i2




def calc_twoway_attenuation(depth_array, temp_array, attenuation_model):
    """
    calc_twoway_attenuation creates a dictionary containing the depth, temperature, and attenuation value at every layer of a 1-D temperature profile
    The attenuation values are calculated using one of theree attenuation models. 
    Attenuaiton models are based on material assumtptions of ice shell that impact electrical conductivity. Called from calc_electrical_conductivity

    param depth_array: numpy array containing the depth values at each layer (in meters)
    param temp_array: numpy array containing the temperature values at each layer (in kelvin)
    param attenuation_model: categorical variable(Either "high_loss", "mid_loss", or "low_loss") that determines which model to use for 2-way attenuation calculation
    return: dictionary contaning depth, temperature, and attenuation values
    """

    # sort depth/temp arrays in ascending order, so layer 0 is the surface layer
    # Assume that the surface layer has depth of 0.0 meters
    if depth_array[0] != 0.0:
        depth_array = np.flip(depth_array)
        temp_array = np.flip(temp_array)


    depth_array = list(depth_array)
    temp_array = list(temp_array)


    # columns 
    layer_array = []
    thickness_array = []
    if attenuation_model is not list:
        attenuation_model_array = [attenuation_model for model_name in range(len(depth_array))]
    
    for depth in depth_array:
        layer_num = depth_array.index(depth) + 1 
        layer_array.append(layer_num)
        if not layer_num == len(depth_array):
            thickness_array.append(abs(depth_array[layer_num-1] - depth_array[layer_num]))
        else:
            thickness_array.append(abs(depth_array[-2] - depth_array[-1]))
    


    summand_list = []
    attenuation_array = []
    for layer_num in layer_array:
        index = layer_num - 1

        # one-way attenuation calculaiton has layer thickness in kilometers
        if attenuation_model_array[index] == "low_loss":
            summand = 2*(0.914*calc_electrical_conductivity(temp_array[index], 0))*(thickness_array[index]/1000)#    summand = 2*(0.914*calc_electrical_conductivity(layer['temp'], 0))*(layer['thickness']/1000)
        elif attenuation_model_array[index] == "mid_loss":
            summand = 2*(0.914*calc_electrical_conductivity(temp_array[index], 150))*(thickness_array[index]/1000)#    summand = 2*(0.914*calc_electrical_conductivity(layer['temp'], 150))*(layer['thickness']/1000) 
        elif attenuation_model_array[index] == "high_loss":
            summand = 2*(0.914*calc_electrical_conductivity(temp_array[index], 300))*(thickness_array[index]/1000)#    summand = 2*(0.914*calc_electrical_conductivity(layer['temp'], 300))*(layer['thickness']/1000)
        else:
           raise Exception(f"The attenuation model {attenuation_model} is not supported. Please use one of the following: ['low_loss', 'mid_loss', 'high_loss']")
        
        summand_list.append(summand)
        attenuation_array.append(float(sum(summand_list)))

    layers_dict = {
            'layerNum':layer_array,
            'depth':depth_array,
            'thickness':thickness_array,
            'temp':temp_array,
            'attenuation_model':attenuation_model_array,
            'twoway_loss':attenuation_array
    }

   
    return layers_dict


def dict_to_df(layers_dict):
    """
    dict_to_df takes in the layers_dict dictionary and returns a Pandas DataFrame with the appropriate column names

    param layers_dict: python dictionary containing metadata about simulation layers and their attributes
    return: Pandas DataFrame with the appropriate column names with layer number set to index
    """

    df = pd.DataFrame(layers_dict)
    df.set_index('layerNum', inplace = True)

    return df