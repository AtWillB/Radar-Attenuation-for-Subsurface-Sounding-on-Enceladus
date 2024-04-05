from gaia import Simulation

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.pyplot as plt



def load_sim(caseID, path_to_case_id, time_index = None, triangulation = False):

    s = Simulation(caseid=caseID, sim_dir=f"{path_to_case_id}{caseID}", grid_dir=f"{path_to_case_id}{caseID}", triangulation = False)
    return_list = (s)

    print("Case: " + s.caseid + " Loaded Grid: " + s.grid.gridfile)
    print("Got " + str(s.update()) + " new files, " + str(len(s.files))+ " in total.")

    if time_index:
        times = s.getTimes()

        o = s.getAtTime(times[time_index])
        print("Loaded time " + str(o.time) + " and iteration " + str(s.time_iter[o.time]))

        st=s.getStats()

        pps=str(s.grid.nCells/s.grid.nShells)
        print("cells per depth " + pps)
        return_list = (s, o)

    
    return return_list


def subplot_temp_composition_with_radius(o, R0, T0, scaleT):
    fig, ax = plt.subplots(1, 2, figsize=(20,10), sharey=True)
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1))

    oT = o.T
    Tprof = oT.getProfile()

    if hasattr(o, 'C'):
        oC = o.C
        Cprof = oC.getProfile()

    scaleD = R0/Tprof[-1,0]

    ax[0].plot(Tprof[:,2]*scaleT + T0,Tprof[:,0]*scaleD, color='blue',linewidth=4)

    ax[0].set_ylabel('Radius [km]', fontsize = 24.0) # Y label
    ax[0].set_xlabel('Temperature [K]', fontsize = 24) # X label
    ax[0].tick_params(axis="x", labelsize=20)
    ax[0].tick_params(axis="y", labelsize=20)
    ax[0].grid()
    ax[0].set_ylim(1560.8-scaleD,1560.8)

    ax[1].plot(Cprof[:,2],Cprof[:,0]*scaleD, color='blue',linewidth=4)

    ax[1].set_ylabel('Radius [km]', fontsize = 24.0) # Y label
    ax[1].set_xlabel('Composition', fontsize = 24) # X label
    ax[1].tick_params(axis="x", labelsize=20)
    ax[1].tick_params(axis="y", labelsize=20)
    ax[1].grid()
    ax[1].set_ylim(1560.8-scaleD,1560.8)

    plt.show()
    return


# Writes a file with 2d info from GAIA
# Data include z, y, and z coordinates, temperature in Kelvin, and Salt concentration
def write_2d(s, o, caseID, T0, scaleT, scaleD, scaleC, save_folder):

    filename = save_folder + caseID + "_2D_data_test.txt"
    file1 = open(filename,"w")
    Tprof = o.T.getProfile()

    file_header = ["# Case: \n", 
                   "# Column list: \n",
                   "# Column 1: x-coordinate in km \n",
                   "# Column 2: y-coordinate in km \n",
                   "# Column 3: z-coordinate in km \n",
                   "# Column 4: Temperature in K \n",
                   "# Column 5: Salt concentration in wt% \n"]

    file1.writelines(file_header)

    for i in range(len(o.T.data)):
        Tscaled = o.T.data[i]*scaleT + T0
        Xscaled = Tscaled*0.0
        if hasattr(o, 'C'):
            Xscaled = o.C.data[i]*scaleC/920*100*1.5
        x_coord = s.grid.coords[i][0] * scaleD
        y_coord = s.grid.coords[i][1] * scaleD
        z_coord = s.grid.coords[i][2] * scaleD

        rad = (x_coord**2 + y_coord**2 + z_coord**2)**(1/2)

        if(rad > Tprof[-1,0]*scaleD):
            x_coord = x_coord * Tprof[-1,0]*scaleD / rad
            y_coord = y_coord * Tprof[-1,0]*scaleD / rad
            z_coord = z_coord * Tprof[-1,0]*scaleD / rad

        if(rad < Tprof[0,0]*scaleD):
            x_coord = x_coord * Tprof[0,0]*scaleD / rad
            y_coord = y_coord * Tprof[0,0]*scaleD / rad
            z_coord = z_coord * Tprof[0,0]*scaleD / rad

        line=str(x_coord) + "\t" + str(y_coord) + "\t" + str(z_coord) + "\t" + str(Tscaled[0]) + "\t" + str(Xscaled[0]) + "\n"
        file1.write(str(line))

    file1.close()

    return


def write_1d_from_2d(o, caseID, T0, scaleT, scaleD, scaleC):
    filename = "../Data/" + caseID + "_1D_profile_test.txt"
    file1 = open(filename,"w")

    Tprof = o.T.getProfile()
    Cprof = o.C.getProfile()

    file_header = ["# Case: \n", 
                   "# Column list: \n","# Column 1: Radius in km \n","# Column 2: Temperature in K \n",
                   "# Column 4: Salt concentration in wt% \n"]

    file1.writelines(file_header)

    for i in range(len(Tprof[:,0])):
        rad = Tprof[i,0]*scaleD
        Tscaled = Tprof[i,2]*scaleT + T0
        Xscaled = Cprof[i,2]*scaleC/920*100*1.5

        line=str(rad) + "\t" + str(Tscaled) + "\t" + str(Xscaled) + "\n"
        file1.write(str(line))

    file1.close()
    return

def graph_2d_attenuation(o, graph_array, title, upperlim, lowerlim = 0, residual = False):
    x=o.grid.coords[:,0].reshape([o.grid.nShells, o.grid.nCellsPerShell[0]])
    y=o.grid.coords[:,1].reshape([o.grid.nShells, o.grid.nCellsPerShell[0]])

    # 
    fig = plt.figure( figsize=(2,7), dpi=500)

    if residual == True:
        count = 0
        for idx, val in np.ndenumerate(graph_array):
            graph_array[idx[0], idx[1]] =  (graph_array[idx[0], :].mean() - graph_array[idx[0], idx[1]])/graph_array[idx[0], :].mean()

        upperlim = 5
        lowerlim = -5

    T = graph_array
    cmap = sns.color_palette("rocket", as_cmap = True)


    

    conjtourplot = plt.contourf(x, y, T, np.linspace(lowerlim, upperlim,256), cmap=cmap, extend="both")
    ticksstep = upperlim/4
    cbarticks = np.arange(T.min() - (T.min()%ticksstep), T.max() - (T.max()%ticksstep) + ticksstep, ticksstep)
    cbarticks = np.arange(lowerlim, upperlim+ticksstep, ticksstep)
    cb = fig.colorbar(conjtourplot, orientation="horizontal", pad=0.02, shrink=0.8, ticks=cbarticks)
    cb.set_label(label=title, fontsize=8) #22 # 28
    plt.setp(cb.ax.xaxis.get_ticklabels(), fontsize=6.5) #20 # 26
    plt.axis('off')

    plt.show()

    return

def graph_2d(o, scaleT, T0, title, upperlim, lowerlim = 0):
    x=o.grid.coords[:,0].reshape([o.grid.nShells, o.grid.nCellsPerShell[0]])
    y=o.grid.coords[:,1].reshape([o.grid.nShells, o.grid.nCellsPerShell[0]])

    fig = plt.figure( figsize=(2,7), dpi=500)
    T = o.T.formatted_data[:,:,0]
    T = T * scaleT + T0
    cmap = sns.color_palette("inferno", as_cmap = True)


    conjtourplot = plt.contourf(x, y, T, np.linspace(lowerlim, upperlim,256), cmap = cmap, extend="both")
    ticksstep = 40
    cbarticks = np.arange(T.min() - (T.min()%ticksstep), T.max() - (T.max()%ticksstep) + ticksstep, ticksstep)
    cbarticks = np.arange(lowerlim, upperlim, ticksstep)
    cb = fig.colorbar(conjtourplot, orientation="horizontal", pad=0.02, shrink=0.8, ticks=cbarticks)
    cb.set_label(label="Temperature [K]", fontsize = 8) #22 # 28
    # plt.title(title)
    plt.setp(cb.ax.xaxis.get_ticklabels(), fontsize = 6.5) #20 # 26
    plt.axis('off')

    # plt.savefig("../data/test1.png")


def attenuation_vs_angratio(big_df, lim = False, a_r = 15):
    big_df_copy = big_df.copy(deep = True)

    big_df_copy['column'] = big_df_copy['column'].str.replace("column_", "")
    big_df_copy['column'] = big_df_copy['column'].astype(int)
    big_df_copy.set_index("column")

    num_col = len(big_df_copy['column'].unique())

    big_df_copy['a_r'] = big_df_copy['column'] * (a_r)/num_col


    if lim is False:
        print("ran false")
        twoway_loss = big_df_copy.groupby("a_r")['twoway_loss'].max().values
        a_r = big_df_copy.groupby("a_r")['twoway_loss'].max().index.array

        plt.plot(a_r, twoway_loss)
        plt.gca().ticklabel_format(useOffset=False)
        plt.title("Total Attenuation accross Angular Ratio")
        plt.xlabel("Angular Ratio [degrees]")
        plt.ylabel("Total 2-way Attenuation [dB]")
        plt.show()

    elif lim is True:
        depth_to_100 = big_df_copy.query("twoway_loss >= 60 and twoway_loss <= 100").groupby('a_r')['depth'].max().values
        a_r = big_df_copy.query("twoway_loss >= 60 and twoway_loss <= 100").groupby('a_r')['depth'].max().index.array

        # print(depth_to_100)
        print(type(depth_to_100))
        depth_to_100 = depth_to_100/1000
        depth_to_100 = depth_to_100.astype(float)
        depth_to_100 = np.round(depth_to_100, decimals=3)


        plt.plot(a_r, depth_to_100)
        plt.gca().ticklabel_format(useOffset=False)
        plt.title("Depth to 100 dB accross Angular Ratio")
        plt.xlabel("Angular Ratio [degrees]")
        plt.ylabel("Depth to 100 dB [km]")
        plt.show()
    return

def twod_profiles(df, proftype = 'temp'):
    df_copy = df.copy(deep = True)
    df_copy['layer_int'] = df_copy['layer'].str.replace('layer',"").astype(int)


    mins = df_copy.groupby('layer_int')[proftype].min().values
    mins_depths = df_copy.groupby('layer_int')['depth'].min().values.astype(float)/1000

    means = df_copy.groupby('layer_int')[proftype].mean().values
    means_depths = df_copy.groupby('layer_int')['depth'].min().values.astype(float)/1000

    maxs = df_copy.groupby('layer_int')[proftype].max().values
    maxs_depths = df_copy.groupby('layer_int')['depth'].min().values.astype(float)/1000

    min_label = r'$T_{Min}$'
    mean_label = r'$T_{Avg}$'
    max_label  = r'$T_{Max}$'
    x_label = 'Temperature [K]'
    if proftype == 'twoway_loss':
        min_label = r'$dB_{Min}$'
        mean_label = r'$dB_{Avg}$'
        max_label  = r'$dB_{Max}$'
        x_label = '2-way Attenuation [dB]'



    plt.plot(mins, mins_depths, label = min_label)
    plt.plot(means, means_depths, label = mean_label)
    plt.plot(maxs, maxs_depths, label = max_label)
    plt.legend(loc='lower left')
    plt.xlabel(x_label)
    plt.ylabel('Depth [km]')
    plt.grid()

    plt.gca().invert_yaxis()
    
    plt.show()
    return
