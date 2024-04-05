import numpy as np
import warnings
warnings.filterwarnings('ignore')


### PROBLEM
## YOU NEVER EVEN USE RHO_WATER 
## OR RHO_SILICATE ONCE DURING THESE PROFILES



### NOTE
## These functions originate from the ConductiveProfile.ipynb file created by Ina Plesa

#   rho_solid = 920 # density at 270K from Hobbs 1974
#   rho_salt = 2160 # NaCl
#   concentration_salt = 0.05 # mass fraction
#   rho_mix = ((1-concentration_salt)/rho_solid+concentration_salt/rho_salt)**(-1) # mixture density
#   rhoc = 8000 # density of core, three-layer Model B (Schubert et al. 2009, Europa book)
#   rhom = 3500 # density of silicate mantle, three-layer Model B (Schubert et al. 2009, Europa book)#

#G = 6.674e-11 # Gravitational constant

#  radii & thickness
# Ro = 1560.8e3 # radius outer boundary
# Rc = 437e3 # radius of core, three-layer Model B (Schubert et al. 2009, Europa book)
# Rm = 1427e3 # radius of core, three-layer Model B (Schubert et al. 2009, Europa book)
# Ri = Ro-thickness # radius inner boundary

def calc_pressure(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Ro):
    """Function Description: Calulates and prints the pressure at ice-ocean interface assuming three-layered structure(Ice, Ocean, Silicate Core?)
    returns pressure(maybe kg s^-2 m^-1)
       param 1: rho_ice - Density of ice layer(kg/m^3)
       param 2: rho_water - density of water layer(kg/m^3)
       param 3: rho_silicate - density of silicate layer(kg/m^3)
       param 4: d_ice - depth of ice(meters)
       param 5: d_ocean - depth of ocean(meters)
       param 6: Ro - Radius of outer boundary(whole moon), meters)
       returns: pref - Pressure at ice-ocean boundary (maybe kg s^-2 m^-1)"""
    
    G = 6.674e-11
    g = 0.113 # Unused
    Ri = Ro - d_ice
    Rc = Ro - d_ice - d_ocean # Unused
    mass = (4/3)*np.pi*(Ro**3 - Ri**3)*rho_ice
    pref = G*mass*rho_ice/Ri
    
    # print('Pressure at ice-ocean interface: ', pref*1e-5, ' bar', pref*1e-6, ' MPa')
    return pref


def calc_TempIceOcean(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Ro):
    """Function Description: Calulates and prints the Temperature at ice-ocean interface
    retursn Temperature(k)
       param 1: rho_ice - Density of ice layer(kg/m^3)
       param 2: rho_water - density of water layer(kg/m^3)
       param 3: rho_silicate - density of silicate layer(kg/m^3)
       param 4: d_ice - depth of ice(meters)
       param 5: d_ocean - depth of ocean(meters)
       param 6: Ro - Radius of outer boundary(whole moon, meters)
       returns: Tb - Temperature at ice-ocean boundary (k)"""
    
    pref = calc_pressure(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Ro)
    Tb = 273.15*(pref/(-395.2e6)+1)**(1/9) # temperature lower boundary
    
    # print('Temperature at ice-ocean interface: ', Tb, ' K')
    return Tb


def calc_ConductiveTempProfile(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs, Ts, resolution):
    """Function Description: Calculates the temperature at each layer between the surface and ice/ocean boundary
    Uses Tc to calculate temp at each layer's depth
    returns list of depths(float, meters) and associated temps(k) at each depth
       param 1: rho_ice - Density of ice layer(kg/m^3)
       param 2: rho_water - density of water layer(kg/m^3)
       param 3: rho_silicate - density of silicate layer(kg/m^3)
       param 4: d_ice - depth of ice(meters)
       param 5: d_ocean - depth of ocean(meters)
       param 6: Rs - Surface Radius(meters)
       param 7: Ts - Surface temperature(k)
       returns: depth - depth of each layer from surface to ice-ocean boundary(meters)
       returns: temp - temperature at each layer from surface to ice-ocean boundary(k)"""
    

    Tc = calc_TempIceOcean(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs)

    
    Rc = Rs - d_ice # radius to ice shell
    n_shells = int(np.round((Rs - Rc) / resolution)) # couldnt you replace (Rs - Rc) with d_ice?
    depth = Rs - np.linspace(Rc, Rs, n_shells, endpoint=True) 
    # grab all depths by taking how many layers you will have, and get that many steps evenly in between Rc and Rs
    
    # T(r) = -rho*H*r**2/(6*k) - c1/(k*r) + c2
    # T(Rs) = Ts
    # T(Rc) = Tc
    temp = -1/(Rs - depth) * (Tc - Ts) * Rc * Rs / (Rc - Rs) + Ts + (Tc - Ts) * Rc / (Rc - Rs) # perform calculation for each depth in depth
      
    return temp, depth

"""Function Description: Calculates the temperature at each layer between the surface and ice/ocean boundary using conductivity(ice/regolith)
    Same as calc_ConductiveTempProfile, but instead uses conductivity(regolith for first 250 m, then ice) in temp per layer calculations
    returns list of depths(float, meters) and associated temps(k) at each depth"""
def calc_ConductiveTempProfileVark(klayer1, klayer2, Rinterface, rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs, Ts, resolution):
    """param 1:  klayer1 - Conductivity of regolith (W/mK)
       param 2:  klayer2 - Conduvity of ice (W/mK)
       param 3:  Rinterface - depth of regolith layer (meters)
       param 4: rho_ice - Density of ice layer(kg/m^3)
       param 5: rho_water - density of water layer(kg/m^3)
       param 6: rho_silicate - density of silicate layer(kg/m^3)
       param 7: d_ice - depth of ice(meters)
       param 8: d_ocean - depth of ocean(meters)
       param 9: Rs - Surface Radius(meters)
       param 10: Ts - Surface temperature(k)
       param 11: resolution - resolution of the model in meters
       returns: depth - depth of each layer from surface to ice-ocean boundary(meters)
       returns: temp - temperature at each layer from surface to ice-ocean boundary(k)"""

    Tc = calc_TempIceOcean(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs)
    
    Rc = Rs - d_ice
    n_shells = int(np.round((Rs - Rc) / resolution))
    radius = np.linspace(Rc, Rs, n_shells, endpoint=True)
    depth = Rs - radius
    conductivity = np.zeros(len(depth))
    temp = np.zeros(len(depth))

    # if @ regolith depth, use regolith conductivity. Else, use ice conductivity
    conductivity[radius > Rinterface] = klayer1
    conductivity[radius <= Rinterface] = klayer2
    
    dr = radius[1] - radius[0]
    
    mat = np.zeros((len(depth),len(depth)))
    mat[0][0] = 1
    mat[-1][-1] = 1
    
    rhs = np.zeros(len(depth))
    rhs[0] = Tc
    rhs[-1] = Ts
    
    for i in range(1, len(depth)-1):
        mat[i][i] = -1/(dr**2) * (conductivity[i]* radius[i]**2 + conductivity[i-1]* radius[i-1]**2)
        mat[i][i+1] = conductivity[i]* radius[i]**2 / (dr**2)
        mat[i][i-1] = conductivity[i-1]* radius[i-1]**2 / (dr**2)
        

    temp = np.linalg.solve(mat, rhs)
    
    return temp, depth


def calc_ConductiveTempProfileVarkT(klayer1, k0, Rinterface, rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs, Ts, resolution, steps):
    """Function Description: Calculates the temperature at each layer between the surface and ice/ocean boundary using conductivity(ice/regolith)
    Same as calc_ConductiveTempProfileVark, but instead has dynamic conducitvity change. simmulating increase in conductivity with pressure? What does steps ACTUALLY simulate??
    returns list of depths(float, meters) and associated temps(k) at each depth
       param 1:  klayer1 - Conductivity of regolith (W/mK)
       param 2:  k0 - Conduvity of ice (W/mK)
       param 3:  Rinterface - depth of regolith layer (meters)
       param 4: rho_ice - Density of ice layer(kg/m^3)
       param 5: rho_water - density of water layer(kg/m^3)
       param 6: rho_silicate - density of silicate layer(kg/m^3)
       param 7: d_ice - depth of ice(meters)
       param 8: d_ocean - depth of ocean(meters)
       param 9: Rs - Surface Radius(meters)
       param 10: Ts - Surface temperature(k)
       param 11: resolution - resolution of the model in meters
       param 12: steps - numbers of times you solve for temp and recursively update conductivity profiles as they change with depth

       returns: depth - depth of each layer from surface to ice-ocean boundary(meters)
       returns: temp - temperature at each layer from surface to ice-ocean boundary(k)"""

    Tc = calc_TempIceOcean(rho_ice, rho_water, rho_silicate, d_ice, d_ocean, Rs)
    
    Rc = Rs - d_ice
    n_shells = int(np.round((Rs - Rc) / resolution))
    radius = np.linspace(Rc, Rs, n_shells, endpoint=True)
    depth = Rs - radius
    conductivity = np.zeros(len(depth))
    temp = np.zeros(len(depth)) + Ts
    # bro what - for each radius, if the radius is greater than then inerface, then use klayer1(regolith conductivity relative to regolith)
    # Else, use replace current conductivty profile with k0/temp[radius <= Rinterface]
    # Now, conductivity will change after every step. 
    # Is this supposed to be like more and more detail in conductivity depending on the steps?
    conductivity[radius > Rinterface] = klayer1
    conductivity[radius <= Rinterface] = k0/temp[radius <= Rinterface]
    
    dr = radius[1] - radius[0]
    
    mat = np.zeros((len(depth),len(depth)))
    mat[0][0] = 1
    mat[-1][-1] = 1
    
    rhs = np.zeros(len(depth))
    rhs[0] = Tc
    rhs[-1] = Ts
    
    for i in range(steps):
        for i in range(1, len(depth)-1):
            mat[i][i] = -1/(dr**2) * (conductivity[i]* radius[i]**2 + conductivity[i-1]* radius[i-1]**2)
            mat[i][i+1] = conductivity[i]* radius[i]**2 / (dr**2)
            mat[i][i-1] = conductivity[i-1]* radius[i-1]**2 / (dr**2)
        # solving for every step? why? Some kind of temp aggregate?
        temp = np.linalg.solve(mat, rhs)
        conductivity[radius <= Rinterface] = k0/temp[radius <= Rinterface]
    
    return temp, depth


def write_data(filename, Ro, temp, depth, d_ice, conductivity_model):
    """Function Description: Write conducivity model data to a txt file with associated file header
   Returns nothing
    param 1: filename - path+name of file
    param 2: Ro - radial depth of ocean
    param 3: temp - temp profile of given model
    param 4: depth - depth profile of given model
    param 5: d_ice - depth of ice
    param 6: conductivity_model - follows pattern 'k_regolith = ' +  str(klayer1) + ' W/(mK), k_ice = ' + str(klayer2) + ' W/(mK), regolith thickness = ' + str(d_regolith) + ' m'"""
   
    caseID = 'thermal_40km_conductive'
    file1 = open(filename,"w")
    
    rad = Ro - depth

    file_header = ["# Temperature data: \n",
                   "# Ice shell thickness: " + str(d_ice) + " \n",
                   "# Conductivity model: " + conductivity_model + " \n",
                   "# Column 1: Radius in m \n",
                   "# Column 2: Depth in m \n",
                   "# Column 3: Temperature in K \n"]

    file1.writelines(file_header)

    for i in range(len(depth)):
        line=str(rad[i]) + "\t" + str(depth[i]) + "\t" + str(temp[i]) + "\n"
        file1.write(str(line))
                
    file1.close()


