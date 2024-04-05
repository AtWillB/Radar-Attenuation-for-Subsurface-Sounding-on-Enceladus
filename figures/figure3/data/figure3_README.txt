# The following README.txt file describes the data organization for figure 3a and 3b
## There were three different data storage formats used to generate figure 3a and 3b


#### Ice Shell Depth maps
    ## Description - The csv files in the specified folder contains values from a ice shell depth map originating from Hemingway et al. 2019. 
    ## Folder Name - ice_shell_depth_map
    ## Number of files - 2
    ## File type - .csv
    ## File format - comma delimited table, no column names, no index
    ## Columns
        # Each column contains all ice shell depth values associated with that entire degree of longitude. 
        # This is likewise for rows and latitude
        # Each cell is the ice shell depth calculated at that degree of longitude and latitude
    ## Subfolders - none

#### Regolith Distribution maps
    ## Description - The csv files in the specified folder contains values from a regolith shell depth map created for figures 3a and 3b
    	## Folder Name - regolith_diststribution_maps
    	## Number of files - 3
    	## File type - .csv
    	## File format - comma delimited table, no column names, no index
    	## Subfolders - None
    	## Columns
            # Each column contains all regolith depth values associated with that entire degree of longitude. 
            # This is likewise for rows and latitude
            # Each cell is the regolith depth calculated at that degree of longitude and latitude



#### Attenuation Calculation maps
    ## Description - The csv files in the specified folder contains calculated relative depth peneteration percentages for a variety of ice shell scenarios 
    	## Folder Name - attenuation_calculations_maps
    	## Number of files - 12
    	## File type - .csv
    	## File format - comma delimited table, no column names, no index
    	## Subfolders
            ## Subfolder 1
                ## contains relative depth penetration percentage values for all combinations of attenuation regolith distributions, and depth maps that used the low loss attenuation profle
            ## Subfolder 2
                ## contains relative depth penetration percentage values for all combinations of attenuation regolith distributions, and depth maps that used the high loss attenuation profle
    	## Columns
            # Each column contains all relative depth penetration percentages associated with that entire degree of longitude. 
            # This is likewise for rows and latitude
            # Each cell is the relative depth penetration calculated at that degree of longitude and latitude