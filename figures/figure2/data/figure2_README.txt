# The following README.txt file describes the data organization for figure 2
## There were two different data storage formats used to generate figure 2


#### Temperature Profiles
	## Description - A list of 45 different 1-D temperature profiles, used to generate the 6 3x3 heat maps in figure 2. Each temperature profile has one of the following parameters varied: regolith depth, ice shell depth, regolith thermal conductivity. Filename indicates the profiles unique parameter combination. 
	## Folder Name - temperature_profiles
	## Number of files - 45
	## File type - .txt
	## File format - white-space delimited table
	## Columns
		## Column 1 - Radius, in meters 
		## Column 2 - Depth, in meters
		## Column 3 - Temperature, in Kelvin
	## Subfolders - None





#### Attenuation Calculations
	## Description - A list of 6 files, each containing the data used to generate one of the heat maps in figure 2. File name indicates the ice shell depth and loss type of given profile combinations. Each value(besides the column and the index) contain the depth penetration percentage of that simulations ice shell depth. 100.0 indicates a simulation in which the radar penetrated 100% of the ice shell(in other words, the attenuation never reached 100 dB). 
	## Folder Name - attenuation_calculations
	## Number of files - 6
	## File type - .csv
	## File format - comma delimited table with column names, with index
	## Columns
		Note - Values in columns with name "[...] W/mK" DO NOT indicate profiles regolith thermal conductivity. Those column names are used to label the heat map. All values besides the first column indicate the depth penetration as a percentage of that simulations total depth (either 5, 21, or 35 km)
		## Column 1 - Regolith depth, in meters
		## Column 2 - Depth penetration percentage for regolith thermal conductivity 0.1 W/mK
		## Column 3 - Depth penetration percentage for regolith thermal conductivity 0.025 W/mK 
		## Column 4 - Depth penetration percentage for regolith thermal conductivity 0.001 W/mK 
	## Subfolders - None




