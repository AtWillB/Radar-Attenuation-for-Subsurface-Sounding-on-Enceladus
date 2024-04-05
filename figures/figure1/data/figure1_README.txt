# The following README.txt file describes the data organization for figure 1
## There were two different data storage formats used to generate figure 1


###### Temperature Profiles
	## Description - The .txt files in the specified folder describe the individual layers of 1-demensional geophysical simulation of Enceladu's internal  structure
	## Folder Name - temperature_profiles
	## Number of files - 11
	## File type - .txt
	## file format - white-space delimeted table
	## Columns
		# Note - Calculations start at the ice-ocean interface and go until the surface
		## Column 1 - Radius, in meters
		## Column 2 - Depth, in meters
		## Column 3 - Temperature, in Kelvin
	## Subfolders
		## Subfolder 1 - ice_shell_depth
			## Contains 3 different 1-dimensional simulations, in which the regolith thermal conductivity is 0.025 W/mK, the regolith thickness is 250 meters, but the ice shell depth is varied
		## Subfolder 2 - regolith_conductivity
			## Contains three different 1-dimensional simulations, in which the ice shell depth is 21 kilometers thick, there is 250 meters of regolith, but the regolith thermal conductivity is varied
		## Subfolder 3 - regolith_depth
			## Contains four different 1-dimensional simulations, in which the ice shell is 21 kilometers thick, the regolith thermal conductivity is 0.025 W/mK, but the regolith depth is varied




###### Attenuation Calculations
	## Description - the .csv viles in the specified folder describe the attenuation calculatuions for varied 1-dimensional temperature profiles of Enceladus. 
	## Folder Name - attenuation_calculations
	## Number of files - 6
	## File type - .csv
	## file format - comma delimeted table
	## Columns
		## index - the number of the current simulation layer, starting from the surface. 
		## Column 1 - depth_[m]: depth of current simulation layer in meters
		## Column 2 - thickness_[m]: thickess of current simultion layer in meters
		## Column 3 - temp_[k]: temperature of current simulation layer in kelvin
		## Column 4 - attenuation_model: type of attenuation model(high_loss, low_loss)
		## Column 5 - twoway_loss_[db]: 2-way attenuation of current simulation layer in decibals
		## Column6 - sim_name: name of simulation used. These will contain the name of the 1-dimenstional temperature profile used to perfrom the attenuation calculations
	##Subfolders
		## Subfolder 1 - regolith_condictivigy
			## Contains 2 different sets of attenuation calculations, performed on the different 1-dimensional temperature profile variations done on regolith thermal conductivity. The first set contains the attenuation calculations done for each layer of each regolith conductivity variation using the low loss attenuation profile, and the second set for the high loss profile
			## Varied Values
				## low loss
				## high loss
		## Subfolder 2 - regolith_thickness
			## Contains 2 different sets of attenuation calculations, performed on the different 1-dimensional temperature profile variations done on regolith thickness. The first set contains the attenuation calculations done for each layer of each regolith thickness variation using the low loss attenuation profile, and the second set for the high loss profile
			## Varied Values
				## low loss
				## high loss
		## Subfolder 3 - ice_shell_depth
			## Contains 2 different sets of attenuation calculations, performed on the different 1-dimensional temperature profile variations done on regolith thickness. The first set contains the attenuation calculations done for each layer of each regolith thickness variation using the low loss attenuation profile, and the second set for the high loss profile
			## Varied Values
				## low loss
				## high loss

