from pocra_models import *

# The following simulations are built from scratch.

# weather array inputs are to be given either as python lists
# or as CSV files with apropriate weather parameters as headers.

# Note that there are 3 ways in which et0 may be supplied/derived:
# 1. Provide a list of et0 values for each day of simulation
# 2. Provide equal length (same as the simulation duration)
#	lists of temp_min, temp_avg, temp_max and r_a;
#	where r_a is the solar-radiation
# 3. Provide equal length (same as the simulation duration)
#	lists of temp_min, temp_avg and temp_max and
#	the latitude of the location. The latter will determine
#	the solar-radiation

# example 1:
# give et0 directly
psmm1 = PocraSMModelSimulation(
	soil_texture='clayey', soil_depth_category='deep (50 to 100 cm)',
	lulc_type='cropped in two seasons', slope=3,
	rain=[4]*122+[0]*243,
	et0=[6]*365,
	crop='bajri'
)

# example 2:
# give temperature data and solar-radiations for dynamic et0
psmm2 = PocraSMModelSimulation(
	soil_texture='clayey', soil_depth_category='deep (50 to 100 cm)',
	lulc_type='cropped in two seasons', slope=3,
	rain=[4]*122+[0]*243,
	r_a=[40]*365, # solar radiation
	temp_min=[22]*365, temp_avg=[28]*365, temp_max=[34]*365,
	crop='bajri'
)

# example 3:
# give temperature data and latitude(which determines solar-radiation) for dynamic et0
psmm3 = PocraSMModelSimulation(
	soil_texture='clayey', soil_depth_category='deep (50 to 100 cm)',
	lulc_type='cropped in two seasons', slope=3,
	rain=[4]*122+[0]*243,
	temp_min=[22]*365, temp_avg=[28]*365, temp_max=[34]*365,
	latitude=18.5,
	crop='bajri'
)

# example 4:
# provide a CSV file containing weather data
# (please check the file used in the example; in particular, note
# that the headers should have exactly same labels as the variable name)
psmm4 = PocraSMModelSimulation(
	soil_texture='clayey', soil_depth_category='deep (50 to 100 cm)',
	lulc_type='cropped in two seasons', slope=3,
	weathers=SimulationIO.create_weathers_from_csv_file('weathers.csv'),
	crop='bajri'
)




####### choose the example to run #######
psmm = psmm4 # psmm1 or psmm2 or ... any other simulation that you may build


# run the chosen example simulation
psmm.run()

# check water-component values after running
# print(f'Primary Runoff is : {psmm.pri_runoff}\n\n')
# print(f'Infiltration is : {psmm.infil}\n\n')
# print(f'AET is : {psmm.aet}\n\n')
# print(f'PET is : {psmm.pet}\n\n')
# print(f'Secondary Runoff is : {psmm.sec_runoff}\n\n')
# print(f'GW-recharge is : {psmm.gw_rech}\n\n')
# # check weather components values, if desired
# print(f'Rain is : {psmm.rain}\n\n')
# print(f'ET0 is : {psmm.et0}\n\n')
# print(f'R_a is : {psmm.r_a}\n\n')

SimulationIO.output_water_components_to_csv(
	psmm,
	components=['pri_runoff', 'avail_sm', 'aet'],
	filepath='results.csv'
)