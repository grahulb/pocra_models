"""
This module contains all the components required to run
PoCRA's soil-moisture model. This model is a point-based simulation model,
in the sense that it represents the interplay, over time, of soil-moisture
with other water balance components at a given location (point, in GIS terms).
The model is based on four broad bio-physical entities:
1. Field: The field conditions at the location
2. Crop: Crop sown at the location, if any
3. Weather: The weather conditions at the location
4. Water: The water-components conceptualized by the soil-moisture model.
The module is structured in-line with this to have four corresponding classes.

Note: The processing of equations for model-simulation has been implemented
as static methods to allow their general-purpose (black-box) use not tied to
the corresponding class' structure, and yet allowing them
to be kept within the class' premises for clear conceptual design
that indicates the entity(modelled by its class) which that method simulates.

PoCRA's soil-moisture model has been customized from SWAT's models
and thus contains sub-models within it. Each of these sub-models generally
belongs to one of the four aforementioned bio-physical entities and thus
have been coded into their corresponding class. In particular,
1.	<Field> class includes the 'model' that determines the set-up of
	field properties with respect to soil-moisture.
2.	<Weather> class includes two models. One computes solar-radiation
	for a given latitude and a given day-of-year and the other computes
	reference evapotranspiration.
3. <Water> class includes two models. One computes potential evapotranspiration
	and the other computes all the water-components for the day(time-step).

Although all these models have been implemented as functions within the class,
they have been designated as static-methods and implemented such that their
results depend only on their input parameters and nothing else.
One major implication of this form of implementation is that anyone with
basic skills of writing arithmetic statements in python can modify the 
internals of the models easily without affecting other parts of the code.

(The user is assumed to know the models and meanings of the input-parameters.
Their documentation is available on IIT's PoCRA webpage.)

Finally, the <PocraSMModelSimulation> class provides an API to simulate
PoCRA's soil-moisture model over a duration of days(time-steps).
The siulation is done for as many days(time-steps) as the length of arrays
of the weather parameters (like rain) given as input to the class constructor.
See this class' documentation to know its API.
"""

import os
import csv
import math

import lookups



class Field:
	"""
	Represents the field conditions at the modelled location.
	These field conditions currently include static properties like:
	1. Soil properties
	2. Land-Use-Land-Cover
	3. Terrain slope
	which are assumed to be static for the purpose of
	PoCRA's soil-moisture model.
	"""
	
	def __init__(s, soil_texture, soil_depth_category, lulc_type, slope):
		"""
		Set basic field properties
		Set derived field parameters required in the model
		"""

		s.soil_texture = soil_texture.lower()
		s.soil_depth_category = soil_depth_category.lower()
		s.lulc_type = lulc_type.lower()
		s.slope = slope


		#### Set parameters actually required in soil-moisture model. ####
		
		# Derived from lookups
		soil_texture_properties = lookups.dict_soil_properties[s.soil_texture]
		s.wp = soil_texture_properties['wp']
		s.fc = soil_texture_properties['fc']
		s.sat = soil_texture_properties['sat']
		s.ksat = soil_texture_properties['ksat']
		
		s.cn_val = lookups.dict_lulc_hsg_curveno[
			lookups.dict_lulc[s.lulc_type]
		][soil_texture_properties['hsg']]
		
		s.soil_depth = lookups.dict_soil_depth_category_to_value[s.soil_depth_category]

		# Derived by calculation
		field_setup = Field.pocra_sm_model_field_setup(
			s.wp, s.fc, s.sat, s.soil_depth, s.cn_val, s.slope, s.ksat
		)
		s.smax = field_setup['smax']
		s.w1 = field_setup['w1']
		s.w2 = field_setup['w2']
		s.daily_perc_factor = field_setup['daily_perc_factor']


	@staticmethod
	def pocra_sm_model_field_setup(wp, fc, sat, soil_depth, cn_val, slope, ksat):
		
		# some utility variables
		sat_minus_wp_depth = (sat-wp) * soil_depth * 1000
		fc_minus_wp_depth = (fc-wp) * soil_depth * 1000
		sat_minus_fc_depth = (sat-fc) * soil_depth * 1000

		# smax
		cn3 = cn_val * math.exp( 0.00673 * (100-cn_val) )
		if (slope > 5.0):
			cn_val = ( (
				((cn3 - cn_val) / 3) * ( 1 - 2 * math.exp(-13.86*slope*0.01) )
			) + cn_val )
		cn1_s = ( cn_val - 
			20 * (100-cn_val) / ( 100-cn_val + math.exp(2.533 - 0.0636*(100-cn_val)) )
		)
		cn3_s = cn_val * math.exp(0.00673*(100-cn_val))
		smax = 25.4 * (1000/cn1_s - 10)
		if smax == 0:
			raise Exception('smax is zero')
		
		# w2
		s3 = 25.4 * (1000/cn3_s - 10)
		w2 = ((
			math.log(fc_minus_wp_depth/(1-s3/smax) - fc_minus_wp_depth)
			- math.log (sat_minus_wp_depth/(1-2.54/smax) - sat_minus_wp_depth)
		) / (sat_minus_fc_depth) )

		# w1
		w1 = (
			math.log(fc_minus_wp_depth/(1- s3/smax) - fc_minus_wp_depth)
			+ w2 * fc_minus_wp_depth
		)

		# daily_perc_factor
		TT_perc = sat_minus_fc_depth/ksat
		daily_perc_factor = 1 - math.exp(-24 / TT_perc)

		return {
			'smax': smax,
			'w1': w1,
			'w2': w2,
			'daily_perc_factor': daily_perc_factor
		}



class Crop:
	"""
	Represents those properties of the crop or other vegetation at the
	modelled location, that play a role in simulating the model.
	"""
	
	def __init__(self, name):

		self.name = name

		if self.name in lookups.dict_of_properties_for_crop_and_croplike:
			crop_properties = lookups.dict_of_properties_for_crop_and_croplike[self.name]
		self.kc = crop_properties['kc']
		self.depletion_factor = crop_properties['depletion_factor']
		self.root_depth = crop_properties['root_depth']
		self.is_pseudo_crop = crop_properties['is_pseudo_crop']



class Weather:
	"""
	Represents the weather conditions at the modelled location
	that play a role in simulating the model.
	"""
	
	G_sc = 0.082 # used in estimating R_a(radiation)

	def __init__(self,
		rain, et0=None,
		temp_min=None, temp_avg=None, temp_max=None, latitude=None,
		r_a=None
	):
		
		self.rain = rain
		self.et0 = et0
		self.temp_min = temp_min
		self.temp_avg = temp_avg
		self.temp_max = temp_max
		self.latitude = latitude
		self.r_a = r_a


	@staticmethod
	def get_pocra_radiation_for_latitude_for_day(latitude, day_of_year):
	
		phi = latitude * math.pi/180

		d_r = 1 + 0.033 * math.cos(2*math.pi * day_of_year/365)
		delta = 0.409 * math.sin((2*math.pi * day_of_year/365) - 1.39)
		omega_s = math.acos(-math.tan(phi) * math.tan(delta))
		r_a = (24*60/math.pi) * Weather.G_sc * d_r * (
			omega_s*math.sin(phi)*math.sin(delta)
			+ math.cos(phi)*math.sin(omega_s)*math.cos(delta)
		)

		return r_a


	@staticmethod
	def get_pocra_et0_for_day_from_weather(
		temp_min, temp_avg, temp_max, r_a=None, latitude=None, day_of_year=None
	):
		
		r_a = r_a or Weather.get_pocra_radiation_for_latitude_for_day(latitude, day_of_year)
		et0 = 0.0023 * (temp_avg + 17.28) * ((temp_max-temp_min)**0.5) * r_a * 0.408

		return r_a, et0



class Water:
	"""
	This represents the components of water-balance
	over the duration of model simulation.
	This class also provides the function <run_pocra_sm_model_for_day>
	that simulates PoCRA's soil-moisture model for a single daily time-step.
	"""

	def __init__(self,
		pri_runoff=None, infil=None, aet=None, sec_runoff=None,
		gw_rech=None, avail_sm=None, pet=None
	):

		self.pri_runoff = pri_runoff
		self.infil = infil
		self.aet = aet
		self.sec_runoff = sec_runoff
		self.gw_rech = gw_rech
		self.avail_sm = avail_sm
		self.pet = pet


	@staticmethod
	def get_pocra_pet_for_day(
		kc, et0=None, r_a=None, latitude=None, day_of_year=None,
		temp_min=None, temp_avg=None, temp_max=None
	):

		if et0 is not None:
			pet = kc * et0
		else:
			r_a, et0 = Weather.get_pocra_et0_for_day_from_weather(
				temp_min, temp_avg, temp_max, r_a, latitude, day_of_year
			)
			pet = kc * et0

		return r_a, et0, pet

	
	@staticmethod
	def run_pocra_sm_model_for_day(
		layer_1_thickness, layer_2_thickness, # layer-dimensions
		prev_avail_sm, sm1_frac, sm2_frac, # soil-moisture state at day-start
		wp, fc, sat, smax, w1, w2, daily_perc_factor, # soil-properties
		depletion_factor, # parameter determined only by crop
		rain, # parameter determined only by weather
		pet # parameter determined by weather and crop
	):
		l1 = layer_1_thickness
		l2 = layer_2_thickness

		####### pri_runoff #######
		s_swat = smax * ( 1 - 
			prev_avail_sm / ( prev_avail_sm + math.exp(w1 - w2 * prev_avail_sm) )
		)
		ia_swat = 0.2 * s_swat
		if rain <= ia_swat:
			pri_runoff = 0
		else:
			pri_runoff = ((rain - ia_swat)**2 ) / (rain + 0.8*s_swat)
		
		####### infil #######
		infil = rain - pri_runoff
		
		####### aet #######
		if (sm1_frac < wp):
			ks = 0
		elif ( sm1_frac > (fc * (1-depletion_factor) + depletion_factor * wp) ):
			ks = 1
		else:
			ks = (sm1_frac - wp) / (fc - wp) / (1-depletion_factor)
		aet = ks * pet
		
		# sm1_before r_to_second_layer(in metres) and sm2_before gw_rech
		sm1_before = ((sm1_frac * l1) + ((infil - aet) / 1000)) / l1
		if (sm1_before < fc):
			r_to_second_layer = 0
		elif (sm2_frac < sat):
			r_to_second_layer = min(
				(sat - sm2_frac) * l2,
				(sm1_before - fc) * l1 * daily_perc_factor
			)
		else:
			r_to_second_layer = 0
		sm2_before = (sm2_frac * l2 + r_to_second_layer) / l2
		
		####### sec_runoff #######
		candidate_new_sm1_frac = (sm1_before * l1 - r_to_second_layer) / l1
		candidate_sec_runoff = ( candidate_new_sm1_frac - sat ) * l1 * 1000
		sec_runoff = max(candidate_sec_runoff, 0)
		new_sm1_frac = min(candidate_new_sm1_frac, sat)
		
		####### gw_rech #######
		candidate_gw_rech = (sm2_before - fc) * l2 * daily_perc_factor * 1000
		gw_rech = max(candidate_gw_rech, 0)
		candidate_sm2_frac = (sm2_before * l2 - gw_rech / 1000) / l2
		new_sm2_frac = min(candidate_sm2_frac, sat)

		####### avail_sm #######
		avail_sm = (new_sm1_frac * l1 + new_sm2_frac * l2 - wp * (l1+l2)) * 1000


		return (
			Water(pri_runoff, infil, aet, sec_runoff, gw_rech, avail_sm, pet),
			{'avail_sm': avail_sm, 'sm1_frac': new_sm1_frac, 'sm2_frac': new_sm2_frac},
		)



class PocraSMModelSimulation:
	"""
	This represents the PoCRA's SM Model for a particular location
	and its simulation over days(time-steps).
	
	Usage:
	>>> from pocra_sm_model import *
	>>> psmm = PocraSMModelSimulation(<various input-parameters>)
	>>> psmm.run()
	>>> aet_list = psmm.aet

	In general after <run>ning the simulation, the list (indexed by time-step)
	of any of the following water-component can be obtained
	(like aet's list was obtained above): avail_sm, pri_runoff, infil,
	aet, pet, sec_runoff and gw_rech.
	"""

	def __init__(self,
		# field-related attributes
		soil_texture=None, soil_depth_category=None, lulc_type=None, slope=None, field=None,
		# weather-related attributes
		weathers=None, rain=None, et0=None,
		r_a=None, temp_min=None, temp_avg=None, temp_max=None, latitude=None,
		# crop
		crop=None,
		# attribute determined by crop+weather
		pet=None,
		# attributes setting the soil-moisture state at the location
		model_state_at_start=None, sowing_date_offset=None, sowing_threshold=None
	):
		"""
		TODO: update this __doc__ as per the new code
		Set model entities.
		field, crop, weathers and waters are instances(plurals are <lists> of instances)
		of corresponding classes.
		state should be a <dict> with three model properties:
		1. key 'avail_sm' : available soil-moisture at the beginning of simulation
		2. key 'sm1_frac' : soil-moisture content in layer 1 expressed as a fraction
		3. key 'sm2_frac' : soil-moisture content in layer 2 expressed as a fraction
		"""

		self.field = field or Field(
			soil_texture, soil_depth_category, lulc_type, slope
		)

		if weathers is not None:
			if all(isinstance(weather, Weather) for weather in weathers):
				# <list>-of-<Weather>s option is an API
				# for the typical case of re-use of same weather conditions
				# for multiple locations
				self.weathers = weathers
		else:
			self.weathers = [Weather(rain=r) for r in rain]
			# print(temp_min)
			if et0 is not None:
				for i in range(len(self.weathers)):
					self.weathers[i].et0 = et0[i]
			else:
				for i in range(len(self.weathers)):
					self.weathers[i].temp_min = temp_min[i]
					self.weathers[i].temp_avg = temp_avg[i]
					self.weathers[i].temp_max = temp_max[i]
				if r_a is not None:
					for i in range(len(self.weathers)):
						self.weathers[i].r_a = r_a[i]
				else:
					for i in range(len(self.weathers)):
						self.weathers[i].latitude = latitude
		self.simulation_length = len(self.weathers)

		self.crop = Crop(crop) if isinstance(crop, str) else crop
		self.pet = pet

		self.model_state = model_state_at_start or {
			'avail_sm': lookups.DEFAULT_AVAIL_SM,
			'sm1_frac': self.field.wp, 'sm2_frac': self.field.wp
		}

		self.sowing_date_offset = sowing_date_offset
		self.sowing_threshold = sowing_threshold or lookups.DEFAULT_SOWING_THRESHOLD

			
		self.waters = [None]*self.simulation_length


	def __getattr__(self, name):
		if name in [
			'pri_runoff', 'infil', 'aet', 'sec_runoff', 'gw_rech', 'avail_sm', 'pet'
		]:
			return [getattr(w, name) for w in self.waters]
		if name in ['rain', 'et0'] or (
			name in ['temp_min', 'temp_avg', 'temp_max', 'r_a'] and hasattr(self.weathers[0], name)
		):
			return [getattr(w, name) for w in self.weathers]
	

	def computation_before_iteration(self):
		
		# determine layer_1_thickness, layer_2_thickness
		if (self.field.soil_depth <= self.crop.root_depth): # thin soil layer
			self.layer_1_thickness = self.field.soil_depth - 0.05
			self.layer_2_thickness = 0.05
		else:
			self.layer_1_thickness = self.crop.root_depth
			self.layer_2_thickness = self.field.soil_depth - self.crop.root_depth

		if self.pet is None and self.sowing_date_offset is None:
			if self.crop.is_pseudo_crop:
				self.sowing_date_offset = 0
			else:
				# determine sowing_date_offset based on sowing_threshold logic
				accumulated_rain = 0
				for i in range(365):
					accumulated_rain += self.weathers[i].rain
					if accumulated_rain >= self.sowing_threshold:
						self.sowing_date_offset = i
						break


	def iterate(s):
		f = s.field
		for i in range(len(s.weathers)):
			if s.pet is None:
				if s.sowing_date_offset <= i < (s.sowing_date_offset + len(s.crop.kc)):
					kc = s.crop.kc[i - s.sowing_date_offset]
				else:
					kc = 0
				if s.weathers[i].et0 is not None:
					_, _, pet_for_day = Water.get_pocra_pet_for_day(kc, et0=s.weathers[i].et0)
				elif None not in [
					s.weathers[i].temp_min, s.weathers[i].temp_avg, s.weathers[i].temp_max
				]:
					params_for_pet_for_day = {
						'temp_min': s.weathers[i].temp_min,
						'temp_avg': s.weathers[i].temp_avg,
						'temp_max': s.weathers[i].temp_max,
					}
					if s.weathers[i].r_a is not None:
						params_for_pet_for_day['r_a'] = s.weathers[i].r_a
					elif s.weathers[i].latitude is not None:
						params_for_pet_for_day['latitude'] = s.weathers[i].latitude
						params_for_pet_for_day['day_of_year'] = ((i+1)+151) if (i+1) <= 214 else (215-i) # leap-years may be tackled if year can be provided by user
					R_a_for_day, et0_for_day, pet_for_day = Water.get_pocra_pet_for_day(kc, **params_for_pet_for_day)
					s.weathers[i].r_a = R_a_for_day
					s.weathers[i].et0 = et0_for_day
				else:
					print(i, s.weathers[i].__dict__)

			else:
				pet_for_day = s.pet[i]

			s.waters[i], s.model_state = Water.run_pocra_sm_model_for_day(
				s.layer_1_thickness, s.layer_2_thickness,
				s.model_state['avail_sm'], s.model_state['sm1_frac'], s.model_state['sm2_frac'],
				f.wp, f.fc, f.sat, f.smax, f.w1, f.w2, f.daily_perc_factor,
				s.crop.depletion_factor,
				s.weathers[i].rain, pet_for_day
			)

		s.pet = [w.pet for w in s.waters]
	
	
	def computation_after_iteration(self):
		
		self.crop_end_index = min(self.sowing_date_offset + len(self.crop.kc), 364)



	def run(self):
		self.computation_before_iteration()
		self.iterate()
		self.computation_after_iteration()



class SimulationIO:
	"""
	This class provides input-output facilities
	to build a <PocraSMModelSimulation> instance.
	"""
	
	@staticmethod
	def create_weathers_from_csv_file(filepath):

		if os.path.exists(filepath):
			with open(filepath, newline='') as f:
				return [
					Weather(**{k: float(v) for k, v in row.items()}) for row in csv.DictReader(f)
				]
	

	@staticmethod
	def output_water_components_to_csv(psmm, components=[], filepath='results.csv'):

		if (isinstance(psmm, PocraSMModelSimulation)
			and len(psmm.waters) > 0 and all(isinstance(w, Water) for w in psmm.waters)
		):
			rows = [{
				c: (getattr(psmm.waters[i], c) if hasattr(psmm.waters[i], c)
					else getattr(psmm.weathers[i], c) if hasattr(psmm.weathers[i], c)
					else 'No such component found')
				for c in components
			} for i in range(len(psmm.weathers))]
			
			with open(filepath, 'w', newline='') as f:
				writer = csv.DictWriter(f, fieldnames=components)
				writer.writeheader()
				writer.writerows(rows)