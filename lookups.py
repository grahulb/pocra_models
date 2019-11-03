
DEFAULT_AVAIL_SM = 0
DEFAULT_SOWING_THRESHOLD = 50

########	Lookup Dictionaries Start	########

# Dictionary for generic land use types:
dict_lulc = {
	'forest-forest blank':'scrub',
	'forest-deciduous (dry/moist/thorn)-open': 'scrub',
	'agricultural land-crop land-rabi crop': 'agriculture',
	'forest-scrub forest': 'scrub',	
	'agricultural land-crop land-kharif crop': 'agriculture',
	'agricultural land-fallow-current fallow': 'fallow land',
	'wastelands-scrub land-open scrub': 'wasteland',
	'wastelands-gullied/ravinous land-gullied': 'wasteland',
	'forest-deciduous (dry/moist/thorn)-dense/closed': 'forest',
	'wastelands-scrub land-dense scrub': 'scrub',
	'built up-built up (rural)-built up area (rural)': 'habitation',
	'waterbodies-reservoir/tanks-dry-zaid extent': 'water',
	'waterbodies-reservoir/tanks-dry-rabi extent': 'water',
	'waterbodies-canal/drain-lined': 'water',
	'agricultural land-crop land-zaid crop': 'agriculture',
	'waterbodies-reservoir/tanks-dry-kharif extent': 'water',
	'agricultural land-crop land-two crop area': 'agriculture',
	'built up-built up (urban)-vegetated area': 'habitation',
	'wastelands-barren rocky/stony waste': 'scrub',
	'agricultural land-plantation-agriculture plantation': 'agriculture',
	'agricultural land-crop land-more than two crop': 'agriculture',
	'waterbodies-river/stream-perennial': 'water',
	'built up-built up (urban)-transportation': 'habitation',
	'built up-built up (urban)-recreational':'habitation',
	'built up-built up (urban)-residential': 'habitation',
	'cropped in more than two seasons':'agriculture',
	'cropped in two seasons':'agriculture',
	'rabi':'agriculture',
	'zaid':'agriculture',
	'kharif':'agriculture',
	'agricultural plantation':'agriculture',
	'deciduousdry/ moist/ thorn - dense/ closed':'forest',
	'evergreen/ semi evergreen - dense/ closed':'forest',
	'forest plantation':'forest',
	'tree clad area - dense/ closed':'forest',
	'fallow land':'fallow land',
	'built up - compactcontinuous':'habitation',
	'built up - compact (continuous)': 'habitation',
	'built up - sparsediscontinuous':'habitation',
	'industrial area':'habitation',
	'rural':'habitation',
	'tree clad area - open':'scrub',
	'deciduousdry/ moist/ thorn - open':'scrub',
	'evergreen/ semi evergreen - open':'scrub',
	'scrub':'scrub',
	'ash/ cooling pond/ effluent and other waste':'wasteland',
	'mining - abandoned':'wasteland',
	'mining - active':'wasteland',
	'quarry':'wasteland',
	'barren rocky':'wasteland',
	'gullied/ ravinous land - gullied':'wasteland',
	'scrub land - dense/ closed':'wasteland',
	'scrub land - open':'wasteland',
	'scrub forest':'scrub',
	'vegetated/ open area':'wasteland',
	'reservoir/ tanks - permanent':'water',
	'reservoir/ tanks - seasonal':'water',
	'river - non perennial':'water',
	'river - perennial':'water',
	'canal/ drain':'water',
	'lakes/ ponds - permanent':'water',
	'lakes/ ponds - seasonal':'water',
	'deciduous (dry/ moist/ thorn) - open': 'scrub',
	'deciduous (dry/ moist/ thorn) - dense/ closed' :  'forest',
	'built up - sparse (discontinuous)' : 'habitation'	
}

# Lookup for properties of various soils
soil_properties = [
	'hsg', 'Sand %', 'Clay %', 'Gravel %', 'wp', 'fc',
	'sat', 'ksat', 'Bulk Density',	'AWC'
]
dict_SoilProperties_in_soil_properties_order = {
	'clay loam': ('D', 32, 34, 0, 0.206, 0.341, 0.442, 2.7, 1.48, 0.14),
	'clayey': ('D', 28, 51, 0, 0.303, 0.427, 0.487, 0.52, 1.36, 0.12),
	'gravelly clay': ('D', 23, 48, 10, 0.285, 0.415, 0.488, 0.83, 1.36, 0.12),
	'gravelly clay loam': ('D', 31, 34, 10, 0.206, 0.343, 0.444, 2.32, 1.47, 0.12),
	'gravelly loam': ('B', 41, 17, 10, 0.109, 0.244, 0.408, 10.83, 1.57, 0.12),
	'gravelly sandy clay loam': ('B', 49, 26, 10, 0.16, 0.273, 0.412, 5.83, 1.56, 0.1),
	'gravelly sandy loam': ('B', 63, 10, 10, 0.065, 0.158, 0.402, 33.29, 1.58, 0.08),
	'gravelly silty clay': ('C', 7, 47, 10, 0.277, 0.42, 0.512, 1.7, 1.29, 0.13),
	'gravelly silty loam': ('C', 21, 15, 10, 0.099, 0.282, 0.415, 6.8, 1.55, 0.16),
	'loamy': ('B', 42, 20, 0, 0.126, 0.256, 0.411, 10.2, 1.56, 0.13),
	'loamy sand': ('A', 82, 8, 0, 0.05, 0.106, 0.41, 69.09, 1.56, 0.06),
	'sandy': ('A', 91, 5, 0, 0.03, 0.071, 0.424, 108.06, 1.53, 0.04),
	'sandy clay': ('D', 51, 42, 0, 0.254, 0.364, 0.43, 0.73, 1.51, 0.11),
	'sandy clay loam': ('C', 57, 28, 0, 0.172, 0.271, 0.406, 6.09, 1.57, 0.1),
	'sandy loam': ('A', 65, 11, 0, 0.172, 0.258, 0.399, 6.67, 1.59, 0.09),
	'silty clay': ('D', 9, 46, 0, 0.272, 0.415, 0.506, 1.9, 1.31, 0.14),
	'silty clay loam': ('D', 11, 34, 0, 0.206, 0.371, 0.47, 2.65, 1.41, 0.17),
	'silty loam': ('B', 19, 16, 0, 0.105, 0.291, 0.418, 6.97, 1.54, 0.19),
	'waterbody mask': ('D', 28, 51, 0, 0.303, 0.427, 0.487, 0.52, 1.36, 0.12),
	'habitation mask': ('D', 32, 34, 0, 0.206, 0.341, 0.442, 2.7, 1.48, 0.14)
}
# following dict provides a friendlier API for soil properties
dict_soil_properties = {
	soil_type: dict(zip(soil_properties, soil_property_values))
		for soil_type, soil_property_values in dict_SoilProperties_in_soil_properties_order.items()
}

# Lookup for SCS curve no based on land ussage and HSG: 
dict_lulc_hsg_curveno = {
	'agriculture': {'A': 67, 'B': 78, 'C': 85, 'D': 89},
	'forest': {'A': 30, 'B': 55, 'C': 70, 'D': 77},
 	'fallow land': {'A': 77, 'B': 86, 'C': 91, 'D': 94},
	'habitation': {'A': 77, 'B': 85, 'C': 90, 'D': 92},
	'scrub': {'A': 49, 'B': 69, 'C': 79, 'D': 84},
	'wasteland': {'A': 68, 'B': 79, 'C': 86, 'D': 89},
	'water': {'A': 100, 'B': 100, 'C': 100, 'D': 100},
	'current fallow': {'A': 72, 'B': 81, 'C': 88, 'D': 91}
	}

# Lookup for soil depth with respect to given soil depth category:
dict_soil_depth_category_to_value = {
	'deep (50 to 100 cm)': 1,
	'habitation mask': 0.1,
	'shallow (10 to 25 cm)': 0.25,
	'very deep (> 100 cm)': 1.5,
	'deep to very deep (> 50 cm)': 1,
	'waterbody mask': 0.1,
	'moderately deep (25 to 50 cm)': 0.5,
	'shallow to very shallow (< 25 cm)': 0.25,
	'very shallow (< 10 cm)': 0.1
}

# Lookup for crop properties: KC, depletion factor, root-depth
crop_properties = ['kc', 'depletion_factor', 'root_depth']
dict_crop = {
	'rice':			    ([1.15]*30 + [1.23]*30 + [1.14]*80 + [1.02]*40,	    0.2,	0.75),
	'bajri':		 	([0.34]*13 +  [0.67]*21 + [1.05]*34 + [0.62]* 22, 	0.55, 	1.0),
	'banana':		 	([0.53]*112 + [1.17]*84 + [1.06]*112 + [1.06]*7,	0.35,	0.5),
	'brinjal': 			([0.51]*44 + [0.84]*58 + [1.29]*58 + [0.9]*30, 		0.45,	0.7),
	'cauliflower': 		([0.63]*14 + [1.05]*18 + [1.46]*43 + [1.25]*10, 	0.45, 	0.4),
	'citrus':			([0.7]*60 + [0.65]*90 + [0.7]*120 + [0.7]*95, 		0.5,	1.1),
	'cotton':			([0.51]*30 + [0.85]*50 + [1.3]*55 + [0.85]* 45,		0.65, 	1.0),
	'fodder_crop': 		([0.35]*14 + [0.7]*25 + [1.01]*29 + [0.61]*22, 		0.55,	0.8),
	'grapes': 			([0.44]*30 + [1.52]*61 + [0.73]*183 + [0.73]*91, 	0.35,	1.0),
	'groundnut': 		([0.47]*23 + [0.79]*32 + [1.1]*42 + [0.74]*23, 		0.5,	0.5),
	'maize': 			([0.56]*14 + [1.11]*25 + [1.6]*29 + [0.97]*22,	 	0.55,	0.9),
	'mirchi': 			([0.44]*40 + [0.87]*55 + [1.31]*63 + [1.12]* 32, 	0.3,	0.5),
	'moong': 			([0.57]*8 + [0.95]*12 + [1.4]*24 + [0.63]*16, 		0.4,	0.6),
	'mosambi': 			([0.7]*60 + [0.65]*90 + [0.7]*120 + [0.7]*95, 		0.5,	1.1),
	'onion': 			([0.53]*12 + [0.75]*19 + [1.07]*54 + [1.07]*30, 	0.35,	0.3),
	'orange': 			([0.7]*60 + [0.65]*90 + [0.7]*120 + [0.7]*95, 		0.5,	1.1),
	'pomegranate': 		([0.46]*21 + [0.26]*77 + [0.56]*56 + [0.7]*211, 	0.5,	1.1),
	'potato': 			([0.62]*25 + [1.03]*30 + [1.58]*30 + [1.16]*20, 	0.35,	0.4),
	'small_vegetables': ([0.73]*20 + [0.97]*20 + [1.61]*15 + [1.45]*5, 		0.3,	0.3),
	'sorghum':			([0.34]*20 + [0.72]*30 + [1.06]*40 + [0.63]*30, 	0.55,	1.0),
	'soyabean': 		([0.33]*16 + [0.7]*23 + [1.03]*47 + [0.56]*19, 		0.5,	0.6),
	'sugarcane': 		([0.51]*28 + [1.58]*48 + [0.95]*151 + [0.95]*138, 	0.5,	1.2),
	'sunflower': 		([0.36]*17 + [0.78]*29 + [1.19]*38 + [0.57]*21, 	0.45,	0.8),
	'tomato': 			([0.58]* 30 + [0.97]*41 + [1.49]*41 + [1.04]*25, 	0.4,	0.7),
	'tur': 				([0.43]*28 + [0.72]*46 + [1.1]*50 + [0.72]*41, 		0.65,	1.0),
	'turmeric': 		([0.59]*40 + [0.98]*67 + [1.51]*73 + [0.98]*60, 	0.65,	1.0),
	'udid': 			([0.41]*11 + [0.69]*17 + [1.01]*33 + [0.46]*22, 	0.4,	0.6),
	'vegetables': 		([0.53]*24 + [0.89]*33 + [1.36]*33 + [0.94]*20, 	0.35,	0.4)
 }
# Lookup for proxy values for crop-properties 
# for vegetation on various LULC types and for other situations:.
# Various types of 'scrub' have been pooled together as 'scrub'.
dict_lulc_pseudo_crop =	{
	'forest':					(	[0.3]*45 + [1.15]*60 + [0.7]*90 + [0.1]*170,	0.8,	3	),
	'wasteland':				(	[0.5]*120 + [0.25]*60 + [0.15]*120 + [0.1]*65,	0.5,	0.5	),
	'scrub':					(	[0.3]*30 + [0.7]*60 + [0.5]*60 + [0.2]*215,		0.6,	1.5	),
	'current fallow crop':		(	[0.2]*60 + [0.3]*62 ,		0.5,	0.9		),
	'permanent fallow crop':	(	[0.2]*60 + [0.3]*62 ,		0.5,	0.9		),
	'non agri':					(	[0.2]*60 + [0.3]*62 ,		0.5,	0.9		),
}

# following dict provides a friendlier API for properties of crops and crop-likes
dict_of_properties_for_crop_and_croplike = {
	crop: dict(zip(crop_properties+['is_pseudo_crop'], crop_property_values+(False,)))
		for crop, crop_property_values in dict_crop.items()
}
dict_of_properties_for_crop_and_croplike.update({
	crop_like: dict(zip(crop_properties+['is_pseudo_crop'], crop_property_values+(True,)))
		for crop_like, crop_property_values in dict_lulc_pseudo_crop.items()
})

# Rabi crop with their fixed PET values
dict_rabi_crop ={
	'rabi_sorghum' : 425,
	'gram' : 250,
	'rabi_maize' : 500,
	'rabi_onion' : 350,
	'rabi_mirchi' : 800,
	'rabi_tomato' : 400,
	'rabi_brinjal' : 650,
	'rabi_vegetables' : 350,
	'rabi_cauliflower' : 350,
	'rabi_potato' : 500,
	'rabi_okra' : 550,
	'rabi_groundnut' : 650,
	'rabi_sunflower' : 400,
	'rabi_fodder' : 350,
	'rabi_wheat' : 500
	}

# Long Kharif crops list
long_kharif_crops = ['cotton','orange','grapes','tur']

########	Lookup Dictionaries End		########