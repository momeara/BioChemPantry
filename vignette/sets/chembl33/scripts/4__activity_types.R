
# help parse chembl activities
library(plyr)
library(dplyr)

activity_unit_transforms = c(
	# SI molar units
	'fM' = 1e-6,
	'pM' = 1e-3,
	'nM' = 1,
	'uM' = 1e3,
	'mM' = 1e6,
	'M' = 1e9,
	
	# SI moles per ml
	'fmol/ml' = 1e-3,
	'pmol/ml' = 1,
	'nmol/ml' = 1e3,
	'umol/ml' = 1e9,
	'mmol/ml' = 1e12,
	
	# per mole
	'M-1' = 1e-9)
activity_units <- names(activity_unit_transforms)

minus_log <- function(x) 10^(9-x)
plus_log <- function(x) 10^(9+x)
identity <- function(x) x

activity_type_transforms = c(
	# Kd
	'Kb' = identity,
	'pKb' = minus_log,
	'-Log Kb' = minus_log,
	'Log 1/Kb' = minus_log,
	'log(1/Kb)' = minus_log,
	'Log Kb' = plus_log,
	'logKb' = plus_log,

	# Kd
	'Kd' = identity,
	'pKd' = minus_log,
	'-Log Kd' = minus_log,
	'Log 1/Kd' = minus_log,
	'log(1/Kd)' = minus_log,
	'Log Kd' = plus_log,
	'logKd' = plus_log,

	'-Log KD' = minus_log,

	# Ki
	'Ki' = identity,
	'pKi' = minus_log,
	'-Log Ki' = minus_log,
	'Log 1/Ki' = minus_log,
	'log(1/Ki)' = minus_log,
	'Log Ki' = plus_log,
	'logKi' = plus_log,


	# Concentration required to elicit a 50% response in an in vitro assay
	'AC50' = identity,
	'-Log AC50' = minus_log,
	'Log 1/AC50' = minus_log,
	'log(1/AC50)' = minus_log,
	'Log AC50' = plus_log,
	'logAC50' = plus_log,

	# anesthesia dose
	'AD50' = identity,
	'-Log AD50' = minus_log,
	'Log 1/AD50' = minus_log,
	'log(1/AD50)' = minus_log,
	'Log AD50' = plus_log,
	'logAD50' = plus_log,

	 # cytotoxcity
	'CC50' = identity,
	'-Log CC50' = minus_log,
	'Log 1/CC50' = minus_log,
	'log(1/CC50)' = minus_log,
	'Log CC50' = plus_log,
	'logCC50' = plus_log,

	# effective concentration (for antagonists)
	'EC50' = identity,
	'pEC50' = minus_log,
	'-Log EC50' = minus_log,
	'Log 1/EC50' = minus_log,
	'log(1/EC50)' = minus_log,
	'Log EC50' = plus_log,
	'logEC50' = plus_log,

	'EC90' = identity,
	'pEC90' = minus_log,
	'-Log EC90' = minus_log,
	'Log 1/EC90' = minus_log,
	'log(1/EC90)' = minus_log,
	'Log EC90' = plus_log,
	'logEC90' = plus_log,

	# effective dose (for agonists)
	'ED50' = identity,
	'pED50' = minus_log,
	'-Log ED50' = minus_log,
	'Log 1/ED50' = minus_log,
	'log(1/ED50)' = minus_log,
	'Log ED50' = plus_log,
	'logED50' = plus_log,

	# concentration for 50% of maximal inhibition of cell proliferation
	'GI50' = identity,
	'pGI50' = minus_log,
	'-Log GI50' = minus_log,
	'Log 1/GI50' = minus_log,
	'log(1/GI50)' = minus_log,
	'Log GI50' = plus_log,
	'logGI50' = plus_log,

	'IA50' = identity,
	'pIA50' = minus_log,
	'pIA50' = minus_log,
	'-Log IA50' = minus_log,
	'Log 1/IA50' = minus_log,
	'log(1/IA50)' = minus_log,
	'Log IA50' = plus_log,
	'logIA50' = plus_log,

	'IC50' = identity,
	'pIC50' = minus_log,
	'-Log IC50' = minus_log,
	'Log 1/IC50' = minus_log,
	'log(1/IC50)' = minus_log,
	'Log IC50' = plus_log,
	'logIC50' = plus_log,

	'IC60' = identity,
	'pIC60' = minus_log,
	'-Log IC60' = minus_log,
	'Log 1/IC60' = minus_log,
	'log(1/IC60)' = minus_log,
	'Log IC60' = plus_log,
	'logIC60' = plus_log,

	'IC70' = identity,
	'pIC70' = minus_log,
	'-Log IC70' = minus_log,
	'Log 1/IC70' = minus_log,
	'log(1/IC70)' = minus_log,
	'Log IC70' = plus_log,
	'logIC70' = plus_log,

	'IC80' = identity,
	'pIC80' = minus_log,
	'-Log IC80' = minus_log,
	'Log 1/IC80' = minus_log,
	'log(1/IC80)' = minus_log,
	'Log IC80' = plus_log,
	'logIC80' = plus_log,

	'IC90' = identity,
	'pIC90' = minus_log,
	'-Log IC90' = minus_log,
	'Log 1/IC90' = minus_log,
	'log(1/IC90)' = minus_log,
	'Log IC90' = plus_log,
	'logIC90' = plus_log,

	'IC95' = identity,
	'pIC95' = minus_log,
	'-Log IC95' = minus_log,
	'Log 1/IC95' = minus_log,
	'log(1/IC95)' = minus_log,
	'Log IC95' = plus_log,
	'logIC95' = plus_log,

	'IC99' = identity,
	'pIC99' = minus_log,
	'-Log IC99' = minus_log,
	'Log 1/IC99' = minus_log,
	'log(1/IC99)' = minus_log,
	'Log IC99' = plus_log,
	'logIC99' = plus_log,

	'IC100' = identity,
	'pIC100' = minus_log,
	'-Log IC100' = minus_log,
	'Log 1/IC100' = minus_log,
	'log(1/IC100)' = minus_log,
	'Log IC100' = plus_log,
	'logIC100' = plus_log,

	'ID50' = identity,
	'pID50' = minus_log,
	'-Log ID50' = minus_log,
	'Log 1/ID50' = minus_log,
	'log(1/ID50)' = minus_log,
	'Log ID50' = plus_log,
	'logID50' = plus_log,

	 # lethal dose
	'LD50' = identity,
	'-Log LD50' = minus_log,
	'Log 1/LD50' = minus_log,
	'log(1/LD50)' = minus_log,
	'Log LD50' = plus_log,
	'logLD50' = plus_log,

	# maximum inhibitory concentration
	'MIC' = identity,
	'MIC50' = identity,
	'MIC90' = identity)
activity_types <- names(activity_type_transforms)


# there is probably a better way....
normalize_activity_value <- function(type, unit, value) {
	plyr::llply(1:length(type), function(i) {
		tryCatch({
			activity_type_transforms[[type[[i]] ]](value[[i]] ) *
			activity_unit_transforms[[unit[[i]] ]]
		}, error = function(e){
			print(paste0("error on row ", i, "\n", e))
			stop()
		})
	}) |> unlist()
}


