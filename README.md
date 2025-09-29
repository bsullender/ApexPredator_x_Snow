Ben Sullender
Sept 2025

This code and input data accompany the manuscript "Apex predators exploit advantagenous snow conditions across hunting modes," currently in review at Journal of Animal Ecology.

Input data are anonymized, without lat/long coordinates, but can still be used with the code to reproduce our analyses. Please refer to our manuscript for further details.

File descriptions:
1) cougarKills.csv - anonymized cougar kill site location data frame, 2533 rows x 9 cols. For every actual cougar kill (case = 1), we generated 5 random locations (case = 0) for comparison.
2) cougarSteps.rds - anonymized cougar location data frame. 
3) wolfSteps.rds - anonymized wolf location data frame.
4) Final_cougarKills.R - code to reproduce cougar kill site selection analysis.
5) Final_predatorSSFs.R - code to reproduce cougar and wolf movement analysis.

For 1)-3), column names are as follows:
* ID_wyr - individual ID conconatenated with water-year
* date - date (m/d/y)
* case - data represents actual observed/used location (case = 1) or randomly generated available location (case = 0)
* snod - snow depth (m)
* dens - snow density (g/m3)
* ShrX - percent shrub within X cells of location. X defined based on median step length.
* deer - deer RSF score, from 0 (lowest value) to 1 (highest value)
* tri - terrain ruggedness index
* cc - canopy coverage, as percentage

For 2) and 3), additional columns are:
* region - NE = northeast / OK = Okangonan
* wyr = water year
* step_id_ - step ID, used for conditional regression
* slo - slope in degrees
* dem - elevation (m)
* ForX - percent forest within X cells of location. X defined based on median step length.
* OpenX - percent open habitat within X cells of location. X defined based on median step length.
* sl - step length (m)
* ta - turn angle (radians)
* dummy - dummy variable (constant, value = 1), used in Cox Proportional Hazards model 
* indiv_step_f - individual step ID, used for conditional regression
