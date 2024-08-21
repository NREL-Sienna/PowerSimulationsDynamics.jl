"""
Validation PSSE/PIDGOV:
This case study defines a three bus system with an infinite bus, GENROU+SEXS+PIDGOV and a load.
The fault drop the line connecting the infinite bus and GENROU
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

raw_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PIDGOV/ThreeBusMulti.raw")
dyr_file = joinpath(TEST_FILES_DIR, "benchmarks/psse/PIDGOV/ThreeBus_PIDGOV_PowerFlag.dyr")
