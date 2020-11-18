"""
Validation PSSE/GAST:
This case study defines a three bus system with an infinite bus, GENROU and a load.
The GENROU machine has connected an ESAC1A Excitation System.
The fault drop the line connecting the infinite bus and GENROU.
"""

##################################################
############### SOLVE PROBLEM ####################
##################################################

include(joinpath(dirname(@__FILE__), "data_tests/data_14_buses_GAST.jl"))
