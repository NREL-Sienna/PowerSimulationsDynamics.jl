import andes
from andes.utils.paths import get_case

case_path = "11BUS_KUNDUR.raw"
dyr_path =  "11BUS_KUNDUR_TGOV.dyr"
ss = andes.run(case_path, addfile = dyr_path, routine='eig')

import numpy as np

eigs = ss.EIG.mu
eigs_sorted = np.sort_complex(eigs)
np.savetxt("eigs_tgov_andes.csv", eigs_sorted, delimiter = ",")
# An additional post process was done to transform Python complex to Julia complex numbers
