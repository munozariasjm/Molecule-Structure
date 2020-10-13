import sympy as sy
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j

params_general = {
'mu_B': 1.399624494, #MHz/Gauss
'g_S': 2.0023,
'g_L': 1,
'muE_X': 1.9*0.503412, #Debye in MHz/(V/cm)
'muE_A': 0.43*0.503412,
}

params_174X000 = {
'Be': 7348.40053,
'Gamma_SR': -81.150,
'bF': 4.80,
'c': 2.46,
'b': (4.80-2.46/3)
}

params_174X010 = {
'Be': 7348.40053,
'Gamma_SR': -81.150,
'bF': 4.80,
'c': 2.46,
'b': (4.80-2.46/3),
'q_lD': -10
}

params_174A000 = {
'Be': 7586.3,
'ASO': 4.047*10**7,
'a': 0,         # extrapolated from YbF
'bF': 0.07,     # extrapolated from YbF
'c': -0.18,     # extrapolated from YbF
'p+2q': -13133,
'g_lp': -0.865
}

params_173X000 = { # all units MHz except for muE
'Be': 7353.01,
'Gamma_SR': -81.064,
'bFYb': -1886,
'cYb': -102.0,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3328,
}

params_173X010 = { # all units MHz except for muE
'Be': 7353.01,
'Gamma_SR': -81.064,
'bFYb': -1886,
'cYb': -102.0,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3328,
'q_lD': -10
}

params_173A000 = {
'Be': 7353.01,
'ASO': 4.047*10**7,
'h1/2Yb': -138,
'dYb': -265.3,
'bFH': 0.07,     # extrapolated from YbF
'cH': -0.18,     # extrapolated from YbF
'e2Qq0': -1940,
'p+2q': -13146
}


all_params = {
'174X000':{**params_174X000,**params_general},
'174X010':{**params_174X010,**params_general},
'173X000':{**params_173X000,**params_general},
'173X010':{**params_173X010,**params_general},
'174A000':{**params_174A000,**params_general},
'173A000':{**params_173A000,**params_general}
}
