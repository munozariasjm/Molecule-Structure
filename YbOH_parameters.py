import sympy as sy
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j

params_general = {
'mu_B': 1.399624494, #MHz/Gauss
'g_S': 2.0023,
'g_L': 1,
'muE_X': 1.09*0.503412, #Debye in MHz/(V/cm)
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
'Gamma_Prime': 0,
'bF': 4.80,
'c': 2.46,
'b': (4.80-2.46/3),
'q_lD': 13
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
'Be': 7351.24,
'Gamma_SR': -81.06,
'bFYb': -1883.21,
'cYb': -81.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3318.70,
}

params_171X000 = {
'Be': 7359.81,
'Gamma_SR': -80.85,
'bFYb': 6823.58,
'cYb': 233.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': 0,
}

params_173X010 = { # all units MHz except for muE
'Be': 7351.2,
'Gamma_SR': -81.064,
'Gamma_Prime':0,
'bFYb': -1883.21,
'cYb': -81.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3318.7,
'q_lD': -10
}

params_173A000 = {
'Be': 7590.30,
'ASO': 4.047*10**7,
'h1/2Yb': -126.51,
'dYb': -261.72,
'bFH': 0.07,     # extrapolated from YbF
'cH': -0.18,     # extrapolated from YbF
'e2Qq0': -1924.66,
'p+2q': -13141.63
}

params_171A000 = {
'Be': 7597.79,
'ASO': 4.047*10**7,
'h1/2Yb': 443.69,
'dYb': 959.04,
'bFH': 0.07,     # extrapolated from YbF
'cH': -0.18,     # extrapolated from YbF
'e2Qq0': 0,
'p+2q': -13150.91
}

all_params = {
'174X000':{**params_174X000,**params_general},
'174X010':{**params_174X010,**params_general},
'173X000':{**params_173X000,**params_general},
'173X010':{**params_173X010,**params_general},
'174A000':{**params_174A000,**params_general},
'173A000':{**params_173A000,**params_general},
'171A000':{**params_171A000,**params_general},
'171X000':{**params_171X000,**params_general}
}
