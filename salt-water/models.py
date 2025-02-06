import numpy as np
from scipy.optimize import brentq

t0, S0, T0, M0, rho0 = 15, 35, 273.15, 1.0, 1025.878259
MNaCl = 58.4428e-3 # molar mass of NaCl in kg/L

################################################################################################
# Density rho from salinity S and temperature t at 1 atm
# Adapted from H.T. El-Dessouky and H.M. Ettouney.
# Fundamentals of Salt Water Desalination.
# Chemical, Petrochemical & Process. Elsevier, 2002.
# 0 ‰ < S < 160 ‰ and 10 < t < 180 °C

Crho = np.array([[1008.05475, 57.6565,  0.163],
                 [ -54.0995,   1.571,  -0.423],
                 [  -6.1235,   1.74,   -0.009],
                 [    0.346,  -0.087,  -0.053]])

def DensityFromSalinityTemperature(S=S0, t=t0, target=0.0):
    return np.polynomial.chebyshev.chebval2d((t-100)/80, (S-75)/75, Crho)-target   # kg/m³

def SalinityFromDensityTemperature(rho=rho0, t=t0):
    return brentq(DensityFromSalinityTemperature, 0, 360.0, (t, rho))   # ‰

# Molarity M as a function of Salinity S and vice-versa at temperature t

def MolarityFromSalinityTemperature(S=S0, t=t0, target=0.0):
    return 1e-6*S*DensityFromSalinityTemperature(S=S, t=t)/MNaCl-target # Mol/kg

def SalinityFromMolarityTemperature(M=M0, t=t0):
    return brentq(MolarityFromSalinityTemperature, 0, 360, (t, M))   # ‰


################################################################################################
# Dynamic viscosity from salinity and temperature at 1 atm only
# Adapted from H.T. El-Dessouky and H.M. Ettouney.
# Fundamentals of Salt Water Desalination.
# Chemical, Petrochemical & Process. Elsevier, 2002.
# 0 ‰ < S < 130 ‰ and 10 < t < 130 °C

def DynamicViscosityFromSalinityTemperature(S=S0, t=t0):
    A = 1.474e-3 + 1.5e-5*t - 3.927e-8*t**2
    B = 1.0734e-5 - 8.5e-8*t + 2.23e-10*t**2
    return((1 + A*S + B*S**2)*np.exp(-3.79418 + 604.129/(139.18+t))/1000)  # kg/m.s or Pa.s


################################################################################################
# Conducivity from salinity and temperature or from molarity and temperature
# Model from data as a (Type 1) Chebyshev polynomial
# Dada from Marija Beˇster-Rogaˇc, R. Neueder, and J. Barthel.
# Conductivity of sodium chloride in water +
# 1,4-dioxane mixtures from 5 to 35°c ii.
# concentrated solution. Journal of Solution Chemistry, 29:51–61, 01 2000.

CLambda = np.array([[ 74.172862, -30.060878, 0.178265, -1.720127, 0.755938],
                    [ 23.480611, -10.081103, 0.838633, -0.412552, 0.272767],
                    [  0.588656,  -0.299116, 0.066510,  0.010868, 0.073651]])

tm, dt = 20.00000, 15.00000
rm, dr = 1.099268, 0.994007

def MolarConductivityFromSalinityTemperature(S=S0, t=t0):
    m = S/(MNaCl*1e3)  # Molality
    return np.polynomial.chebyshev.chebval2d((t-tm)/dt, (np.sqrt(m)-rm)/dr, CLambda)

def ConductivityFromSalinityTemperature(S=S0, t=t0):
    M = S/(MNaCl*1e3)*DensityFromSalinityTemperature(S=S, t=t)/1000
    return MolarConductivityFromSalinityTemperature(S=S, t=t)*M/10

def ConductivityFromMolarityTemperature(M=M0, t=t0):
    return ConductivityFromSalinityTemperature(S=SalinityFromMolarityTemperature(M=M, t=t), t=t)

def ConductivityFromDensityTemperature(rho=rho0, t=t0):
    return ConductivityFromSalinityTemperature(S=SalinityFromDensityTemperature(rho=rho, t=t), t=t)
