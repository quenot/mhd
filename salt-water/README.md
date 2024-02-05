# mhd / salt-water
The model.py file contains a number of fuctions for linking the salinity (S), the temperature (t), the Molarity (M), the volumic mass (rho), the conductivity (sigma), the dynamic viscosity (mu) and the kinematic viscosity (nu) in salt water.

When S and t are given, rho, mu, nu and sigma can be obtained by:

import models <br>
rho = models.DensityFromSalinityTemperature(S=S, t=t) <br>
mu = models.DynamicViscosityFromSalinityTemperature(S=S, t=t) <br>
nu = mu/rho <br>
sigma = models.ConductivityFromSalinityTemperature(S=S, t=t)

See domains of validity in the code.
