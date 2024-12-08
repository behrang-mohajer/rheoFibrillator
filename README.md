# \texttt{rheoFibrillator}
This package was introduced to simulate droplet fibrillation in the spunbonding process.
Many features are still under development, but fill free to run Allrun script and give it a try.

This Python library simulates spherical droplets distribution and deformation within a steady-state flow. It inputs
 - ```input_data.py``` for process condition,
 - an OpenFOAM (```XXX.foam```) or a VTK (```YYY.vtk```), and
 - the initial particle release over a volume or a patch that can repeat or alternate over time.

It requires these open-source packages to run: PyVitsa, Numpy, Jupyter Notebook (or similar), and optionally Numba for parallelization.

The users pick where and how often droplets are released in a random or designated initial distribution (mass and location). Then they pick the direction of droplet elongation, and this package simulates what happens in terms of:
-elongation into ellipsoids and fibrils
-distribution (radial forces are under development)
-breakup 
-coalescence
