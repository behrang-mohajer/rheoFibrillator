# SB_droplet_distribution
This python library simulate droplets distribution within a steady-state flow. It inputs
 - ```input_data.py``` for process condition,
 - an OpenFOAM file ```XXX.foam```, and
 - the initial particle relase over a volume or a patch that can repeat or alternate over time.

It requires these open-source packages to run: PyVitsa, Numpy, jupyet notebook (or similar), and optionally Numexpr for parallalization.
