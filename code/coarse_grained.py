import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random

# constants
diff_coeff = 1.49e-7  # diffusion coefficient
ActivationEnergy_Cu = 134.5e3
atomic_distance = 2.54e-10 # = 2.54 Angstrom
sample_thickness = 20e-6 # = 20 micrometers
al_concentration = 0.1

# parameters
N = 100
L = sample_thickness/N

class Box():
    def __init__(self, sites, al_atoms):
        self.sites = sites
        self.al_atoms = al_atoms
        self.cu_atoms = sites - al_atoms
        self.al_concentration = al_atoms/sites
        self.above = None
        self.below = None