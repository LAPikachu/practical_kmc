import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random

# constants
diff_coeff = 1.49e-7  # diffusion coefficient
ActivationEnergy_Cu = 134.5e3
atomic_distance = 2.54e-10 # = 2.54 Angstrom
sample_thickness = atomic_distance*78740 # \approx 20 micrometers
sample_width = atomic_distance*19.685 # 1/4 sample_thickness
al_concentration = 0.1

# parameters
N = 10
L = sample_thickness/N
atoms_per_box = int((L/atomic_distance) * (sample_width/atomic_distance))

# time parameters
total_time = 10e6

def make_array(al_concentration, atoms_per_box, N):
    array = np.zeros(N, dtype=object)
    array[:N//2] = int(al_concentration*atoms_per_box)
    for i, entry in enumerate(array):
        al_atoms = entry
        cu_atoms = atoms_per_box - al_atoms
        array[i] = {'Al': al_atoms, 'Cu': cu_atoms, 'possible_jumps': None}
    for i, cell in enumerate(array):
        if i == 0:
            cell['possible_jumps'] = min(cell['Al'], array[i+1]['Cu'])
        elif i == N-1:
            cell['possible_jumps'] = min(cell['Al'], array[i-1]['Cu'])
        else:
            cell['possible_jumps'] = max(min(cell['Al'], array[i+1]['Cu']),
                                         min(cell['Al'], array[i-1]['Cu']))
    return array
    


if __name__ == '__main__':
    array = make_array(al_concentration, atoms_per_box, N)
    print(array)
    print(array[0]['possible_jumps'])
    print(array[1]['possible_jumps'])
    print(array[2]['possible_jumps'])
    print(array[3]['possible_jumps'])
    print(array[4]['possible_jumps'])
    print(array[5]['possible_jumps'])
    print(array[6]['possible_jumps'])
    print(array[7]['possible_jumps'])
    print(array[8]['possible_jumps'])
    print(array[9]['possible_jumps'])
    print(array[9]['Al'])
    print(array[9]['Cu'])

    


