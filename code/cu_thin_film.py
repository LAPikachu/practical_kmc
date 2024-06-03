import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd

# constants
diff_coeff = 1.49e-7  # diffusion coefficient
ActivationEnergy_Cu = 134.5e3  
distance = 2.54e-10 
R = 8.314  # general gas constant J/(mol*K)
cu_thickness = 100e-9 

number_of_atoms_in_y = cu_thickness / distance 

# parameters
grid_dim_y = int(number_of_atoms_in_y)
grid_dim_x = grid_dim_y // 10
total_time = 10 
initial_concentration = 0.1  
al_concentration = 0.2  
initial_temperature = 298 
end_temperature = 600  
num_steps = 100
time = 0 



def diffusion_coefficient_cu(temperature):
    return diff_coeff * np.exp(-ActivationEnergy_Cu / (R * temperature))

def generate_cu_al_grid(percent_aluminum):
    
    grid = np.zeros((grid_dim_y, grid_dim_x))

    grid[:grid_dim_y//2,:] = np.random.choice([0, 1], size=(grid_dim_y//2, grid_dim_x), p=[1-percent_aluminum, percent_aluminum])

    return grid


def get_coordinates_and_type_of_neighbors(grid, row, col):
    nearby_atoms = []
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)] # down, up, right, left
    for direction in directions:
        i = row + direction[0]
        j = col + direction[1]
        if 0 <= i < grid_dim_y and 0 <= j < grid_dim_x:
            nearby_atoms.append([(i, j), grid[i, j]])
    return nearby_atoms


def get_atom_jump_site_matrix(grid):
    jump_site_list = []
    new_grid = np.zeros((grid_dim_y, grid_dim_x), dtype=object)
    for i in range(grid_dim_y):
        for j in range(grid_dim_x):
            current_type = grid[i, j]
            neighbors = get_coordinates_and_type_of_neighbors(grid, i, j)
            possible_jump_sites = [neighbor[0] for neighbor in neighbors if neighbor[1] != current_type]
            for possible_jump_site in possible_jump_sites: # make sure jump sites are unique
                if possible_jump_site not in jump_site_list: 
                    jump_site_list.append(possible_jump_site)
            no_possible_jump_sites = len(possible_jump_sites)
            new_grid[i, j] = [current_type, possible_jump_sites, no_possible_jump_sites]
            
    return new_grid, jump_site_list

def get_type_matrix(grid):
    values = np.zeros((grid_dim_y, grid_dim_x))
    for i in range(grid_dim_y):
        for j in range(grid_dim_x):
            values[i, j] = grid[i, j][0]
    return values

def get_updated_matrix(atom_jump_site_matrix, jump_from, jump_to):
    jump_from_type = atom_jump_site_matrix[jump_from[0], jump_from[1]][0]
    jump_to_type = atom_jump_site_matrix[jump_to[0], jump_to[1]][0]
    atom_jump_site_matrix[jump_from[0], jump_from[1]][0] = jump_to_type
    atom_jump_site_matrix[jump_to[0], jump_to[1]][0] = jump_from_type
    coordinates_to_update = [jump_from, jump_to]
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for direction in directions:
        i = jump_from[0] + direction[0]
        j = jump_from[1] + direction[1]
        if 0 <= i < grid_dim_y and 0 <= j < grid_dim_x:
            if (i, j) not in coordinates_to_update:
                coordinates_to_update.append((i, j))
    for direction in directions:
        i = jump_to[0] + direction[0]
        j = jump_to[1] + direction[1]
        if 0 <= i < grid_dim_y and 0 <= j < grid_dim_x:
            if (i, j) not in coordinates_to_update:
                coordinates_to_update.append((i, j))
    return atom_jump_site_matrix


def make_jumps(atom_jump_site_matrix, jump_site_list, no_of_jumps):
    while no_of_jumps > 0:
        jump_from = random.choice(jump_site_list)
        jump_to_list = atom_jump_site_matrix[jump_from[0], jump_from[1]][1]
        jump_to = random.choice(jump_to_list)
        if jump_to in jump_site_list:
            jump_site_list.remove(jump_to)
        if jump_from in jump_site_list:
            jump_site_list.remove(jump_from)
        atom_jump_site_matrix = get_updated_matrix(atom_jump_site_matrix, jump_from, jump_to)

        no_of_jumps -= 1
        type_matrix = get_type_matrix(atom_jump_site_matrix)
    return atom_jump_site_matrix

def kmc_sim(time ,total_time, temperature, end_temperature, lattice, num_steps=100):
    current_temperature = temperature
    current_time = time
    
    atom_jump_site_matrix, jump_site_list = get_atom_jump_site_matrix(lattice)

    no_possible_jumps = int(len(jump_site_list)//2)

    jump_rate = 4*diffusion_coefficient_cu(temperature) / distance**2

    random_nr_1 = random.uniform(0, 1)

    no_actual_jumps = int(no_possible_jumps * random_nr_1)

    
    random_nr_2 = random.uniform(0, 1)
    t_ij = -np.log(random_nr_2) / (no_actual_jumps*jump_rate)

    if t_ij < total_time/num_steps: # check if the time is small enough
        current_time += t_ij 
        atom_jump_site_matrix = make_jumps(atom_jump_site_matrix, jump_site_list, no_actual_jumps)

    else: # if the time is too large, perform no jump at all
        no_actual_jumps = 0
        t_ij = total_time/num_steps     
        current_time += t_ij

    # update temperature
    current_temperature = initial_temperature + ((end_temperature - initial_temperature) * (current_time / total_time))
    print('current time: ', current_time, 'current temperature: ', current_temperature, 'jumps executed: ', no_actual_jumps, 'time step: ', t_ij)
    lattice = get_type_matrix(atom_jump_site_matrix)
    return current_time, current_temperature, lattice


if __name__ == '__main__':
    # Generate CuAl grid
    lattice = generate_cu_al_grid(al_concentration)
    time = 0
    save_list = []
    temperature = 293.0
    counter = 0
    while temperature < end_temperature:
        current_time, current_temperature, lattice = kmc_sim(time, total_time, temperature, end_temperature,
                                                            lattice, num_steps)
        time = current_time
        temperature = current_temperature
        counter += 1
        if counter % 10 == 0:
            save_list = [time, temperature, lattice]
            lattice.tofile(f'dump/lattice_t_{time:.2f}_temp_{temperature:.2f}.bin')
    