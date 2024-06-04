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
total_time = 5 
al_concentration = 0.1  
initial_temperature = 500
end_temperature = 600  
num_steps = 10
time = 0 



def diffusion_coefficient_cu(temperature):
    return diff_coeff * np.exp(-ActivationEnergy_Cu / (R * temperature))

def generate_cu_al_grid(percent_aluminum):
    
    grid = np.zeros((grid_dim_y, grid_dim_x))

    grid[:grid_dim_y//2,:] = np.random.choice([0, 1], size=(grid_dim_y//2, grid_dim_x), p=[1-percent_aluminum, percent_aluminum])

    return grid


def get_coordinates_of_neighbors(row, col):
    nearby_atoms = []
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)] # down, up, right, left
    for direction in directions:
        i = row + direction[0]
        j = col + direction[1]
        if 0 <= i < grid_dim_y and 0 <= j < grid_dim_x:
            nearby_atoms.append((i, j))
    return nearby_atoms


def get_atom_jump_site_matrix(grid):
    jump_site_list = []
    new_grid = np.zeros((grid_dim_y, grid_dim_x), dtype=object)
    for i in range(grid_dim_y):
        for j in range(grid_dim_x):
            neighbors = []
            possible_jump_sites = []
            no_possible_jump_sites = 0
            current_type = grid[i, j]
            neighbors = get_coordinates_of_neighbors(i, j)
            possible_jump_sites = [neighbor for neighbor in neighbors if grid[neighbor]  != current_type]
            for possible_jump_site in possible_jump_sites: # make sure jump sites are unique
                if possible_jump_site not in jump_site_list: 
                    jump_site_list.append(possible_jump_site)
            no_possible_jump_sites = len(possible_jump_sites)
            new_grid[i, j] = {'type': current_type, 'possible_jump_sites_coordinates': possible_jump_sites,
                              'no_possible_jump_sites': no_possible_jump_sites,'neighbors_coordinates': neighbors}
            
    return new_grid, jump_site_list

def get_type_matrix(grid):
    values = np.zeros((grid_dim_y, grid_dim_x))
    for i in range(grid_dim_y):
        for j in range(grid_dim_x):
            values[i, j] = grid[i, j]['type']
    return values

def get_updated_matrix(atom_jump_site_matrix, jump_site_list, jump_from, jump_to):
    #make the jump
    jump_from_type = atom_jump_site_matrix[jump_from]['type']
    jump_to_type = atom_jump_site_matrix[jump_to]['type']
    atom_jump_site_matrix[jump_from[0], jump_from[1]]['type'] = jump_to_type
    atom_jump_site_matrix[jump_to[0], jump_to[1]]['type'] = jump_from_type
    # create a list of coordinates to update
    coordinates_to_update = [jump_from, jump_to]
    for coordinates in atom_jump_site_matrix[jump_from]['neighbors_coordinates']:
        if coordinates not in coordinates_to_update:
            coordinates_to_update.append(coordinates)
    for coordinates in atom_jump_site_matrix[jump_to]['neighbors_coordinates']:
        if coordinates not in coordinates_to_update:
            coordinates_to_update.append(coordinates)
    # upate coordinates 
    for coordinates in coordinates_to_update:
        neighbors = []
        possible_jump_sites = []
        no_possible_jump_sites = 0
        current_type = atom_jump_site_matrix[coordinates]['type']
        neighbors = atom_jump_site_matrix[coordinates]['neighbors_coordinates']
        possible_jump_sites = [neighbor for neighbor in neighbors if atom_jump_site_matrix[neighbor]['type'] != current_type]
        no_possible_jump_sites = len(possible_jump_sites)
        atom_jump_site_matrix[coordinates]['possible_jump_sites_coordinates'] = possible_jump_sites
        atom_jump_site_matrix[coordinates]['no_possible_jump_sites'] = no_possible_jump_sites
        for possible_jump_site in possible_jump_sites:
            if possible_jump_site not in jump_site_list:
                jump_site_list.append(possible_jump_site)
    return atom_jump_site_matrix, jump_site_list

def make_sure_to_find_tuple(jump_from, jump_site_list):
    if type(jump_from) != tuple:
        jump_from = random.choice(jump_site_list)
        make_sure_to_find_tuple(jump_from, jump_site_list)
    return jump_from

def get_jump_to_list(jump_from, jump_site_list, atom_jump_site_matrix):
    jump_to_list = atom_jump_site_matrix[jump_from]['possible_jump_sites_coordinates']
    if len(jump_to_list) < 1:
        jump_site_list.remove(jump_from)
        jump_from = random.choice(jump_site_list)
        jump_to_list, jump_site_list = get_jump_to_list(jump_from, jump_site_list, atom_jump_site_matrix)

    if len(jump_site_list) < 1:
        raise ValueError('System is stuck, no more jumps possible.')
    return jump_to_list, jump_site_list

def make_jumps(atom_jump_site_matrix, jump_site_list, no_of_jumps):
    while no_of_jumps > 0:
        jump_from = random.choice(jump_site_list)
        jump_to_list, jump_site_list = get_jump_to_list(jump_from, jump_site_list, atom_jump_site_matrix) 
        
        jump_to = random.choice(jump_to_list)
        if jump_to in jump_site_list:
            jump_site_list.remove(jump_to)
        if jump_from in jump_site_list:
            jump_site_list.remove(jump_from)
        atom_jump_site_matrix, jump_site_list = get_updated_matrix(atom_jump_site_matrix, jump_site_list, jump_from, jump_to)

        no_of_jumps -= 1
    return atom_jump_site_matrix, jump_site_list

def kmc_sim(time ,total_time, temperature, end_temperature, atom_jump_site_matrix, jump_site_list , num_steps=10):
    current_temperature = temperature
    current_time = time
    
    no_possible_jumps = int(len(jump_site_list)//2)

    jump_rate = 4*diffusion_coefficient_cu(temperature) / distance**2

    random_nr_1 = random.uniform(0, 1)

    no_actual_jumps = int(no_possible_jumps * random_nr_1)

    random_nr_2 = random.uniform(0, 1)
    t_ij =-np.log(random_nr_2) / (no_actual_jumps*jump_rate) #calculate time step

    if t_ij < total_time/num_steps: # check if the time is small enough
        current_time += t_ij 
        atom_jump_site_matrix, jump_site_list = make_jumps(atom_jump_site_matrix, jump_site_list, no_actual_jumps)

    else: # if the time is too large, perform no jump at all
        no_actual_jumps = 0
        t_ij = total_time/num_steps     
        current_time += t_ij

    # update temperature
    current_temperature = initial_temperature + ((end_temperature - initial_temperature) * (current_time / total_time))
    print('current time: ', current_time, 'current temperature: ', current_temperature, 'jumps executed: ', no_actual_jumps, 'time step: ', t_ij)
    return current_time, current_temperature, atom_jump_site_matrix, jump_site_list


if __name__ == '__main__':
    # Generate CuAl grid
    lattice = generate_cu_al_grid(al_concentration)
    time = 0
    save_list = []
    temperature = 293.0
    counter = 0
    while temperature < end_temperature:
        atom_jump_site_matrix, jump_site_list = get_atom_jump_site_matrix(lattice)
        current_time, current_temperature, atom_jump_site_matrix, jump_site_list = kmc_sim(time, total_time, temperature, end_temperature,
                                                            atom_jump_site_matrix, jump_site_list, num_steps)
        time = current_time
        temperature = current_temperature
        counter += 1
        if counter % 30 == 0:
            lattice = get_type_matrix(atom_jump_site_matrix)
            save_list = [time, temperature, lattice]
            lattice.tofile(f'dump/lattice_t_{time:.2f}_temp_{temperature:.2f}.bin')
    