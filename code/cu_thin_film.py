import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd
import json

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
al_concentration = 0.1  
initial_temperature = 298 
end_temperature = 600  
num_steps = 100
time = 0 



def diffusion_coefficient_cu(temperature):
    return diff_coeff * np.exp(-ActivationEnergy_Cu / (R * temperature))


def generate_cu_al_grid(percent_aluminum):
    num_aluminum_atoms = int(percent_aluminum * grid_dim_x*grid_dim_y) # percentage of aluminum atoms
    
    grid = np.zeros((grid_dim_y, grid_dim_x))

    grid[:grid_dim_y//2,:] = np.random.choice([0, 1], size=(grid_dim_y//2, grid_dim_x), p=[1-percent_aluminum, percent_aluminum])

    return grid

def get_coordinates_of_neighbor_type(grid, row, col, type=0):
    nearby_atoms = []
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for direction in directions:
        i = row + direction[0]
        j = col + direction[1]
        if 0 <= i < grid_dim_y and 0 <= j < grid_dim_x and grid[i, j] == type:
            nearby_atoms.append((i, j))
    return nearby_atoms


def get_coordinates_atoms_type(grid, type = 1):
    coordinates = []
    for i in range(0,grid_dim_y):
        for j in range(0, grid_dim_x):
            if grid[i, j] == type:
                coordinates.append((i, j))
    return coordinates

def get_dataframe_with_jump_rates(lattice):
    al_coordinates = get_coordinates_atoms_type(lattice, type=1)
    possible_jump_sites = []
    no_possible_jump_sites = []
    for atom_coordinates in al_coordinates:
        row, col = atom_coordinates
        current_possible_jump_sites = get_coordinates_of_neighbor_type(lattice, row, col, type=0)
        no_possible_jump_sites.append(len(current_possible_jump_sites))
        possible_jump_sites.append(current_possible_jump_sites)
    
    total_possible_jump_sites = sum(no_possible_jump_sites)

    al_atoms_dict = {'coordinates':al_coordinates ,
                     'jump_sites_coordinates': possible_jump_sites, 
                     'no_jump_sites': no_possible_jump_sites}    
    al_atoms_df = pd.DataFrame(al_atoms_dict)
    return al_atoms_df

def make_jumps(lattice, al_atoms_df, jumps):
        jump_sites = al_atoms_df[al_atoms_df['no_jump_sites'] > 0] # only jump sites with possible jumps
        jump_site = jump_sites.sample(1) # choose a random jump site
        jump_from = jump_site['coordinates'].values
        jump_to_list = jump_site['jump_sites_coordinates'].values
        jump_to = random.choice(jump_to_list)
        lattice[jump_from[0][0],jump_from[0][1]] = 0
        lattice[jump_to[0][0],jump_to[0][1]] = 1
        jumps -=1
        return lattice, jumps



def kmc_sim(time ,total_time, temperature, end_temperature, lattice, num_steps=100):
    current_temperature = temperature
    current_time = time

    neighbor_number = 4  # for 2D

    al_atoms_df = get_dataframe_with_jump_rates(lattice)

    total_possible_jump_sites = al_atoms_df['no_jump_sites'].sum()

    jump_rate = 4*diffusion_coefficient_cu(temperature) / distance**2

    random_nr_1 = random.uniform(0, 1)

    jumps = total_possible_jump_sites * random_nr_1

    lattice, jumps = make_jumps(lattice, al_atoms_df, jumps) # first jump
    
    random_nr_2 = random.uniform(0, 1)
    t_ij = -np.log(random_nr_2) / (total_possible_jump_sites*jump_rate)

    if t_ij < total_time/num_steps: # check if the time is small enough
        current_time += t_ij 
        jumps_executed = jumps
        while jumps >0: # consecutive jumps
            al_atoms_df = get_dataframe_with_jump_rates(lattice) # always update the dataaframe, not to jump a site twice
            lattice, jumps = make_jumps(lattice, al_atoms_df, jumps) 
            jumps -= 1
    else: # if the time is too large, perform no jump at all
        jumps_executed = 0
        t_ij = total_time/num_steps     
        current_time += t_ij

    # update temperature
    current_temperature = end_temperature * (current_time / total_time)
    print('current time: ', current_time, 'current temperature: ', current_temperature, 'current jumps: ', jumps_executed, 'time step: ', t_ij)
    
    return current_time, current_temperature, lattice


if __name__ == '__main__':
    # Generate CuAl grid
    lattice =  np.zeros((grid_dim_y,grid_dim_x)) #generate_cu_al_grid(al_concentration) # 0: Cu, 1: Al
    lattice[:grid_dim_y//2,:] = np.random.choice([0, 1], size=(grid_dim_y//2, grid_dim_x), p=[1-al_concentration, al_concentration])
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis('off')
    time = 0
    save_list = []
    temperature = initial_temperature
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
    
    fig, axs = plt.subplots(ncols=10, figsize=(20, 2))
    for i, time, temperature, lattice in enumerate(save_list):
        ax.axis('off')
        ax.set_title('t = {time} s T = {temperature}'.format(time=time, temperature=temperature))
        sns.heatmap(lattice, cmap='Blues' , cbar=False, square=True, ax=axs[i], yticklabels=True, xticklabels=True)