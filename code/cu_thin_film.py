import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd

# parameters
grid_dim_y = 4
grid_dim_x = 4
TotalTime = 10 
initial_concentration = 0.1  
cu_thickness = 100e-9  
al_concentration = 0.9  
initial_temperature = T = 298 
end_temperature = 600  
num_steps = 1000  
t = 0 

# constants
diff_coeff = 1.49e-7  # diffusion coefficient
ActivationEnergy_Cu = 134.5e3  
distance = 2.54e-10 
R = 8.314  # general gas constant J/(mol*K)

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

def kmc_sim(TotalTime, initial_concentration, cu_thickness, al_concentration, initial_temperature, end_temperature, num_steps=1000):
    current_time = 0
    T = initial_temperature
    current_temperature = initial_temperature
    
    # Generate CuAl grid
    lattice =  np.zeros((grid_dim_y,grid_dim_x)) #generate_cu_al_grid(al_concentration) # 0: Cu, 1: Al
    lattice[:grid_dim_y//2,:] = 1
    fig, ax = plt.subplots(figsize=(6, 6))
    #sns.heatmap(lattice, cmap='Blues' , cbar=False, square=True, ax=ax, yticklabels=True, xticklabels=True)
    #plt.show()
    
    neighbor_number = 4  # for 2D
    al_coordinates = get_coordinates_atoms_type(lattice, type=1)
    possible_jump_sites = []
    no_possible_jump_sites = []
    for atom_coordinates in al_coordinates:
        row, col = atom_coordinates
        current_possible_jump_sites = get_coordinates_of_neighbor_type(lattice, row, col, type=0)
        no_possible_jump_sites.append(len(current_possible_jump_sites))
        possible_jump_sites.append(current_possible_jump_sites)
    
    total_possible_jump_sites = sum(no_possible_jump_sites)    
    



            
    # update time
    t_ij = -np.log(random.uniform(0, 1)) / np.sum(k_cu_mat)
    t += t_ij
    vT = R * T  
    T += vT * t_ij
        
    # update temperature
    current_temperature = initial_temperature + (end_temperature - initial_temperature) * (current_time / TotalTime)
    current_time += delta_t
    
    time = np.linspace(0, TotalTime, num_steps+1)
    return time, atoms


if __name__ == '__main__':
    time, atoms = kmc_sim(TotalTime, initial_concentration, cu_thickness, al_concentration, initial_temperature, end_temperature, num_steps)

    # Ergebnisse plotten
    plt.figure(figsize=(10, 6))
    plt.plot(time, concentration_cu, color='pink', linestyle='--', label='Cu concentration')
    plt.axhline(y=initial_concentration * 100, color='red', linestyle='--', label='Anfangskonzentration Cu')
    plt.title('Entwicklung der Kupferkonzentration in CuAl-Dünnschicht')
    plt.legend()
    plt.show()
