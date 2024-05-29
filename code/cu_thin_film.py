import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import pandas as pd

# parameters
grid_dim_y = 10
grid_dim_x = 10
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

# Funktion zur Berechnung des Diffusionskoeffizienten für Cu
def diffusion_coefficient_cu(temperature):
    return diff_coeff * np.exp(-ActivationEnergy_Cu / (R * temperature))


def generate_cu_al_grid(N, percent_aluminum):
    num_aluminum_atoms = int(percent_aluminum * N) # percentage of aluminum atoms
    
    grid = np.zeros((grid_dim_y, grid_dim_x))

    grid[:grid_dim_y//2,:] = np.random.choice([0, 1], size=(grid_dim_y//2, grid_dim_x), p=[1-percent_aluminum, percent_aluminum])

    return grid

def get_coordinates_of_nearby_type(grid, row, col, type=1):
    nearby_atoms = []
    for i in range(row-1, row+2):
        for j in range(col-1, col+2):
            if 0 <= i < 10 and 0 <= j < 10 and (i != row or j != col):
                if grid[i, j] == type:
                    nearby_atoms.append((i, j))
    return nearby_atoms


def get_coordinates_atoms_type(grid, type = 0):
    coordinates = []
    for i in range(grid_dim_y):
        for j in range(grid_dim_x):
            if grid[i, j] == 0:
                coordinates.append((i, j))
    return coordinates

def kmc_diff_sim(TotalTime, initial_concentration, cu_thickness, al_concentration, initial_temperature, end_temperature, num_steps=1000):
    current_time = t = 0
    T = initial_temperature
    current_temperature = initial_temperature
    current_concentration_cu = np.ones(num_steps+1) * initial_concentration * 100  # 
    delta_t = TotalTime / num_steps
    
    # Generate CuAl grid
    lattice = generate_cu_al_grid(100, al_concentration) # 0: Cu, 1: Al
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(lattice, cmap='Blues' , cbar=False, square=True, ax=ax, yticklabels=True, xticklabels=True)
    plt.show()
    
    neighbor_number = 4  # for 2D

    for i in range(num_steps):
        cu_coordinates = get_coordinates_atoms_type(lattice)
        # number of Cu atoms 
        number_cu_neighbor = np.zeros(len(cu_coordinates))
        for j, (row, col) in enumerate(cu_coordinates):
            surroundings = get_coordinates_of_nearby_type(lattice, row, col, type=1)
            number_cu_neighbor[j] = neighbor_number - len(surroundings)
        
        # Berechne Sprungraten für Cu-Atome
        k_cu = diffusion_coefficient_cu(current_temperature)
        k_cu_mat = number_cu_neighbor * k_cu
        
        # Führe Diffusion oder Rückkehr für jedes Cu-Atom durch
        for j, (row, col) in enumerate(cu_coordinates):
            # Berechne Sprungraten für Diffusion und Rückkehr
            k_diffusion = k_cu_mat[j]
            k_reversion = neighbor_number * diffusion_coefficient_cu(current_temperature)
            r_diffusion = random.uniform(0, 1)
            r_reversion = random.uniform(0, 1)
            
            # Überprüfe, ob das aktuelle Cu-Atom diffundiert
            if r_diffusion < k_diffusion / (k_diffusion + k_reversion):
                jump_direction = random.choice([(1, 0), (-1, 0), (0, 1), (0, -1)])
                new_position = ((row + jump_direction[0]) % 10, (col + jump_direction[1]) % 10)
                
                lattice[row, col] = 'Al'
                lattice[new_position] = 'Cu'
                # Aktualisiere die Konzentration nur für das diffundierte Atom
                current_concentration_cu[i+1] = current_concentration_cu[i]
            elif r_reversion < k_reversion / (k_diffusion + k_reversion):
                lattice[row, col] = 'Al'
                # Aktualisiere die Konzentration für Nicht-Schritte
                current_concentration_cu[i+1] = current_concentration_cu[i] + 0.01
        
        # Zeit und Temperatur aktualisieren
        tij = -np.log(random.uniform(0, 1)) / np.sum(k_cu_mat)
        t += tij
        vT = R * T  # vT hier berechnen
        T += vT * tij
        
        # Temperatur aktualisieren
        current_temperature = initial_temperature + (end_temperature - initial_temperature) * (current_time / TotalTime)
        current_time += delta_t
    
    time = np.linspace(0, TotalTime, num_steps+1)
    return time, current_concentration_cu


# Simulation durchführen
time, concentration_cu = kmc_diff_sim(TotalTime, initial_concentration, cu_thickness, al_concentration, initial_temperature, end_temperature, num_steps)

# Ergebnisse plotten
plt.figure(figsize=(10, 6))
plt.plot(time, concentration_cu, color='blue', linestyle='-', label='Cu-Konzentration in CuAl-Dünnschicht')
plt.axhline(y=initial_concentration * 100, color='red', linestyle='--', label='Anfangskonzentration Cu')
plt.xlabel('Zeit (s)')
plt.ylabel('Cu-Konzentration (%)')
plt.title('Entwicklung der Kupferkonzentration in CuAl-Dünnschicht')
plt.legend()
plt.show()
