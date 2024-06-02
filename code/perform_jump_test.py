import numpy as np
import random

grid_dim_y = 10
grid_dim_x = 5
percent_aluminum = 0.5

def generate_cu_al_grid(percent_aluminum):
    num_aluminum_atoms = int(percent_aluminum * grid_dim_x*grid_dim_y) # percentage of aluminum atoms
    
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

test_grid = generate_cu_al_grid(percent_aluminum)
new_grid = np.zeros((grid_dim_y, grid_dim_x), dtype=object)

def get_atom_jump_site_matrix(grid):
    jump_site_list = []
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
    pass # go through the 8 elements that need to be updated


def make_jumps(atom_jump_site_matrix, jump_site_list, no_of_jumps):
    while no_of_jumps > 0:
        jump_from = random.choice(jump_site_list, k=1)
        jump_to_list = atom_jump_site_matrix[jump_from[0], jump_from[1]][1]
        jump_to = random.choice(jump_to_list, k=1)
        jump_site_list.remove(jump_to)
        jump_site_list.remove(jump_from)
        atom_jump_site_matrix = get_updated_matrix(atom_jump_site_matrix, jump_from, jump_to)

        no_of_jumps -= 1


atom_jump_site_matrix, jump_site_list = get_atom_jump_site_matrix(test_grid)


print(jump_site_list)

values_matrix = get_type_matrix(atom_jump_site_matrix)

print(values_matrix)