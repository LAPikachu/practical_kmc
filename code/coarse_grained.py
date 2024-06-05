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
R = 8.314  # general gas constant J/(mol*K)

# parameters
N = 20
L = sample_thickness/N
atoms_per_box = int((L/atomic_distance) * (sample_width/atomic_distance))

# time parameters
total_time = 1e6
default_time_step = total_time/1000

# temperature parameters
end_temperature = 600

def diffusion_coefficient_cu(temperature):
    return diff_coeff * np.exp(-ActivationEnergy_Cu / (R * temperature))


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

def get_jump_site_array(array):
    jump_site_array = []
    for i, cell in enumerate(array):
        if cell['possible_jumps'] > 0:
            jump_site_array.append(i)
    return jump_site_array

def get_total_possible_jumps(array, jump_site_array):
    total_possible_jumps = 0
    for index in jump_site_array:
        total_possible_jumps += array[index]['possible_jumps']
    return total_possible_jumps

def update_indices(array, indices):
    for index in indices:
        if index < 0 or index >= N-1:
            continue
        elif index == 0:
            array[index]['possible_jumps'] = min(array[index]['Al'], array[index+1]['Cu'])
        elif index == N-1:
            array[index]['possible_jumps'] = min(array[index]['Al'], array[index-1]['Cu'])
        else:
            array[index]['possible_jumps'] = max(min(array[index]['Al'], array[index+1]['Cu']),
                                         min(array[index]['Al'], array[index-1]['Cu']))
    return array

def perfom_jump(array, index, jump_index, jump_site_array):
    if array[jump_index]['Al'] == 0: # add to jump_site_array, if it was empty before
        jump_site_array.append(jump_index) 
    array[index]['Al'] -= 1
    array[index]['Cu'] += 1
    array[jump_index]['Al'] += 1
    array[jump_index]['Cu'] -= 1
    if array[index]['Al'] == 0: # remove from jump_site_array, if it is empty now
        jump_site_array.remove(index)
    indices = [index, jump_index]
    if index < jump_index:
        indices.append(index-1)
        indices.append(jump_index+1)
    else:
        indices.append(index+1)
        indices.append(jump_index-1)
    array = update_indices(array, indices)
    return array, jump_site_array

def make_jumps(array, jump_site_array, jumps):
    while jumps > 0:
        random_index = int(random.choice(jump_site_array))
        if random_index == 0:
            random_direction = int(random.choice([-1, 1]))
            if random_direction == -1:
                continue
            jump_index= int(random_index + 1)
            array, jump_site_array = perfom_jump(array, random_index, jump_index, jump_site_array)
            jumps -= 1
        elif random_index == N-1:
            random_direction = int(random.choice([-1, 1]))
            if random_direction == 1:
                continue
            jump_index = int(random_index-1)
            array, jump_site_array = perfom_jump(array, random_index, jump_index, jump_site_array)
            jumps -= 1
        else:
            random_direction = int(random.choice([-1, 1]))
            jump_index = int(random_index + random_direction)
            array, jump_site_array = perfom_jump(array, random_index, jump_index, jump_site_array)
            jumps -= 1
    return array, jump_site_array

def coarse_grained_kmc(array, jump_site_array, current_time, total_time, temprature):
    no_of_possible_jumps_total = get_total_possible_jumps(array, jump_site_array)
    jump_rate = 4*diffusion_coefficient_cu(temprature)/L**2
    random_number_1 = random.uniform(0, 1)
    no_actual_jumps = int(no_of_possible_jumps_total * random_number_1)
    random_number_2 = random.uniform(0, 1)
    t_ij = -np.log(random_number_2)/(jump_rate*no_of_possible_jumps_total)

    if current_time + t_ij < total_time:
        current_time += t_ij
        jumps = no_actual_jumps
        array, jump_site_array = make_jumps(array, jump_site_array, jumps)
    else:
        no_actual_jumps = 0
        current_time += default_time_step
    
    temperature = start_temperature + ((end_temperature - start_temperature) * current_time / total_time)
    print(f'Current time: {current_time}, Current temperature: {temperature}, jumps: {no_actual_jumps}')
    
    return array, current_time, temperature



if __name__ == '__main__':
    array = make_array(al_concentration, atoms_per_box, N)
    jump_site_array = get_jump_site_array(array)
    print(array[1])
    current_time = 0
    start_temperature = 293
    temperature = start_temperature
    counter = 0
    while temperature < 600:
        array, current_time, temperature = coarse_grained_kmc(array, jump_site_array, current_time, total_time, temperature)
        counter += 1
        if counter % 100 == 0:
            al_array =  [cell['Al'] for cell in array]
            with open(f'coarse_dump_02/list_t_{current_time:.2f}_temp_{temperature:.2f}.txt', 'w') as f:
                for item in al_array:
                    f.write(f"{item} \n") 
    print(f'steps: {counter} \n')
    print('cutoff preventions: ')


    


