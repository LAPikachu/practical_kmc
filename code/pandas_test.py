import pandas as pd
al_atoms_dict = {'atom_id': [0, 1, 2, ], 
                     'coordinates': [(0, 0), (0, 1), (1, 0), (1, 1)],
                     'jump_sites_coordinates': [[], [(0,3)], [(2,0)], [(1,2),(2,1)] ], 
                     'no_jump_sites': [0, 1, 1, 2] }   
al_atoms_df = pd.DataFrame(al_atoms_dict)

def find_row_by_coordinates(df, coordinates):
    row = df[df['coordinates'] == coordinates]
    return row

indexes_to_update = df.query['B'].apply(lambda x: value_to_remove in x).tolist() 
print(indexes_to_update)
for index in indexes_to_update:
    jump_sites = df.loc[index, 'B']
    print(jump_sites)
    jump_sites.remove(value_to_remove)
    df.at[index, 'B'] = jump_sites

print(df)


