import pandas as pd

df = pd.DataFrame({'A': [(0,0), (0,1), (1,0), (1,1) ],
                   'B': [[(0,3), (0,1)], [(0,3), (0,1)], [(0,0), (0,0)], [(0,0), (0,3)]],
                   'C':[2,2,2,2]})
print(df)
value_to_remove = (0,4)

indexes_to_update = df.query['B'].apply(lambda x: value_to_remove in x).tolist() 
print(indexes_to_update)
for index in indexes_to_update:
    jump_sites = df.loc[index, 'B']
    print(jump_sites)
    jump_sites.remove(value_to_remove)
    df.at[index, 'B'] = jump_sites

print(df)


