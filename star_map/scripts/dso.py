import pandas as pd

df = pd.read_csv('star_map/data/dso.csv')
# https://github.com/astronexus/HYG-Database/tree/main/misc

condition = (
    (df['mag'] < 6.5) & (df['mag'].notna()) |  # Objects brighter than mag 6.5 and not NaN
    df['name'].notna() |  # Objects where name is not NaN
    ((df['cat1'] == 'M') | (df['cat2'] == 'M'))  # Objects categorized as Messier in cat1 or cat2
)

filtered_df = df[condition]
filtered_df = filtered_df.sort_values(by='mag', ascending=True)
final_subset = filtered_df[['ra', 'dec', 'type', 'const', 'mag', 'name', 'id', 'id1', 'cat1', 'id2', 'cat2']]
final_subset.to_csv('star_map/data/dso_subset.csv', index=False)

# File manually cleaned then saved as dso_clean.csv
