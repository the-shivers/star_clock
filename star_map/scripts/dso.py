import pandas as pd

df = pd.read_csv('star_map/data/dso.csv')
subset = df[df.mag < 14]
subset.to_csv('dso_subset.csv', index=False)