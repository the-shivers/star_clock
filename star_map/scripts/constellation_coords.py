# Produces a CSV file for constellation labels.
# Just averages the positions of the stars (lol)

import os
import numpy as np
import pandas as pd
from star_map_classes import (
    StarHolder,
    ConstellationParser
)

# Project directory
dir = os.path.expanduser('~/star_clock')

# Star and constellation config
# constellations_json = f'{dir}/star_map/data/western_iau_sky_culture.json'
constellations_json = f'{dir}/star_map/data/sky_cultures/snt/western_snt_sky_culture.json'
mag_limit = 7.5 # For limiting size of stars list.
star_data_loc = f'{dir}/star_map/data/stars/athyg_24_reduced_m10.csv'

starholder = StarHolder(star_data_loc, 7.5)
parser = ConstellationParser(constellations_json, starholder)
constellations = parser.parse()

coords_dict = {}

for constellation in constellations:
    stars = constellation.get_unique_stars()
    sin_sum = cos_sum = dec_total = 0
    for star in stars:
        ra_radians = star.ra * 2 * np.pi / 24
        sin_sum += np.sin(ra_radians)
        cos_sum += np.cos(ra_radians)
        dec_total += star.dec
    count = len(stars)
    mean_ra_radians = np.arctan2(sin_sum / count, cos_sum / count)
    mean_ra = (mean_ra_radians * 24 / (2 * np.pi)) % 24
    mean_dec = dec_total / count
    print(f'{constellation.latin_name}: ra:{mean_ra}, dec:{mean_dec}')
    coords_dict[constellation.latin_name] = {'ra': mean_ra, 'dec': mean_dec, 'rot':0}

df = pd.DataFrame.from_dict(coords_dict, orient='index').reset_index()
df.columns = ['latin_name', 'ra', 'dec', 'rot']
csv_filename = f'{dir}/star_map/data/sky_cultures/snt/constellation_coords_raw.csv'
df.to_csv(csv_filename, index=False)