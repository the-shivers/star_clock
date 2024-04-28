import pandas as pd
import json

"""
The ultimate goal of this file is to return two SVG files
given three parameters (hemisphere, maximum apparent magnitude, constellation skyculture).
It returns one SVG file for stars, and one SVG file for constellation lines.
"""

# CONSTANTS

CONSTELLATIONS_LOC = 'star_map/data/sky_cultures/western_rey.json'
MAG_LIMIT = 9 # 6.5 is naked eye mag limit
STAR_DATA_LOC = 'star_map/data/star_databases/athyg_24_reduced_m10.csv'
DATA_TYPES = {
    'hip': 'Int64', # Nullable integer, for stars with no HIPPARCOS ID
    'proper': str,
    'ra': float,
    'dec': float,
    'mag': float,
    'ci': float
}


# FUNCTIONS

def get_star_data(star_data_loc, data_types):
    return pd.read_csv(star_data_loc, usecols = data_types.keys(), dtype=data_types)

def filter_star_data(star_data, mag_limit):
    filtered_data = star_data.dropna(subset=['hip'])
    filtered_data = filtered_data[filtered_data['mag'] <= mag_limit]
    return filtered_data

def get_stars_dict(star_data_loc, data_types, mag_limit):
    star_df = pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)
    filtered_star_df = star_df.dropna(subset=['hip'])
    filtered_star_df = filtered_star_df[filtered_star_df['mag'] <= mag_limit]
    stars_dict = {int(row['hip']): Star(int(row['hip']), row['proper'], row['ra'], row['dec'], row['mag'], row['ci'])
                for row in filtered_star_df.to_dict('records')}
    return stars_dict


# CLASSES

class Star:
    def __init__(self, hip, proper, ra, dec, mag, ci):
        self.hip = hip
        self.proper = proper
        self.ra = ra
        self.dec = dec
        self.mag = mag
        self.ci = ci

    def __repr__(self):
        return (
            f"Star(hip={self.hip!r}, proper={self.proper!r}, "
            f"ra={self.ra!r}, dec={self.dec!r}, mag={self.mag!r}, ci={self.ci!r})"
        )

    def __str__(self):
        formatted_ra = f"{self.ra:.2f}°"
        formatted_dec = f"{self.dec:+.2f}°"
        star_info = f"{self.proper} (HIP {self.hip})" if self.proper else f"HIP {self.hip}"
        return f"{star_info}: Mag {self.mag:.2f}, RA {formatted_ra}, Dec {formatted_dec}"


class Constellation:
    def __init__(self, abbrv, common_name, latin_name, seg_count, seg_list=None):
        self.abbrv = abbrv
        self.common_name = common_name
        self.latin_name = latin_name
        self.seg_count = seg_count
        self.seg_list = seg_list if seg_list is not None else []
        self._unique_stars = None

    def get_unique_stars(self):
        if self._unique_stars is None:  # Calculate only if needed
            stars = [star for segment in self.seg_list for star in segment]
            self._unique_stars = set(stars)
        return self._unique_stars
    
    def __str__(self):
        unique_stars = self.get_unique_stars()
        return (
            f"{self.common_name} ({self.abbrv}): Unique Stars: {len(unique_stars)}, "
            f"Segments: {len(self.seg_list)}"
        )
    
    def __repr__(self):
        return (
            f"Constellation(abbrv={self.abbrv!r}, common_name={self.common_name!r}, "
            f"latin_name={self.latin_name!r}, seg_count={self.seg_count!r}, seg_list={self.seg_list!r}"
        )

    
class ConstellationParser:
    def __init__(self, filepath, stars_dict):
        self.filepath = filepath
        self.stars_dict = stars_dict

    def parse(self):
        with open(self.filepath, 'r') as file:
            data = json.load(file)['constellations']

        constellations = []
        for item in data:
            abbrv = item['id'].split()[-1]
            if len(abbrv) >= 4: # Indicates a subconstellation, e.g. ursa major paws
                continue
            common_name = item['common_name']['english']
            latin_name = item['common_name'].get('native', None)
            seg_list = []
            for line in item['lines']:
                seg_list.extend([(self.stars_dict.get(line[i]), self.stars_dict.get(line[i + 1]))
                                 for i in range(len(line) - 1)])
            seg_count = len(seg_list)
            constellations.append(Constellation(abbrv, common_name, latin_name, seg_count, seg_list))
        return constellations
    

class Constellationship:
    def __init__(self, constellations, name):
        self.constellations = constellations
        self.name = name
        self._unique_stars = None

    def get_unique_stars(self):
        if self._unique_stars is None:  # Calculate only if needed
            unique_stars = set()
            for constellation in self.constellations:
                unique_stars.update(constellation.get_unique_stars())
            self._unique_stars = unique_stars
        return self._unique_stars
    
    def get_dimmest_star(self):
        unique_stars = self.get_unique_stars()
        dimmest_star = None
        min_mag = float('-26') # The sun, brightest apparent star.
        for star in unique_stars:
            print(star)
            if star.mag > min_mag:
                min_mag = star.mag
                dimmest_star = star
        return dimmest_star


if __name__ == '__main__':
    stars_dict = get_stars_dict(STAR_DATA_LOC, DATA_TYPES, MAG_LIMIT)
    parser = ConstellationParser(CONSTELLATIONS_LOC, stars_dict)
    constellations = parser.parse()
    constellationship = Constellationship(constellations, 'rey')
    print(constellationship.get_dimmest_star())

