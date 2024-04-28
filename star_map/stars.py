import pandas as pd

"""
The ultimate goal of this file is to return two SVG files
given three parameters (hemisphere, maximum apparent magnitude, constellation skyculture).
It returns one SVG file for stars, and one SVG file for constellation lines.
"""

# CONSTANTS

CONSTELLATIONS_LOC = 'star_map/data/sky_cultures/rey/constellationship.fab'
CONSTELLATION_NAMES_LOC = 'star_map/data/sky_cultures/rey/constellation_names.eng.fab'
MAG_LIMIT = 8 # 6.5 is naked eye mag limit
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

def load_constellation_names(filepath):
    names_dict = {}
    with open(filepath, 'r', encoding='utf-8') as file:
        for line in file:
            if line.strip().startswith('#'): # skip comments
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2: # skip malformed lines
                continue
            abbrv = parts[0].strip('.')
            common_name = parts[1].strip('"')
            latin_name = parts[2].strip('_(")').rstrip('")')
            names_dict[abbrv] = (common_name, latin_name)
    return names_dict


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
            f"Star(hip={self.hip}, proper={self.proper!r}, "
            f"ra={self.ra}, dec={self.dec}, mag={self.mag}, ci={self.ci})"
        )

    def __str__(self):
        formatted_ra = f"{self.ra:.2f}°"
        formatted_dec = f"{self.dec:+.2f}°"
        star_info = f"{self.proper} (HIP {self.hip})" if self.proper else f"HIP {self.hip}"
        return f"{star_info}: Mag {self.mag:.2f}, RA {formatted_ra}, Dec {formatted_dec}"


class Constellation:
    def __init__(self, abbrv, common_name, latin_name, seg_count, seg_list=None, sub_constellations=None):
        self.abbrv = abbrv
        self.common_name = common_name
        self.latin_name = latin_name
        self.seg_count = seg_count
        self.seg_list = seg_list if seg_list is not None else []
        self.sub_constellations = sub_constellations if sub_constellations is not None else []
        self._unique_stars = None  # Initialize to None

    def add_sub_constellation(self, sub_constellation):
        self.sub_constellations.append(sub_constellation)

    def get_unique_stars(self):
        if self._unique_stars is None:  # Calculate only if needed
            stars = [star for segment in self.seg_list for star in segment]
            self._unique_stars = set(stars)
        return self._unique_stars
    
    def __str__(self):
        unique_stars = self.get_unique_stars()
        sub_constellation_names = ', '.join(sub.common_name for sub in self.sub_constellations)
        return (
            f"{self.common_name} ({self.abbrv}): Unique Stars: {len(unique_stars)}, "
            f"Segments: {len(self.seg_list)}, "
            f"Sub-Constellations: [{sub_constellation_names}]"
        )
    
    def __repr__(self):
        return (
            f"Constellation(abbrv={self.abbrv!r}, common_name={self.common_name!r}, "
            f"latin_name={self.latin_name!r}, seg_count={self.seg_count!r}, seg_list={self.seg_list!r}, "
            f"sub_constellations={[sub.abbrv for sub in self.sub_constellations]!r})"
        )


class ConstellationParser:
    def __init__(self, filepath, names_dict, stars_dict):
        self.filepath = filepath
        self.names_dict = names_dict
        self.stars_dict = stars_dict

    def parse(self):
        constellations = []
        current_main = None # current main const. for dealing with sub-constellations
        with open(self.filepath, 'r') as file:
            for line in file:
                if line.startswith('#'): # skip comments
                    continue
                if not line.startswith('.'): # main constellation logic
                    constellation = self.parse_line(line)
                    if constellation:
                        constellations.append(constellation)
                        current_main = constellation
                else: # Sub-constellation logic (e.g. ursa major paws)
                    sub_constellation = self.parse_line(line)
                    if sub_constellation and current_main:
                        current_main.add_sub_constellation(sub_constellation)
        return constellations
    
    def parse_line(self, line):
        parts = line.strip().split()
        if len(parts) < 3:
            return None
        abbrv = parts[0].strip('.')
        seg_count = int(parts[1])
        seg_list = [(self.stars_dict.get(int(hip1)), self.stars_dict.get(int(hip2))) for hip1, hip2 in zip(parts[2::2], parts[3::2])]
        return Constellation(
            abbrv = abbrv,
            common_name = self.names_dict.get(abbrv)[0],
            latin_name = self.names_dict.get(abbrv)[1],
            seg_count = seg_count,
            seg_list = seg_list,
            sub_constellations = None
        )


if __name__ == '__main__':
    names_dict = load_constellation_names(CONSTELLATION_NAMES_LOC)
    stars_dict = get_stars_dict(STAR_DATA_LOC, DATA_TYPES, MAG_LIMIT)
    parser = ConstellationParser(CONSTELLATIONS_LOC, names_dict, stars_dict)
    constellations = parser.parse()
    print(constellations[0])



# data = data[pd.notnull(data['hip'])]
# data['hip'] = data['hip'].astype(int)
# data = data[data['mag'] <= MAG_LIMIT]
# data = data[COLS]
# data['theta'] = np.radians(data['ra'] * 15)
# data['r'] = abs(data['dec'] / 90)
# north = data[data['dec'] >= 0]







