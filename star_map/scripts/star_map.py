"""
The ultimate goal of this file is to return two SVG files
given three parameters (hemisphere, maximum apparent magnitude, constellation skyculture).
It returns one SVG file for stars, and one SVG file for constellation lines.
"""

import pandas as pd
import numpy as np
import svgwrite
import json

CONSTELLATIONS_LOC = 'star_map/data/wester_iau_sky_culture.json'
MAG_LIMIT = 5 # 6.5 is naked eye mag limit
STAR_DATA_LOC = 'star_map/data/athyg_24_reduced_m10.csv'
DATA_TYPES = {
    'hip': 'Int64', # Nullable integer, for stars with no HIPPARCOS ID
    'proper': str,
    'ra': float,
    'dec': float,
    'mag': float,
    'ci': float
}

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
            if star.mag > min_mag:
                min_mag = star.mag
                dimmest_star = star
        return dimmest_star


def get_stars_dict(star_data_loc, data_types, mag_limit):
    """Generate a dictionary for stars with a HIP number within a magnitude limit."""
    star_data = pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)
    filtered_data = star_data.dropna(subset=['hip'])
    filtered_data = filtered_data[filtered_data['mag'] <= mag_limit]
    return {int(row['hip']): Star(int(row['hip']), row['proper'], row['ra'], row['dec'], row['mag'], row['ci'])
        for row in filtered_data.to_dict('records')}

def is_point_inside_circle(point, center, radius):
    x, y = point
    x0, y0 = center
    return (x - x0)**2 + (y - y0)**2 <= radius**2

def find_line_circle_intersections(center, radius, start, end):
    (x0, y0) = center
    (x1, y1) = start
    (x2, y2) = end
    dx, dy = x2 - x1, y2 - y1
    a = dx**2 + dy**2
    b = 2 * (dx * (x1 - x0) + dy * (y1 - y0))
    c = (x1 - x0)**2 + (y1 - y0)**2 - radius**2
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        if is_point_inside_circle(start, center, radius) and is_point_inside_circle(end, center, radius):
            return [start, end]  
        return [] 
    t1 = (-b + np.sqrt(discriminant)) / (2 * a)
    t2 = (-b - np.sqrt(discriminant)) / (2 * a)
    intersections = [(x1 + t * dx, y1 + t * dy) for t in [t1, t2] if 0 <= t <= 1]
    return intersections

def truncate_line_to_circle(center, radius, start, end):
    intersections = find_line_circle_intersections(center, radius, start, end)
    start_inside = is_point_inside_circle(start, center, radius)
    end_inside = is_point_inside_circle(end, center, radius)
    if start_inside and end_inside:
        return [start, end]  # Both points inside
    elif start_inside or end_inside:
        endpoint_inside = start if start_inside else end
        if intersections:
            return [endpoint_inside, intersections[0]]
    elif intersections:
        if len(intersections) == 2:
            return intersections
    return None


def bv_to_color(bv):
    # Define color ranges and corresponding B-V index breakpoints
    colors = ['#94B6FF', '#99B9FF', '#C9D9FF', '#ECEEFF', '#FFFFFF', '#FFF2EE', '#FFE7D2', '#FFCC98']
    breakpoints = [-0.33, -0.15, 0.0, 0.25, 0.58, 0.81, 1.4, float('inf')]  # Open-ended for values beyond 1.4
    # Determine the segment the bv value falls into
    for i in range(len(breakpoints) - 1):
        if breakpoints[i] <= bv < breakpoints[i + 1]:
            # Calculate the normalized position of bv between the current and next breakpoint
            lower_bound = breakpoints[i]
            upper_bound = breakpoints[i + 1]
            t = (bv - lower_bound) / (upper_bound - lower_bound) if upper_bound != float('inf') else 1.0
            # Linear interpolation of hex colors
            color_start = tuple(int(colors[i][j:j+2], 16) for j in range(1, 7, 2))
            if i+1 < len(colors):
                color_end = tuple(int(colors[i+1][j:j+2], 16) for j in range(1, 7, 2))
            else:
                color_end = color_start  # No interpolation beyond the last specified color
            r = int(color_start[0] + t * (color_end[0] - color_start[0]))
            g = int(color_start[1] + t * (color_end[1] - color_start[1]))
            b = int(color_start[2] + t * (color_end[2] - color_start[2]))
            return f"#{r:02x}{g:02x}{b:02x}"
    return '#ffffff'

def line_segment_length(start, end):
    """Calculate the Euclidean distance between two points."""
    (x1, y1) = start
    (x2, y2) = end
    length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return length

def get_northern_hemisphere_cartesian(star, graphic_size, dec_degrees, star_circle_dia):
    theta_radians = star.ra / 24 * np.pi * 2
    hypotenuse = ((90 - star.dec) / dec_degrees) * star_circle_dia/2
    x = (graphic_size - star_circle_dia) / 2 + star_circle_dia/2 + np.sin(theta_radians) * hypotenuse
    y = 100 + star_circle_dia/2 - np.cos(theta_radians) * hypotenuse
    return (x, y)


stars_dict = get_stars_dict(STAR_DATA_LOC, DATA_TYPES, 100) # No mag limit for now.
parser = ConstellationParser(CONSTELLATIONS_LOC, stars_dict)
constellations = parser.parse()
constellationship = Constellationship(constellations, 'iau')

# Okay, let's say we want to parameterize this.
size = 1000 # This is the size of the ENTIRE ILLUSTRATION. 
full_circle_dia = 900 # The circle containing the starscape AND months, tickmarks.
star_circle_dia = 800 # This is the circle containing our stars
dec_degrees = 108 # Number of degrees of declination to map. 90 would be one celestial hemisphere.

# First let's try transforming the SVG and just... see what happens lol.
drawing = svgwrite.Drawing(filename="star_map.svg", size=(f"{size}px", f"{size}px"))
full_circle = drawing.add(drawing.circle(center=(size/2, size/2), r=full_circle_dia/2, fill="#ffffff", stroke='blue'))
star_circle = drawing.add(drawing.circle(center=(size/2, size/2), r=star_circle_dia/2, fill="#112233", stroke='black'))
line1 = drawing.add(drawing.line(start=(size/2, 0), end=(size/2, size), stroke='white', stroke_width=0.1))
line2 = drawing.add(drawing.line(start=(0, size/2), end=(size, size/2), stroke='white', stroke_width=0.1))
line3 = drawing.add(drawing.line(start=(0, 0), end=(size, size), stroke='white', stroke_width=0.1))
line4 = drawing.add(drawing.line(start=(0, size), end=(size, 0), stroke='white', stroke_width=0.1))
max_mag = 5

# Constellations
for constellation in constellationship.constellations:
    for stars in constellation.seg_list:
        # Segment is a tuple of two star objects
        # So we could just get the x, y of each star, draw a line between them.
        start = get_northern_hemisphere_cartesian(stars[0], size, dec_degrees, star_circle_dia)
        end = get_northern_hemisphere_cartesian(stars[1], size, dec_degrees, star_circle_dia)
        # if (start[0] >= 0 and end[0] >= 0) or (end[0] >= 0 and end[0] >= 0):
        if is_point_inside_circle(start, (size/2, size/2), star_circle_dia/2) or is_point_inside_circle(end, (size/2, size/2), star_circle_dia/2):
            endpoints = truncate_line_to_circle((size/2, size/2), star_circle_dia/2, start, end)
            if endpoints:
                if line_segment_length(*endpoints) > 200:
                    print(stars)
                a = drawing.add(drawing.line(endpoints[0], endpoints[1], stroke='white', stroke_width=0.1))

# Stars
for hip, star in stars_dict.items():
    if star.dec < 90 - dec_degrees or star.mag > max_mag:
        continue
    x, y = get_northern_hemisphere_cartesian(star, size, dec_degrees, star_circle_dia)
    star_radius = max_mag - star.mag
    fill = bv_to_color(star.ci)
    a = drawing.add(drawing.circle(center=(x, y), r=star_radius, fill=fill))

def add_combined_paths_to_svg(drawing, transformed_paths):
    combined_path_string = ""
    for path_data in transformed_paths:
        path_string = f"M {path_data['M'][0]},{path_data['M'][1]} "  # Move to the start point
        for bezier in path_data['beziers']:
            path_string += f"C {bezier[1][0]},{bezier[1][1]} {bezier[2][0]},{bezier[2][1]} {bezier[3][0]},{bezier[3][1]} "
        combined_path_string += path_string  # Combine into one path string
    # Create a single SVG path element with the combined path string
    drawing.add(drawing.path(d=combined_path_string, fill="#FFFFFF", fill_opacity=0.1, stroke="none", stroke_width=0, fill_rule="evenodd"))

# Usage
add_combined_paths_to_svg(drawing, transformed_data)


drawing.save()

def create_svg_hemisphere(hemisphere='north', mag_limit=6.5, output='output.svg'):
    pass

