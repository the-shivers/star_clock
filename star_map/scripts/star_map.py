import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
import re
import svgwrite
import json

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

class SVGHemisphere:
    def __init__(self, size, full_circle_dia, star_circle_dia, dec_degrees, filename="star_map.svg", is_north=True):
        self.size = size
        self.full_circle_dia = full_circle_dia
        self.star_circle_dia = star_circle_dia
        self.dec_degrees = dec_degrees
        self.is_north = is_north
        self.drawing = svgwrite.Drawing(filename=filename, size=(f"{size}px", f"{size}px"))
        self.elements = {}

    # SVG-specific geometry helper functions
    def get_hemisphere_cartesian(self, ra, dec):
        """Given RA and DEC of some stellar object, produces cartesian coordinates for the azimuthal equidistant projection."""
        if self.is_north:
            theta_radians = ra / 24 * np.pi * 2
            hypotenuse = ((90 - dec) / self.dec_degrees) * self.star_circle_dia/2
        else:
            theta_radians = 2 * np.pi - (ra / 24 * np.pi * 2)
            hypotenuse = ((90 + dec) / self.dec_degrees) * self.star_circle_dia/2
        x = self.size / 2 + np.sin(theta_radians) * hypotenuse
        y = self.size/2 - np.cos(theta_radians) * hypotenuse
        return (x, y)
    
    def is_point_inside_circle(self, point):
        x, y = point
        x0, y0 = self.size/2, self.size/2
        return (x - x0)**2 + (y - y0)**2 <= (self.star_circle_dia/2)**2
    
    def find_line_circle_intersections(self, start, end):
        (x0, y0) = self.size/2, self.size/2
        (x1, y1) = start
        (x2, y2) = end
        dx, dy = x2 - x1, y2 - y1
        a = dx**2 + dy**2
        b = 2 * (dx * (x1 - x0) + dy * (y1 - y0))
        c = (x1 - x0)**2 + (y1 - y0)**2 - (self.star_circle_dia/2)**2
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            if self.is_point_inside_circle(start) and self.is_point_inside_circle(end):
                return [start, end]  
            return [] 
        t1 = (-b + np.sqrt(discriminant)) / (2 * a)
        t2 = (-b - np.sqrt(discriminant)) / (2 * a)
        intersections = [(x1 + t * dx, y1 + t * dy) for t in [t1, t2] if 0 <= t <= 1]
        return intersections
    
    def truncate_line(self, start, end, abs_amount):
        # Removes absolute amount of length from line segment ends.
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        line_length = np.sqrt(dx**2 + dy**2)
        truncation_ratio = abs_amount / line_length
        new_start_x = start[0] + dx * truncation_ratio
        new_start_y = start[1] + dy * truncation_ratio
        new_end_x = end[0] - dx * truncation_ratio
        new_end_y = end[1] - dy * truncation_ratio
        return ((new_start_x, new_start_y), (new_end_x, new_end_y))
    
    def truncate_line_to_circle(self, start, end):
        # If line segment exceeds circle bounds, truncate line to circle.
        # Note: should only be run if EXACTLY ONE endpoint inside the circle!
        intersections = self.find_line_circle_intersections(start, end)
        start_inside = self.is_point_inside_circle(start)
        endpoint_inside = start if start_inside else end
        return [endpoint_inside, intersections[0]]
    
    # def get_spherical_midpoint(self, star1, star2):
    #     # Use modulo to handle wrap-around by shifting the midpoint calculation
    #     # Note, this could error if RA coords are diametrically opposed (=12h apart)
    #     midpoint_ra = (star1.ra + (star2.ra - star1.ra + 24) % 24 / 2) % 24
    #     midpoint_dec = (star1.dec + star2.dec) / 2
    #     return (midpoint_ra, midpoint_dec)

    def get_unit_sphere_cartesian(self, star):
        ra_rad = star.ra * (2 * np.pi / 24)
        # Convert DEC from degrees to radians
        dec_rad = star.dec * (np.pi / 180)
        # Calculate Cartesian coordinates
        x = np.cos(dec_rad) * np.cos(ra_rad)
        y = np.cos(dec_rad) * np.sin(ra_rad)
        z = np.sin(dec_rad)
        return (x, y, z)
    
    def slerp(self, star1, star2, num_points=5):
        # Get unit sphere Cartesian coordinates for both stars
        p1 = np.array(self.get_unit_sphere_cartesian(star1))
        p2 = np.array(self.get_unit_sphere_cartesian(star2))
        # Compute the cosine of the angular distance between the points
        cos_theta = np.dot(p1, p2)
        theta = np.arccos(np.clip(cos_theta, -1, 1))  # Clip to avoid precision issues causing NaNs
        # Slerp formula implementation
        results = []
        for i in range(1, num_points + 1):
            t = i / (num_points + 1)  # Fractional positions along the path
            sin_theta = np.sin(theta)
            if sin_theta == 0:
                # Both points are the same
                interpolated = p1
            else:
                a = np.sin((1 - t) * theta) / sin_theta
                b = np.sin(t * theta) / sin_theta
                interpolated = a * p1 + b * p2
            results.append(interpolated)
        return results
    
    def unit_sphere_cartesian_to_ra_dec(self, cartesian):
        x, y, z = cartesian
        # Calculate declination in radians and convert to degrees
        dec = np.arcsin(z) * (180 / np.pi)
        # Calculate right ascension in radians and convert to hours
        ra = np.arctan2(y, x) * (12 / np.pi)
        if ra < 0:
            ra += 24  # Normalize RA to be in the range [0, 24) hours
        return ra, dec
    
    # def add_smooth_bezier_path(self, points, stroke_width=0.2):
    #     print("\nadding smoothed bezier!! with points: \n")
    #     for point in points:
    #         print(point)
    #     path = svgwrite.path.Path(d=self.make_smooth_path_d(points), stroke="white", fill="none", stroke_width=stroke_width)
    #     self.drawing.add(path)

    # def make_smooth_path_d(self, points):
    #     n = len(points)
    #     if n < 2:
    #         return ""
    #     coords = [self.get_hemisphere_cartesian(p.ra, p.dec) for p in points]
    #     # Calculate tangents: Using a symmetric difference approach
    #     tangents = [0.5 * (np.array(coords[i+1]) - np.array(coords[i-1])) if 0 < i < n-1 else np.array(coords[1]) - np.array(coords[0]) if i == 0 else np.array(coords[-1]) - np.array(coords[-2]) for i in range(n)]
    #     # Modify tangent length based on desired curvature smoothness
    #     tangent_lengths = [np.linalg.norm(t) / 4 for t in tangents]  # adjust divisor for tightness of curve
    #     d = f"M {coords[0][0]}, {coords[0][1]}"
    #     for i in range(n-1):
    #         # Adjust control points to ensure smooth transitions
    #         ctrl1 = np.array(coords[i]) + tangents[i] * tangent_lengths[i]
    #         ctrl2 = np.array(coords[i+1]) - tangents[i+1] * tangent_lengths[i+1]
    #         d += f" C {ctrl1[0]}, {ctrl1[1]}, {ctrl2[0]}, {ctrl2[1]}, {coords[i+1][0]}, {coords[i+1][1]}"
    #     return d

    # def make_smooth_path_d(self, points):
    #     if len(points) < 2:
    #         return  # Need at least two points to draw a curve
    #     # Convert RA/DEC to Cartesian coordinates
    #     coords = [self.get_hemisphere_cartesian(p.ra, p.dec) for p in points]
    #     # Start the path at the first point
    #     d = f"M {coords[0][0]} {coords[0][1]} "
    #     # Use 'T' commands to create a smooth curve through all the points
    #     for coord in coords[1:]:
    #         d += f"T {coord[0]} {coord[1]} "
    #     # Add the path to the SVG
    #     path = svgwrite.path.Path(d=d, fill="none", stroke="white", stroke_width=0.2)
    #     self.drawing.add(path)

    def make_smooth_path_d(self, points, stroke_width):
        n = len(points)
        if n < 2:
            return ""
        coords = [np.array(self.get_hemisphere_cartesian(p.ra, p.dec)) for p in points]

        # Start the path at the first coordinate
        d = f"M {coords[0][0]} {coords[0][1]}"

        # Generate control points and segments
        for i in range(1, n):
            # Calculate control points one third and two thirds along the segment
            ctrl1 = coords[i-1] + (coords[i] - coords[i-1]) / 3
            ctrl2 = coords[i-1] + 2 * (coords[i] - coords[i-1]) / 3
            d += f" C {ctrl1[0]} {ctrl1[1]} {ctrl2[0]} {ctrl2[1]} {coords[i][0]} {coords[i][1]}"
        path = svgwrite.path.Path(d=d, fill="none", stroke="white", stroke_width=stroke_width)
        self.drawing.add(path)


    
    def transform_paths(self, paths_data, x_dim, y_dim):
        """Given some paths data, uses information from SVGHemisphere to transform into azimuthal equidistant cartesian coords."""
        transformed_paths = []
        for path_data in paths_data:
            transformed_path = {'M': None, 'beziers': []}
            # Transform the 'M' starting point
            m_ra, m_dec = get_equirect_coords(path_data['M'][0], path_data['M'][1], x_dim=x_dim, y_dim=y_dim)
            transformed_path['M'] = self.get_hemisphere_cartesian(m_ra, m_dec)
            # Transform each bezier segment
            for bezier in path_data['beziers']:
                transformed_bezier = []
                for point in bezier:
                    ra, dec = get_equirect_coords(point[0], point[1], x_dim=x_dim, y_dim=y_dim)
                    transformed_bezier.append(self.get_hemisphere_cartesian(ra, dec))
                transformed_path['beziers'].append(transformed_bezier)
            transformed_paths.append(transformed_path)
        return transformed_paths
    
    # SVG-specific coloration helper functions
    def create_star_gradient(self, star):
        # Create a radial gradient from white at the center to the specified outer color at the edges
        gradient = self.drawing.defs.add(self.drawing.radialGradient(id=star.hip))
        inner_color = '#FFFFFF'
        outer_color = bv_to_color(star.ci)
        gradient.add_stop_color(offset='10%', color=inner_color, opacity='1')
        gradient.add_stop_color(offset='70%', color=outer_color, opacity='1')
        gradient.add_stop_color(offset='100%', color=outer_color, opacity='0')
        return gradient

    # Drawing Functions
    def add_star_circle(self):
        # This is the night-sky circle containing all stars and constellations.
        self.elements['star_circle'] = self.drawing.add(
            self.drawing.circle(
                center=(self.size / 2, self.size / 2),
                r=self.star_circle_dia / 2,
                fill="#112233",
                stroke='black'
            )
        )

    def add_milky_way_svg(self, source_svg, x_dim, y_dim):
        paths_data = extract_and_structure_paths(source_svg)
        transformed_data = self.transform_paths(paths_data, x_dim, y_dim)
        combined_path_string = ""
        for path_data in transformed_data:
            path_string = f"M {path_data['M'][0]},{path_data['M'][1]} "  # Move to the start point
            for bezier in path_data['beziers']:
                path_string += f"C {bezier[1][0]},{bezier[1][1]} {bezier[2][0]},{bezier[2][1]} {bezier[3][0]},{bezier[3][1]} "
            combined_path_string += path_string  # Combine into one path string
        # Create a single SVG path element with the combined path string
        self.elements[source_svg] = self.drawing.add(
            self.drawing.path(
                d=combined_path_string, 
                fill="#FFFFFF", 
                fill_opacity=0.1, 
                stroke="none", 
                stroke_width=0, 
                fill_rule="evenodd"
            )
        )

    # def add_constellation_lines(self, constellationship, truncation_amount=0.1, stroke_width=0.2):
    #     for constellation in constellationship.constellations: # constellation = [(star, star), (star, star), ...]
    #         for stars in constellation.seg_list: # stars = (star, star)
    #             start = self.get_hemisphere_cartesian(stars[0].ra, stars[0].dec)
    #             end = self.get_hemisphere_cartesian(stars[1].ra, stars[1].dec)
    #             start, end = self.truncate_line(start, end, truncation_amount)
    #             if self.is_point_inside_circle(start) and self.is_point_inside_circle(end): # Both inside
    #                 self.elements[f'{stars[0].hip}-{stars[1].hip}'] = self.drawing.add(
    #                     self.drawing.line(
    #                         start, end, stroke='white', stroke_width=stroke_width)
    #                 )
    #             elif self.is_point_inside_circle(start) or self.is_point_inside_circle(end): # One inside
    #                 endpoints = self.truncate_line_to_circle(start, end)
    #                 self.elements[f'{stars[0].hip}-{stars[1].hip}'] = self.drawing.add(
    #                     self.drawing.line(
    #                         endpoints[0], endpoints[1], stroke='white', stroke_width=stroke_width)
    #                 )

    def add_constellation_lines(self, constellationship, stroke_width=0.2):
        for constellation in constellationship.constellations: # constellation = [(star, star), (star, star), ...]
            print('\nconstellation:\n', constellation, '\n\n')
            for stars in constellation.seg_list: # stars = (star, star)
                interpolated = svg_north.slerp(stars[0], stars[1], num_points=5)
                converted = [svg_north.unit_sphere_cartesian_to_ra_dec(cartesian) for cartesian in interpolated]
                points_list = [stars[0]]
                print('\n\n',stars[0])
                for counter, coords in enumerate(converted):
                    star = Star(counter-8, str(counter-8), coords[0], coords[1], 2, -0.5)
                    print(star)
                    points_list.append(star)
                points_list.append(stars[1])
                print(stars[1])
                self.make_smooth_path_d(points_list, stroke_width)
                # self.add_smooth_bezier_path(points_list, stroke_width=stroke_width)

    def add_stars(self, stars_dict, mag_limit, min_radius=0.5, max_radius=5, scale_type=1):
        for star in stars_dict.values():
            if self.is_north:
                if star.dec < 90 - self.dec_degrees or star.mag > mag_limit:
                    continue
            else:
                if star.dec > -90 + self.dec_degrees or star.mag > mag_limit:
                    continue
            x, y = self.get_hemisphere_cartesian(star.ra, star.dec)
            star_radius = mag_to_radius(star.mag, min_radius=min_radius, max_radius=max_radius, scale_type=scale_type, max_mag=mag_limit, min_mag=-1.46)
            self.create_star_gradient(star)
            self.elements[star.hip] = self.drawing.add(self.drawing.circle(center=(x, y), r=star_radius, fill=f'url(#{star.hip})'))

    def save_drawing(self):
        self.drawing.save()


# Other helper functions
def get_stars_dict(star_data_loc, data_types, mag_limit):
    """Generate a dictionary for stars with a HIP number within a magnitude limit."""
    star_data = pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)
    filtered_data = star_data.dropna(subset=['hip'])
    filtered_data = filtered_data[filtered_data['mag'] <= mag_limit]
    return {int(row['hip']): Star(int(row['hip']), row['proper'], row['ra'], row['dec'], row['mag'], row['ci'])
        for row in filtered_data.to_dict('records')}

def mag_to_radius(mag, min_radius=1, max_radius=10, scale_type=1, max_mag=6.5, min_mag=0):
    if mag > max_mag:
        return min_radius  # If the star is dimmer than the max visible magnitude, use the smallest radius
    if mag < min_mag:
        mag = min_mag  # Cap the magnitude at the brightest star's magnitude to avoid negative radius values
    relative_magnitude = (max_mag - mag) / (max_mag - min_mag)
    scaled_magnitude = relative_magnitude ** scale_type
    radius = min_radius + (max_radius - min_radius) * scaled_magnitude
    return radius

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

def extract_and_structure_paths(svg_file):
    tree = ET.parse(svg_file)
    root = tree.getroot()
    all_paths = []
    for path in root.findall('{http://www.w3.org/2000/svg}path'):
        d_attr = re.sub(r'\s*[Zz]\s*', '', path.get('d').strip())
        subpaths = []
        current_path = None
        last_point = None
        commands = re.findall(r'([MC][^MC]*)', d_attr)
        for command in commands:
            command_type = command[0]
            points = [(float(x), float(y)) for x, y in re.findall(r'(\d+\.?\d*),(\d+\.?\d*)', command[1:].strip())]
            if command_type == 'M':
                current_path = {'M': points[0], 'beziers': []}
                subpaths.append(current_path)
                last_point = points[0]
            elif command_type == 'C' and current_path is not None:
                for i in range(0, len(points), 3):
                    if i + 2 < len(points) and last_point is not None:
                        bezier = [last_point, points[i], points[i+1], points[i+2]]
                        current_path['beziers'].append(bezier)
                        last_point = points[i+2]
        all_paths.extend(subpaths)
    return all_paths

def get_equirect_coords(x, y, x_dim, y_dim):
    """Converts equirectangular pixel coordinates to RA and DEC. x_dim should be width in pixels, y_dim should be height in pixels."""
    ra = 24 - 24 * (x / x_dim)
    dec = 90 - 180 * (y / y_dim)
    return ra, dec

if __name__ == '__main__':
    # Constants
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

    # SVG Information
    size = 2000 # This is the size of the ENTIRE ILLUSTRATION. 
    full_circle_dia = 1800 # The circle containing the starscape AND months, tickmarks.
    star_circle_dia = 1600 # This is the circle containing our stars
    dec_degrees = 108 # Number of degrees of declination to map. 90 would be one celestial hemisphere. 180 is both hemispheres, but it gets hella distorted!

    # Milky Way SVG Information
    svg_files = ['star_map/data/mw_1.svg', 'star_map/data/mw_2.svg', 'star_map/data/mw_3.svg']
    x_dim = 8998
    y_dim = 4498

    # Star and Constellation Information
    stars_dict = get_stars_dict(STAR_DATA_LOC, DATA_TYPES, 100) # No mag limit for now.
    parser = ConstellationParser(CONSTELLATIONS_LOC, stars_dict)
    constellations = parser.parse()
    constellationship = Constellationship(constellations, 'iau')
    
    # Drawing Curved Paths
    svg_north = SVGHemisphere(size, full_circle_dia, star_circle_dia, dec_degrees, filename="star_map.svg", is_north=True)
    svg_north.add_star_circle()
    for file in svg_files:
        svg_north.add_milky_way_svg(file, x_dim, y_dim)
    svg_north.add_constellation_lines(constellationship, stroke_width=0.5)
    svg_north.add_stars(stars_dict, mag_limit=7.5, min_radius=0.3, max_radius=11, scale_type=1.5)
    svg_north.save_drawing()

    # straight paths
    # svg_north.add_star_circle()
    # for file in svg_files:
    #     svg_north.add_milky_way_svg(file, x_dim, y_dim)
    # svg_north.add_constellation_lines(constellationship, truncation_amount=2.5, stroke_width=0.2)
    # svg_north.add_stars(stars_dict, mag_limit=7.5, min_radius=0.1, max_radius=5, scale_type=1.5)
    # svg_north.save_drawing()

    # straight paths south
    # svg_south = SVGHemisphere(size, full_circle_dia, star_circle_dia, dec_degrees, filename="star_map_s.svg", is_north=False)
    # svg_south.add_star_circle()
    # for file in svg_files:
    #     svg_south.add_milky_way_svg(file, x_dim, y_dim)
    # svg_south.add_constellation_lines(constellationship, truncation_amount=2.5, stroke_width=0.2)
    # svg_south.add_stars(stars_dict, mag_limit=7.5, min_radius=0.1, max_radius=5, scale_type=1.3)
    # svg_south.save_drawing()