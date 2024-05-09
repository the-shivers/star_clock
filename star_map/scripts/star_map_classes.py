import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from PIL import ImageFont, ImageDraw, Image
import re
import svgwrite
import json
import operator


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
    

class StarHolder:
    def __init__(self, csv_loc, mag_limit=20):
        self.stars = self.get_stars(csv_loc, mag_limit)

    def get_stars(self, csv_loc, mag_limit):
        """Generate a dictionary for stars with a HIP number within a magnitude limit."""
        dtypes = {
            'hip': 'Int64', # Nullable integer, for stars with no HIPPARCOS ID
            'proper': str,
            'ra': float,
            'dec': float,
            'mag': float,
            'ci': float
        }
        star_data = pd.read_csv(csv_loc, usecols=dtypes.keys(), dtype=dtypes)
        filtered_data = star_data.dropna(subset=['hip'])
        filtered_data = filtered_data[filtered_data['mag'] <= mag_limit]
        filtered_data.sort_values('mag', inplace=True, ascending=False)
        return [Star(int(row['hip']), row['proper'], row['ra'], row['dec'], row['mag'], row['ci'])
            for row in filtered_data.to_dict('records')]

    def get_star_by_name(self, name):
        for star in self.stars:
            if star.proper == name:
                return star

    def get_star_by_hip(self, hip):
        for star in self.stars:
            if star.hip == hip:
                return star

    def get_filtered_stars(self, filters):
        """
        filters should be in the format:
            [('numerical_attribute', 'comparison', value), ...]
        So, for example:
            [('mag', '<', 1), ('mag', '>=', 0), ('ci', '<', -0.2)] # Find bluish stars between mag 0 and 1.
        """
        comps = {'<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge, '==': operator.eq, '!=': operator.ne}
        filtered_stars = []
        for star in self.stars:
            match = True
            for attr, comp, value in filters:
                if not comps[comp](getattr(star, attr), value):
                    match = False
                    break
            if match:
                filtered_stars.append(star)
        return filtered_stars
    
    def head(self, n):
        for star in self.stars[0:n]:
            print(star)

    def tail(self, n):
        for star in self.stars[-n:]:
            print(star)
    
    def __repr__(self):
        return f"StarHolder with {len(self.stars)} stars"


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
    def __init__(self, filepath, starholder):
        self.filepath = filepath
        self.starholder = starholder

    def parse(self):
        with open(self.filepath, 'r') as file:
            data = json.load(file)['constellations']
        constellations = []
        for item in data:
            abbrv = item['id'].split()[-1]
            if len(abbrv) >= 4: # Indicates a subconstellation, e.g. ursa major paws for rey star culture
                continue
            common_name = item['common_name']['english']
            latin_name = item['common_name'].get('native', common_name)
            seg_list = []
            for line in item['lines']:
                for i in range(len(line) - 1):
                    # Error correction
                    # if self.starholder.get_star_by_hip(line[i]) is None or self.starholder.get_star_by_hip(line[i + 1] is None):
                    #     print(f"Error: {item['common_name']['english']} could not find stars: {line[i]}, {line[i + 1]}")
                    seg_list.append((self.starholder.get_star_by_hip(line[i]), self.starholder.get_star_by_hip(line[i + 1])))
            seg_count = len(seg_list)
            constellations.append(Constellation(abbrv, common_name, latin_name, seg_count, seg_list))
        return constellations
    
    def __repr__(self):
        return (
            f"ConstellationParser(filepath={self.filepath!r}, starholder={self.starholder})"
        )
    

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
    
    def __repr__(self):
        return (
            f'Constellationship(constellations={self.constellations[0:3]}, name={self.name})'
        )
    
    def __str__(self):
        return (
            f'Constellationship "{self.name}" with {len(self.constellations)} constellations.'
        )


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
    
    def is_point_inside_star_circle(self, point):
        x, y = point
        x0, y0 = self.size/2, self.size/2
        return (x - x0)**2 + (y - y0)**2 <= (self.star_circle_dia/2)**2

    def get_unit_sphere_cartesian(self, star):
        ra_rad = star.ra * (2 * np.pi / 24)
        dec_rad = star.dec * (np.pi / 180)
        x = np.cos(dec_rad) * np.cos(ra_rad)
        y = np.cos(dec_rad) * np.sin(ra_rad)
        z = np.sin(dec_rad)
        return (x, y, z)
    
    def slerp(self, star1, star2, num_points=5):
        p1 = np.array(self.get_unit_sphere_cartesian(star1))
        p2 = np.array(self.get_unit_sphere_cartesian(star2))
        cos_theta = np.dot(p1, p2)
        theta = np.arccos(np.clip(cos_theta, -1, 1))  # Clip to avoid precision issues causing NaNs
        results = []
        for i in range(1, num_points + 1):
            t = i / (num_points + 1)
            sin_theta = np.sin(theta)
            if sin_theta == 0:
                interpolated = p1
            else:
                a = np.sin((1 - t) * theta) / sin_theta
                b = np.sin(t * theta) / sin_theta
                interpolated = a * p1 + b * p2
            results.append(interpolated)
        return results
    
    def unit_sphere_cartesian_to_ra_dec(self, cartesian):
        x, y, z = cartesian
        dec = np.arcsin(z) * (180 / np.pi)
        ra = np.arctan2(y, x) * (12 / np.pi)
        if ra < 0:
            ra += 24
        return ra, dec

    def make_smooth_path_d(self, points, stroke_color, stroke_width, mask_id='star-masks'):
        n = len(points)
        if n < 2:
            return ""
        coords = [np.array(self.get_hemisphere_cartesian(p.ra, p.dec)) for p in points]
        d = f"M {coords[0][0]} {coords[0][1]}"
        for i in range(1, n):
            # Calculate control points one third and two thirds along the segment
            ctrl1 = coords[i-1] + (coords[i] - coords[i-1]) / 3
            ctrl2 = coords[i-1] + 2 * (coords[i] - coords[i-1]) / 3
            d += f" C {ctrl1[0]} {ctrl1[1]} {ctrl2[0]} {ctrl2[1]} {coords[i][0]} {coords[i][1]}"
        path = svgwrite.path.Path(d=d, fill="none", stroke=stroke_color, stroke_width=stroke_width, mask=f'url(#{mask_id})')
        self.elements[f'{points[0].hip}-{points[-1].hip}'] = self.drawing.add(path)

    def get_equirect_coords(self, x, y, x_dim, y_dim):
        """Converts equirectangular pixel coordinates to RA and DEC. x_dim should be width in pixels, y_dim should be height in pixels."""
        ra = 24 - 24 * (x / x_dim)
        dec = 90 - 180 * (y / y_dim)
        return ra, dec

    def transform_paths(self, paths_data, x_dim, y_dim):
        """Given some paths data, uses information from SVGHemisphere to transform into azimuthal equidistant cartesian coords."""
        transformed_paths = []
        for path_data in paths_data:
            transformed_path = {'M': None, 'beziers': []}
            m_ra, m_dec = self.get_equirect_coords(path_data['M'][0], path_data['M'][1], x_dim=x_dim, y_dim=y_dim)
            transformed_path['M'] = self.get_hemisphere_cartesian(m_ra, m_dec)
            for bezier in path_data['beziers']:
                transformed_bezier = []
                for point in bezier:
                    ra, dec = self.get_equirect_coords(point[0], point[1], x_dim=x_dim, y_dim=y_dim)
                    transformed_bezier.append(self.get_hemisphere_cartesian(ra, dec))
                transformed_path['beziers'].append(transformed_bezier)
            transformed_paths.append(transformed_path)
        return transformed_paths
    
    # SVG-specific coloration/sizing helper functions
    def bv_to_color(self, bv):
        colors = ['#94B6FF', '#99B9FF', '#C9D9FF', '#ECEEFF', '#FFFFFF', '#FFF2EE', '#FFE7D2', '#FFCC98']
        breakpoints = [-0.33, -0.15, 0.0, 0.25, 0.58, 0.81, 1.4, float('inf')]  # Open-ended for values beyond 1.4
        for i in range(len(breakpoints) - 1):
            if breakpoints[i] <= bv < breakpoints[i + 1]:
                lower_bound = breakpoints[i]
                upper_bound = breakpoints[i + 1]
                t = (bv - lower_bound) / (upper_bound - lower_bound) if upper_bound != float('inf') else 1.0
                color_start = tuple(int(colors[i][j:j+2], 16) for j in range(1, 7, 2))
                if i+1 < len(colors):
                    color_end = tuple(int(colors[i+1][j:j+2], 16) for j in range(1, 7, 2))
                else:
                    color_end = color_start
                r = int(color_start[0] + t * (color_end[0] - color_start[0]))
                g = int(color_start[1] + t * (color_end[1] - color_start[1]))
                b = int(color_start[2] + t * (color_end[2] - color_start[2]))
                return f"#{r:02x}{g:02x}{b:02x}"
        return '#ffffff'

    def create_star_gradient(self, star):
        gradient = self.drawing.defs.add(self.drawing.radialGradient(id=star.hip))
        inner_color = '#FFFFFF'
        outer_color = self.bv_to_color(star.ci)
        gradient.add_stop_color(offset='10%', color=inner_color, opacity='1')
        gradient.add_stop_color(offset='70%', color=outer_color, opacity='1')
        gradient.add_stop_color(offset='100%', color=outer_color, opacity='0')
        return gradient
    
    def mag_to_radius(self, mag, min_radius=1, max_radius=10, scale_type=1, max_mag=6.5, min_mag=0):
        if mag > max_mag:
            return min_radius  # If the star is dimmer than the max visible magnitude, use the smallest radius
        if mag < min_mag:
            mag = min_mag  # Cap the magnitude at the brightest star's magnitude to avoid negative radius values
        relative_magnitude = (max_mag - mag) / (max_mag - min_mag)
        scaled_magnitude = relative_magnitude ** scale_type
        radius = min_radius + (max_radius - min_radius) * scaled_magnitude
        return radius
    
    def add_star_mask(self, constellationship, truncation_rate=2, mask_id = 'star-masks'):
        # Note: must be performed AFTER stars have been added to the drawing.
        self.elements[mask_id] = self.drawing.defs.add(self.drawing.mask(id=mask_id, maskUnits="userSpaceOnUse"))
        self.elements[mask_id].add(self.drawing.circle((self.size/2, self.size/2), r=self.star_circle_dia/2, fill="white"))
        constellationship.get_unique_stars()
        for star in constellationship._unique_stars:
            if (self.is_north and star.dec > 90 - self.dec_degrees) or (not self.is_north and star.dec < self.dec_degrees - 90):
                x = self.elements[star.hip].attribs['cx']
                y = self.elements[star.hip].attribs['cy']
                star_radius = self.elements[star.hip].attribs['r']
                self.elements[mask_id].add(self.drawing.circle(center=(x, y), r=star_radius + truncation_rate, fill="black"))
    
    # Other helper functions
    def extract_and_structure_paths(self, svg_file):
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

    # Drawing Functions
    def add_star_circle(self, fill="#112233", stroke="black"):
        self.elements['star_circle'] = self.drawing.add(
            self.drawing.circle(
                center=(self.size / 2, self.size / 2),
                r=self.star_circle_dia / 2,
                fill=fill,
                stroke=stroke
            )
        )

    def add_azimuthal_axes(self, n=8, stroke_color='#FFFFFF', stroke_width=1, ticks=True, tick_degs=10, tick_width=10, dec_rings=True):
        """Adds azimuthal axes. n is the number of slices we divide the pizza into."""
        ticks = int((self.dec_degrees - 0.1) // tick_degs)
        max_tick_distance = tick_degs * (ticks + 1) / self.dec_degrees * self.star_circle_dia / 2
        tick_distances = [j / (ticks + 1) * max_tick_distance for j in range(1, ticks+1)]
        if self.is_north:
            labels = [f'{90 - int(tick_degs*(j+1))}°' for j in range(len(tick_distances))]
        else:
            labels = [f'{-90 + int(tick_degs*(j+1))}°' for j in range(len(tick_distances))]
        tick_dict = dict(zip(labels, tick_distances))
        for i in range(n):
            unit_x = np.sin(2*i/n*np.pi)
            unit_y = np.cos(2*i/n*np.pi)
            self.elements[f'axis_{i}'] = self.drawing.add(
                self.drawing.line(
                    start=(self.size/2, self.size/2),
                    end=(self.size/2+self.star_circle_dia/2*unit_x, self.size/2+self.star_circle_dia/2*unit_y),
                    stroke_width=stroke_width,
                    stroke=stroke_color
                )
            )
            if ticks:
                for key, value in tick_dict.items():
                    start_x = self.size/2+value*unit_x
                    start_y = self.size/2+value*unit_y
                    x1 = start_x + unit_y * tick_width/2
                    y1 = start_y - unit_x * tick_width/2
                    x2 = start_x - unit_y * tick_width/2
                    y2 = start_y + unit_x * tick_width/2
                    self.elements[f'axis_{i}_tick_{key[:-1]}'] = self.drawing.add(
                        self.drawing.line(
                            start=(x1, y1),
                            end=(x2, y2),
                            stroke_width=stroke_width,
                            stroke=stroke_color
                        )
                    )
            if dec_rings:
                for key, value in tick_dict.items():
                    self.elements[f'dec_ring_{key[:-1]}'] = self.drawing.add(
                        self.drawing.circle(
                            center=(self.size/2, self.size/2), r=value, fill='none', stroke=stroke_color, stroke_width=stroke_width
                        )
                    )

    def add_equator(self, stroke_color='#FFFFFF', stroke_width=0.2):
        if self.dec_degrees < 90:
            return
        radius = 90 / self.dec_degrees * self.star_circle_dia / 2
        self.elements['equator'] = self.drawing.add(
            self.drawing.circle(
                center=(self.size/2, self.size/2), r=radius, fill='none', stroke=stroke_color, stroke_width=stroke_width
            )
        )

    def add_milky_way_svg(self, source_svg, x_dim, y_dim, color='#FFFFFF', opacity=0.1):
        paths_data = self.extract_and_structure_paths(source_svg)
        transformed_data = self.transform_paths(paths_data, x_dim, y_dim)
        combined_path_string = ""
        for path_data in transformed_data:
            path_string = f"M {path_data['M'][0]},{path_data['M'][1]} "
            for bezier in path_data['beziers']:
                path_string += f"C {bezier[1][0]},{bezier[1][1]} {bezier[2][0]},{bezier[2][1]} {bezier[3][0]},{bezier[3][1]} "
            combined_path_string += path_string
        self.elements[source_svg] = self.drawing.add(
            self.drawing.path(
                d=combined_path_string, 
                fill=color, 
                fill_opacity=opacity, 
                stroke="none", 
                stroke_width=0, 
                fill_rule="evenodd"
            )
        )

    def add_constellation_lines_curved(self, constellationship, stroke_color='#FFFFFF', stroke_width=0.2, mask_id='star-masks'):
        for constellation in constellationship.constellations: # constellation = [(star, star), (star, star), ...]
            for stars in constellation.seg_list: # stars = (star, star)
                start = self.get_hemisphere_cartesian(stars[0].ra, stars[0].dec)
                end = self.get_hemisphere_cartesian(stars[1].ra, stars[1].dec)
                if not self.is_point_inside_star_circle(start) and not self.is_point_inside_star_circle(end):
                    continue
                interpolated = self.slerp(stars[0], stars[1], num_points=5)
                converted = [self.unit_sphere_cartesian_to_ra_dec(cartesian) for cartesian in interpolated]
                points_list = [stars[0]] + [Star(0, '', coords[0], coords[1], 2, -0.5) for coords in converted] + [stars[1]]
                self.make_smooth_path_d(points_list, stroke_color, stroke_width, mask_id)
        
    def add_constellation_lines_straight(self, constellationship, stroke_color='#FFFFFF', stroke_width=0.2, mask_id='star-masks'):
        for constellation in constellationship.constellations: # constellation = [(star, star), (star, star), ...]
            for stars in constellation.seg_list: # stars = (star, star)
                start = self.get_hemisphere_cartesian(stars[0].ra, stars[0].dec)
                end = self.get_hemisphere_cartesian(stars[1].ra, stars[1].dec)
                if not self.is_point_inside_star_circle(start) and not self.is_point_inside_star_circle(end):
                    continue
                self.elements[f'{stars[0].hip}-{stars[1].hip}'] = self.drawing.add(
                    self.drawing.line(
                        start, end, stroke=stroke_color, stroke_width=stroke_width, mask=f'url(#{mask_id})')
                )

    def add_stars(self, starholder, constellationship, mag_limit, min_radius=0.5, max_radius=5, scale_type=1, gradient=True, glow=False):
        constellation_stars = constellationship.get_unique_stars()
        for star in starholder.stars:
            if self.is_north:
                if star.dec < 90 - self.dec_degrees or (star.mag > mag_limit and star not in constellation_stars):
                    continue
            else:
                if star.dec > -90 + self.dec_degrees or (star.mag > mag_limit and star not in constellation_stars):
                    continue
            x, y = self.get_hemisphere_cartesian(star.ra, star.dec)
            star_radius = self.mag_to_radius(star.mag, min_radius=min_radius, max_radius=max_radius, scale_type=scale_type, max_mag=mag_limit, min_mag=-1.46)
            if gradient:
                self.create_star_gradient(star)
                fill = f'url(#{star.hip})'
            else:
                fill = self.bv_to_color(star.ci)
            self.elements[star.hip] = self.drawing.add(self.drawing.circle(center=(x, y), r=star_radius, fill=fill))

    def add_text_centered_rotated(self, text, style, ra, dec, rotation):
        if self.is_north and dec < 90 - self.dec_degrees or not self.is_north and dec > -90 + self.dec_degrees:
            return
        x, y = self.get_hemisphere_cartesian(ra, dec)
        font = ImageFont.truetype(style['src'], style['font_size'])
        dummy_image = Image.new('RGB', (1, 1))
        draw = ImageDraw.Draw(dummy_image)
        # Get width, correct for additional letter spacing
        text_width = 0
        if isinstance(style['letter_spacing'], int):
            for char in text:
                bbox = draw.textbbox((0, 0), char, font=font)
                char_width = bbox[2] - bbox[0]
                text_width += char_width + style['letter_spacing']
        text_width -= style['letter_spacing']  # Remove extra spacing added at the end
        # Now we can get height
        bbox = draw.textbbox((0, 0), text, font=font)
        text_height = bbox[3] - bbox[1]
        centered_x = x - text_width / 2
        centered_y = y + text_height / 2
        text_element_stroke = self.drawing.text(
            text,
            insert=(centered_x, centered_y),
            font_family=style['font_family'],
            font_size=style['font_size'],
            font_weight=style['font_weight'],
            font_style=style['font_style'],
            letter_spacing=style['letter_spacing'],
            transform=f'rotate({rotation}, {x}, {y})',
            fill='none',
            stroke=style['stroke'],
            stroke_width=style['stroke_width']
        )
        text_element = self.drawing.text(
            text,
            insert=(centered_x, centered_y),
            font_family=style['font_family'],
            font_size=style['font_size'],
            font_weight=style['font_weight'],
            font_style=style['font_style'],
            letter_spacing=style['letter_spacing'],
            transform=f'rotate({rotation}, {x}, {y})',
            fill=style['fill'],
            stroke='none'
        )
        self.elements[f'text_{text[0:10].replace(" ","")}_stroke'] = self.drawing.add(text_element_stroke)
        self.elements[f'text_{text[0:10].replace(" ","")}'] = self.drawing.add(text_element)

    def save_drawing(self):
        self.drawing.save()