import xml.etree.ElementTree as ET
import re
import numpy as np

def get_northern_hemisphere_cartesian2(ra, dec, graphic_size=1000, dec_degrees=108, star_circle_dia=800):
    theta_radians = ra / 24 * np.pi * 2
    hypotenuse = ((90 - dec) / dec_degrees) * star_circle_dia/2
    x = (graphic_size - star_circle_dia) / 2 + star_circle_dia/2 + np.sin(theta_radians) * hypotenuse
    y = graphic_size/2 - np.cos(theta_radians) * hypotenuse
    return (x, y)

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

def get_equirect_coords(x, y):
    """Converts equirectangular pixel coordinates to RA and DEC"""
    ra = 24 - 24 * (x / 8998)
    dec = 90 - 180 * (y / 4498)
    return ra, dec

def transform_paths(paths_data):
    transformed_paths = []
    for path_data in paths_data:
        transformed_path = {'M': None, 'beziers': []}
        # Transform the 'M' starting point
        m_ra, m_dec = get_equirect_coords(path_data['M'][0], path_data['M'][1])
        transformed_path['M'] = get_northern_hemisphere_cartesian2(m_ra, m_dec)
        # Transform each bezier segment
        for bezier in path_data['beziers']:
            transformed_bezier = []
            for point in bezier:
                ra, dec = get_equirect_coords(point[0], point[1])
                transformed_bezier.append(get_northern_hemisphere_cartesian2(ra, dec))
            transformed_path['beziers'].append(transformed_bezier)
        transformed_paths.append(transformed_path)
    return transformed_paths

svg_file = 'star_map/data/mw_1.svg'
paths_data = extract_and_structure_paths(svg_file)
transformed_data = transform_paths(paths_data)

svg_file = 'star_map/data/mw_2.svg'
paths_data = extract_and_structure_paths(svg_file)
transformed_data = transform_paths(paths_data)

svg_file = 'star_map/data/mw_3.svg'
paths_data = extract_and_structure_paths(svg_file)
transformed_data = transform_paths(paths_data)


