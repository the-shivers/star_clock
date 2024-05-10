import os
import pandas as pd
from star_map_classes import (
    StarHolder,
    Constellationship,
    ConstellationParser,
    SVGHemisphere
)

# Project directory
dir = os.path.expanduser('~/star_clock')

# Star and constellation config
sky_culture = 'snt'
star_data_loc = f'{dir}/star_map/data/stars/athyg_24_reduced_m10.csv'
constellations_json = f'{dir}/star_map/data/sky_cultures/{sky_culture}/constellationship.json'
const_coords_loc = f'{dir}/star_map/data/sky_cultures/{sky_culture}/constellation_coords.csv'
is_north = True
mag_limit = 5.5 # For limiting size of stars list.
min_radius = 0.5
max_radius = 10
scale_type = 1.3 # >1 = exponential with bigger stars MUCH brighter than dim ones, 1 = linear, <1 = logarithmic
gradient = False # Gradient effect on stars

# SVG configuration
size = 2000 # This is the size of the ENTIRE ILLUSTRATION. 
full_circle_dia = 1780 # The circle containing the starscape AND months, tickmarks.
star_circle_dia = 1600 # This is the circle containing our stars
dec_degrees = 108 # Number of degrees of declination to map. 90 would be one celestial hemisphere. 180 is both hemispheres, but it gets hella distorted!
hemisphere_code = 'n' if is_north else 's'
output_loc = f'{dir}/star_map/svg_output/star_map_{hemisphere_code}.svg'

# Milky Way SVG config
svg_files = [f'{dir}/star_map/data/milky_way/mw_1.svg', f'{dir}/star_map/data/milky_way/mw_2.svg', f'{dir}/star_map/data/milky_way/mw_3.svg']
x_dim = 8998
y_dim = 4498

# Color palette and fonts styles
star_circle_col = '#1E3A56'
star_circle_stroke = '#FFFFFF'

date_circle_col = '#100F20'
date_circle_stroke = 'none'

tick_col = '#ffffff'
stroke_width = 1
maj_tick = 15
min_tick = 6
sev_tick = 50

axes_col = '#3C6893'
axes_n = 24
axes_stroke_width = 1
axes_ticks = True
axes_tick_degs = 10
axes_tick_width = 6

milky_way_color = '#FFFFFF'
milky_way_alpha = 0.08

equator_col = '#FFFFFF'
equator_stroke_width = 0.3

constellation_lines_col = '#86C2FF'
constellation_stroke_width = 1
constellation_trunc_rate = 1.5

styles = {
    'constellation': {
        'font_family': 'Josefin Sans',
        'font_size': 12,
        'font_weight': 300, # Light
        'font_style': 'normal', # 'normal', 'italic
        'letter_spacing': 5, # 'normal', 5, 3
        'fill': '#86C2FF',
        'src': f'{dir}/star_map/fonts/JosefinSans-Light.ttf',
        'stroke': star_circle_col,
        'stroke_width': '5px'
    },
    'small': {
        'font_family': 'Josefin Sans',
        'font_size': 14,
        'font_weight': 300, # Light
        'font_style': 'normal',
        'letter_spacing': 1,
        'font_style': 'normal',
        'fill': '#FFFFFF',
        'src': f'{dir}/star_map/fonts/JosefinSans-Light.ttf'
    },
    'month': {
        'font_family': 'Josefin Sans',
        'font_size': 20,
        'font_weight': 500, # Medium
        'font_style': 'normal', # 'normal', 'italic
        'letter_spacing': 7, # 'normal', 5, 3
        'fill': '#86C2FF',
        'src': f'{dir}/star_map/fonts/JosefinSans-Medium.ttf',
    },
    'ra_dec_lab': {
        'font_family': 'Josefin Sans',
        'font_size': 10,
        'font_weight': 300, # Light
        'font_style': 'normal', # 'normal', 'italic
        'letter_spacing': 1, # 'normal', 5, 3
        'fill': '#3C6893',
        'src': f'{dir}/star_map/fonts/JosefinSans-Light.ttf',
    }
}

if __name__ == '__main__':
    print("Getting stars...")
    starholder = StarHolder(star_data_loc)
    parser = ConstellationParser(constellations_json, starholder)
    print("Parsing constellations...")
    constellations = parser.parse()
    print("Building constellationship...")
    constellationship = Constellationship(constellations, sky_culture)
    print("Getting coordinates of constellations...")
    constellation_coords_df = pd.read_csv(const_coords_loc)
    constellation_coords_dict = constellation_coords_df.set_index('latin_name').to_dict(orient='index')
    
    print("Building SVG...")
    svg_north = SVGHemisphere(size, full_circle_dia, star_circle_dia, dec_degrees, filename=output_loc, is_north=is_north)
    svg_north.add_date_circle(fill=date_circle_col, stroke=date_circle_stroke)
    svg_north.add_dates_and_tickmarks(style=styles['small'], tick_color=tick_col, stroke_width=stroke_width, maj_tick=maj_tick, min_tick=min_tick)
    svg_north.add_star_circle(fill=star_circle_col, stroke=star_circle_stroke)
    svg_north.add_azimuthal_axes(style=styles['ra_dec_lab'], n=axes_n, stroke_color=axes_col, stroke_width=axes_stroke_width, tick_degs=axes_tick_degs, tick_width=axes_tick_width)
    svg_north.add_months(styles['month'])
    
    print("Adding constellation lines...")
    svg_north.add_constellation_lines_curved(constellationship, stroke_width=constellation_stroke_width, stroke_color=constellation_lines_col)
    print("Adding stars...")
    svg_north.add_stars(starholder, constellationship, mag_limit=mag_limit, min_radius=min_radius, max_radius=max_radius, scale_type=scale_type, gradient=gradient, mask_id = 'starfield-mask')
    print("Masking...")
    svg_north.add_segment_mask(constellationship, truncation_rate=constellation_trunc_rate, mask_id='segment-masks')
    svg_north.add_star_and_mw_mask(mask_id='starfield-mask')
    print("Labelling...")
    for key, value in constellation_coords_dict.items():
        for i, subkey in enumerate(key.split(' ')):
            mult = 1 if is_north else -1
            svg_north.add_text_centered_rotated(subkey.upper(), styles['constellation'], value['ra'], value['dec'] + i * 2 * mult, mult * value['ra'] / 24 * 360 + value['rot'])
    print("Milking...")
    for file in svg_files:
        svg_north.add_milky_way_svg(file, x_dim, y_dim, color=milky_way_color, opacity=milky_way_alpha)
    print(f"Saving file! {output_loc}")
    svg_north.save_drawing()

# Core goals
# TODO: DSO symbols and labels
# TODO: Star labels
# TODO: Adjust label positions
# TODO: ecliptic maybe?


# Optioanl goals
# TODO: Figure out why constellation parsing is so damn slow, possibly pickle them.
# TODO: Fix radial gradient performance. We can just do classes for these and assign them more reasonably than having 10000. Especially relevant when we get to glow if we want to add it.
# TODO: Stars: Gradient should be a parameter. Should probably be defined with a dict up front or list or something.
# TODO: Stars: Glow.
# TODO: Maybe glowing lines.
# TODO: Color scheme updates, config dicts, palettes


