import os
from star_map_classes import (
    StarHolder,
    Constellationship,
    ConstellationParser,
    SVGHemisphere
)

# Project directory
dir = os.path.expanduser('~/star_clock')

# Star and constellation config
constellations_json = f'{dir}/star_map/data/wester_iau_sky_culture.json'
mag_limit = 7.5 # For limiting size of stars list.
star_data_loc = f'{dir}/star_map/data/athyg_24_reduced_m10.csv'

# SVG configuration
size = 2000 # This is the size of the ENTIRE ILLUSTRATION. 
full_circle_dia = 1800 # The circle containing the starscape AND months, tickmarks.
star_circle_dia = 1600 # This is the circle containing our stars
dec_degrees = 108 # Number of degrees of declination to map. 90 would be one celestial hemisphere. 180 is both hemispheres, but it gets hella distorted!
output_loc = f'{dir}/star_map/svg_output/star_map.svg'

# Milky Way SVG config
svg_files = [f'{dir}/star_map/data/mw_1.svg', f'{dir}/star_map/data/mw_2.svg', f'{dir}/star_map/data/mw_3.svg']
x_dim = 8998
y_dim = 4498

# Label styles
styles = {
    'big': {
        'font_family': 'Josefin Sans',
        'font_size': 20,
        'font_weight': 400, # Regular
        'font_style': 'normal',
        'fill': '#FFFFFF',
        'src': f'{dir}/star_map/data/JosefinSans-Regular.ttf'
    },
    'small': {
        'font_family': 'Josefin Sans',
        'font_size': 14,
        'font_weight': 300, # Light
        'font_style': 'normal',
        'fill': '#FFFFFF',
        'src': f'{dir}/star_map/data/JosefinSans-Light.ttf'
    }
}

if __name__ == '__main__':
    # Star and Constellation Information
    starholder = StarHolder(star_data_loc, 7.5) # No mag limit for now.
    parser = ConstellationParser(constellations_json, starholder)
    constellations = parser.parse()
    constellationship = Constellationship(constellations, 'iau')
    
    svg_north = SVGHemisphere(size, full_circle_dia, star_circle_dia, dec_degrees, filename=output_loc, is_north=True)
    svg_north.add_star_circle(fill="#1E3A56", stroke="none")
    svg_north.add_azimuthal_axes(n=8, stroke_color='#3C6893', stroke_width=1, ticks=True, tick_degs=10, tick_width=8, dec_rings=False)
    for file in svg_files:
        svg_north.add_milky_way_svg(file, x_dim, y_dim, opacity=0.08)
    # svg_north.add_constellation_lines_straight(constellationship, stroke_width=1, stroke_color='#86C2FF')
    svg_north.add_constellation_lines_curved(constellationship, stroke_width=1, stroke_color='#86C2FF')
    svg_north.add_stars(starholder, mag_limit=6.5, min_radius=0.5, max_radius=10, scale_type=1.3, gradient=False)
    svg_north.add_star_mask(constellationship, truncation_rate=1.3, mask_id='star-masks')
    # svg_north.add_text('URSA MINOR', (940, 1060), fill='#FFFFFF', font_weight=400, font_size=20)
    # svg_north.add_text('Polaris', (1013, 1001), fill='#FFFFFF', font_weight=300, font_size=14)
    # svg_north.add_text_centered_rotated('O', styles['big'], 0, 90, 0)
    svg_north.add_text_centered_rotated('OOOOO', styles['big'], 0, 90, 0)
    svg_north.add_text_centered_rotated('OOOOO', styles['big'], 0, 90, 90)
    svg_north.add_text_centered_rotated('OOOOO', styles['big'], 0, 90, 180)
    svg_north.add_text_centered_rotated('OOOOO', styles['big'], 0, 90, 270)
    # svg_north.add_text_centered_rotated('O Rot', styles['big'], 0, 90, 180)
    # svg_north.add_text_centered_rotated('O', styles['big'], 0, 80, 0)
    svg_north.save_drawing()


# TODO: Fix radial gradient performance. We can just do classes for these and assign them more reasonably than having 10000. Especially relevant when we get to glow if we want to add it.
# TODO: Stars: Gradient should be a parameter. Should probably be defined with a dict up front or list or something.
# TODO: Stars: Glow.
# TODO: Maybe glowing lines.

# TODO: Start labeling! And add DSOs which aer now clean lol.
# TODO: Axes, labels, dates, elliptic, equator, RA/DEC


