import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

class SkyRegion:
    def __init__(self, ra_min, ra_max, dec_min, dec_max):
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
        self.star_count = 0  # Directly track the count of stars

    def __repr__(self):
        return f"SkyRegion(RA: {self.ra_min}-{self.ra_max}h, DEC: {self.dec_min}-{self.dec_max}°, Stars: {self.star_count})"

    def calculate_density(self):
        """Calculate and return the density of stars in this region."""
        # Approximate area calculation considering small region approximation
        area = (self.ra_max - self.ra_min) * (self.dec_max - self.dec_min) * np.cos(np.radians((self.dec_min + self.dec_max) / 2))
        return self.star_count / area

class SkyGrid:
    def __init__(self, ra_range, dec_range, ra_divisions, dec_divisions):
        self.regions = []
        self.ra_divisions = np.linspace(ra_range[0], ra_range[1], ra_divisions + 1)
        self.dec_divisions = np.linspace(dec_range[0], dec_range[1], dec_divisions + 1)
        for i in range(ra_divisions):
            for j in range(dec_divisions):
                ra_min = self.ra_divisions[i]
                ra_max = self.ra_divisions[i+1]
                dec_min = self.dec_divisions[j]
                dec_max = self.dec_divisions[j+1]
                self.regions.append(SkyRegion(ra_min, ra_max, dec_min, dec_max))
        # Create a grid to access regions by their index directly
        self.region_grid = np.array(self.regions).reshape((ra_divisions, dec_divisions))

    def assign_stars(self, star_data, verbose=False):
        """Assigns each star to the appropriate region based on its RA and DEC."""
        total_stars = len(star_data)
        ten_percent_chunk = total_stars // 10  # Determine the size of each 10% chunk of the data
        if verbose:
            print("Starting star assignment...")
        for index, star in enumerate(star_data.itertuples(index=True, name='Star')):
            # Print progress at each 10% interval
            if verbose and index % ten_percent_chunk == 0 and index != 0:
                percent_complete = (index // ten_percent_chunk) * 10
                print(f"{percent_complete}% complete.")
            # Extracting RA and DEC from the star tuple
            ra = getattr(star, 'ra')
            dec = getattr(star, 'dec')
            ra_index = np.searchsorted(self.ra_divisions, ra, side='right') - 1
            dec_index = np.searchsorted(self.dec_divisions, dec, side='right') - 1
            if 0 <= ra_index < len(self.ra_divisions) - 1 and 0 <= dec_index < len(self.dec_divisions) - 1:
                self.region_grid[ra_index, dec_index].star_count += 1
        if verbose:
            print("Star assignment completed.")

    def calculate_densities(self):
        """Calculate densities for all regions."""
        return [region.calculate_density() for row in self.region_grid for region in row]
    
    def summarize_densities(self, n=3):
        """Prints a summary of the top n and bottom n regions by density."""
        # Calculate densities for all regions and pair them with regions
        region_densities = [(region, region.calculate_density()) for region in self.regions]
        # Sort the regions by density
        sorted_regions = sorted(region_densities, key=lambda x: x[1], reverse=True)
        # Select top n and bottom n regions
        top_regions = sorted_regions[:n]
        bottom_regions = sorted_regions[-n:]
        # Create a formatted string for output
        summary_str = "Top {} densest regions:\n".format(n)
        for region, density in top_regions:
            summary_str += "{}, Density: {:.2f}\n".format(region, density)
        summary_str += "\nBottom {} least dense regions:\n".format(n)
        for region, density in bottom_regions:
            summary_str += "{}, Density: {:.2f}\n".format(region, density)
        # Print the summary
        print(summary_str)
    
    def __repr__(self):
        return f"SkyGrid(RA Divisions: {len(self.ra_divisions)-1}, DEC Divisions: {len(self.dec_divisions)-1}, Total Regions: {len(self.regions)})"


# FUNCTIONS

def get_star_data(star_data_loc, data_types):
    """Load star data from a CSV file, retaining all entries."""
    return pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)

def get_stars_dict(star_data, mag_limit):
    """Generate a dictionary for stars with a HIP number within a magnitude limit."""
    filtered_data = star_data.dropna(subset=['hip'])
    filtered_data = filtered_data[filtered_data['mag'] <= mag_limit]
    return {int(row['hip']): Star(int(row['hip']), row['proper'], row['ra'], row['dec'], row['mag'], row['ci'])
        for row in filtered_data.to_dict('records')}

def plot_sky_regions_density(sky_grid):
    """
    Rough plot of sky grid densities. Note that there are some issues here:
    - RA is in reverse for the northern celestial hemisphere.
    - Degree tickmarks are at an ugly slant.
    
    But it gets the job done.
    """
    fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))
    axs[0].set_title('Northern Celestial Hemisphere')
    axs[1].set_title('Southern Celestial Hemisphere')
    cmap = plt.cm.viridis
    norm = plt.Normalize(min(sky_grid.calculate_densities()), max(sky_grid.calculate_densities()))
    for ax in axs:
        ax.set_ylim(0, 1)  # Set the limits of the plot to [0, 1]
        ax.spines['polar'].set_visible(False)  # Remove the polar axes border
    for region in sky_grid.regions:
        theta_start = np.radians((region.ra_min / 24 * 360) + 7.5)  # Shift by half an hour to left align
        theta_end = np.radians((region.ra_max / 24 * 360) + 7.5)
        width = theta_end - theta_start
        r_outer = 1 - (abs(region.dec_max) / 90)
        r_inner = 1 - (abs(region.dec_min) / 90)
        ax_index = 0 if region.dec_min >= 0 else 1
        ax = axs[ax_index]
        density = region.calculate_density()
        color = cmap(norm(density))
        ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor='none')
        print(region)
    # Configure the theta ticks (RA in hours) and radial ticks (Declination in degrees)
    theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]  # Exclude the last tick to avoid duplication of 0h and 24h
    theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
    r_ticks = np.array([1-(2/9), 1-(4/9), 1-(6/9), 1-(8/9)])  # Corrected for 20°, 40°, 60°, 80° declinations
    r_labels = ["20°", "40°", "60°", "80°"]
    for ax in axs:
        ax.set_theta_zero_location('N')  # Set the 0 degree to the top of the plot
        ax.set_xticks(theta_ticks)
        ax.set_xticklabels(theta_labels)
        ax.set_yticks(r_ticks)
        ax.set_yticklabels(r_labels, fontsize=8)  # Set font size smaller for degree labels
        for label in ax.get_yticklabels():
            label.set_rotation('vertical')  # Rotate the degree labels to vertical
            label.set_color('white')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plt.colorbar(sm, ax=axs, orientation='horizontal', pad=0.1, aspect=30, label='Star Density')
    plt.show()




if __name__ == '__main__':
    star_df = get_star_data(STAR_DATA_LOC, DATA_TYPES)
    stars_dict = get_stars_dict(star_df, MAG_LIMIT)
    # parser = ConstellationParser(CONSTELLATIONS_LOC, stars_dict)
    # constellations = parser.parse()
    # constellationship = Constellationship(constellations, 'rey')
    ra_range = (0, 24)  # Total Right Ascension range
    dec_range = (-90, 90) # Total Declination range
    precision = 5 # 10 is one cell per degree of declination and degree of RA (180 * 360 cells). 1 is chunky cells 10 deg tall and 10 deg wide. (18 * 24 cells)
    ra_divisions = 4 * 9 * precision  # Number of divisions along RA
    dec_divisions = 2 * 9 * precision  # Number of divisions along DEC
    sky_grid = SkyGrid(ra_range, dec_range, ra_divisions, dec_divisions)
    sky_grid.assign_stars(star_df, verbose=True)
    # densities = sky_grid.calculate_densities()
    sky_grid.summarize_densities(n=10)
    plot_sky_regions_density(sky_grid)




    

