"""
The original purpose of this file was to map out the milky way so it could be rendered on a map.
The methodology was going to be looking at star density, segmenting the celestial sphere
into "cells" based on declination and right ascension, then if the star density in that cell
met some cutoff, it would be considered part of the milky way. This went okay, but ultimately 
I found some better milky way outlines that I will be using instead. 
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

STAR_DATA_LOC = 'star_map/data/athyg_v24.csv'
DATA_TYPES = {
    'hip': 'Int64', # Nullable integer, for stars with no HIPPARCOS ID
    'proper': str,
    'ra': float,
    'dec': float,
    'mag': float,
    'ci': float
}

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
    
    def smooth_densities(self):
        """Applies a convolutional smoothing to the densities across the grid."""
        smoothed_density_grid = np.zeros_like(self.region_grid, dtype=float)
        ra_size, dec_size = self.region_grid.shape
        total_steps = ra_size  # Total number of rows to process
        print("Starting density smoothing...")
        for i in range(ra_size):
            for j in range(dec_size):
                total_stars = 0
                total_area = 0
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < ra_size and 0 <= nj < dec_size:
                            neighbor = self.region_grid[ni, nj]
                            # Approximate area for each region, considering small region and cosine approximation for declination
                            area = (neighbor.ra_max - neighbor.ra_min) * \
                                (neighbor.dec_max - neighbor.dec_min) * \
                                np.cos(np.radians((neighbor.dec_min + neighbor.dec_max) / 2))
                            total_stars += neighbor.star_count
                            total_area += area
                if total_area > 0:
                    smoothed_density_grid[i, j] = total_stars / total_area
            # Print progress after each row is processed
            progress_percentage = ((i + 1) / total_steps) * 100
            print(f"Density smoothing progress: {progress_percentage:.2f}% complete.")
        print("Density smoothing completed.")
        return smoothed_density_grid
    
    def summarize_densities(self, n=3):
        """Prints a summary of the top n and bottom n regions by density."""
        region_densities = [(region, region.calculate_density()) for region in self.regions]
        sorted_regions = sorted(region_densities, key=lambda x: x[1], reverse=True)
        top_regions = sorted_regions[:n]
        bottom_regions = sorted_regions[-n:]
        summary_str = "Top {} densest regions:\n".format(n)
        for region, density in top_regions:
            summary_str += "{}, Density: {:.2f}\n".format(region, density)
        summary_str += "\nBottom {} least dense regions:\n".format(n)
        for region, density in bottom_regions:
            summary_str += "{}, Density: {:.2f}\n".format(region, density)
        print(summary_str)
    
    def __repr__(self):
        return f"SkyGrid(RA Divisions: {len(self.ra_divisions)-1}, DEC Divisions: {len(self.dec_divisions)-1}, Total Regions: {len(self.regions)})"


# FUNCTIONS

def plot_sky_regions_density(sky_grid):
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
        if ax_index == 0:
            theta_start = 2 * np.pi - theta_start
            theta_end = 2 * np.pi - theta_end
        ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor=color)
        print(region)
    # Configure the theta ticks (RA in hours) and radial ticks (Declination in degrees)
    northern_theta_ticks = np.concatenate([[0], np.linspace(0, 2 * np.pi, 9)[:-1][::-1][:-1]])
    southern_theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]  # Not reversed
    theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
    r_ticks = np.array([1-(2/9), 1-(4/9), 1-(6/9), 1-(8/9)])  # Corrected for 20°, 40°, 60°, 80° declinations
    r_labels = ["20°", "40°", "60°", "80°"]
    for i, ax in enumerate(axs):
        ax.set_theta_zero_location('N')  # Set the 0 degree to the top of the plot
        if i == 0:
            ax.set_xticks(northern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        else:
            ax.set_xticks(southern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        ax.set_yticks(r_ticks)
        ax.set_yticklabels(r_labels, fontsize=8)  # Set font size smaller for degree labels
        for label in ax.get_yticklabels():
            label.set_rotation('vertical')  # Rotate the degree labels to vertical
            label.set_color('white')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # plt.colorbar(sm, ax=axs, orientation='horizontal', pad=0.1, aspect=30, label='Star Density')
    plt.show()


def plot_sky_regions_smoothed_density(sky_grid, density_grid, output_path = 'test.png'):
    print("Producing smoothed figure.")
    fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))
    axs[0].set_title('Northern Celestial Hemisphere')
    axs[1].set_title('Southern Celestial Hemisphere')
    cmap = plt.cm.viridis
    norm = plt.Normalize(np.min(density_grid), np.max(density_grid))
    ra_size, dec_size = sky_grid.region_grid.shape
    for ax in axs:
        ax.set_ylim(0, 1)
        ax.spines['polar'].set_visible(False)
    for i in range(ra_size):
        for j in range(dec_size):
            region = sky_grid.region_grid[i, j]
            theta_start = np.radians((region.ra_min / 24 * 360) + 7.5)  # Shift by half an hour to left align
            theta_end = np.radians((region.ra_max / 24 * 360) + 7.5)
            width = theta_end - theta_start
            r_outer = 1 - (abs(region.dec_max) / 90)
            r_inner = 1 - (abs(region.dec_min) / 90)
            ax_index = 0 if region.dec_min >= 0 else 1
            ax = axs[ax_index]
            density = density_grid[i, j]
            color = cmap(norm(density))
            if ax_index == 0:
                theta_start = 2 * np.pi - theta_start
                theta_end = 2 * np.pi - theta_end
            ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor=color)
    # Configure the theta ticks (RA in hours) and radial ticks (Declination in degrees)
    northern_theta_ticks = np.concatenate([[0], np.linspace(0, 2 * np.pi, 9)[:-1][::-1][:-1]])
    southern_theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]  # Not reversed
    theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
    r_ticks = np.array([1-(2/9), 1-(4/9), 1-(6/9), 1-(8/9)])  # Corrected for 20°, 40°, 60°, 80° declinations
    r_labels = ["20°", "40°", "60°", "80°"]
    for i, ax in enumerate(axs):
        ax.set_theta_zero_location('N')  # Set the 0 degree to the top of the plot
        if i == 0:
            ax.set_xticks(northern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        else:
            ax.set_xticks(southern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        ax.set_yticks(r_ticks)
        ax.set_yticklabels(r_labels, fontsize=8)  # Set font size smaller for degree labels
        for label in ax.get_yticklabels():
            label.set_rotation('vertical')  # Rotate the degree labels to vertical
            label.set_color('white')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # plt.colorbar(sm, ax=axs, orientation='horizontal', pad=0.1, aspect=30, label='Star Density')
    # plt.show()
    print('Saving figure!')
    plt.savefig(output_path, format='png', dpi=300)
    plt.close(fig)


def plot_sky_regions_binary(sky_grid, cutoff_density1, cutoff_density2): # 2 is higher
    fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))
    axs[0].set_title('Northern Celestial Hemisphere')
    axs[1].set_title('Southern Celestial Hemisphere')
    for ax in axs:
        ax.set_ylim(0, 1)
        ax.spines['polar'].set_visible(False)
    for region in sky_grid.regions:
        theta_start = np.radians((region.ra_min / 24 * 360) + 7.5)  # Shift by half an hour to left align
        theta_end = np.radians((region.ra_max / 24 * 360) + 7.5)
        width = theta_end - theta_start
        r_outer = 1 - (abs(region.dec_max) / 90)
        r_inner = 1 - (abs(region.dec_min) / 90)
        ax_index = 0 if region.dec_min >= 0 else 1
        ax = axs[ax_index]
        density = region.calculate_density()
        color = 'blue' if density > cutoff_density1 else 'white'
        color = 'green' if density > cutoff_density2 else color
        if ax_index == 0:
            theta_start = 2 * np.pi - theta_start
            theta_start = 2 * np.pi - theta_end
        ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor=color)
        print(region)
    northern_theta_ticks = np.concatenate([[0], np.linspace(0, 2 * np.pi, 9)[:-1][::-1][:-1]])
    southern_theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]  # Not reversed
    theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
    r_ticks = np.array([1-(2/9), 1-(4/9), 1-(6/9), 1-(8/9)])  # Corrected for 20°, 40°, 60°, 80° declinations
    r_labels = ["20°", "40°", "60°", "80°"]
    for i, ax in enumerate(axs):
        ax.set_theta_zero_location('N')
        if i == 0:
            ax.set_xticks(northern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        else:
            ax.set_xticks(southern_theta_ticks)
            ax.set_xticklabels(theta_labels)
        ax.set_yticks(r_ticks)
        ax.set_yticklabels(r_labels, fontsize=8)
        for label in ax.get_yticklabels():
            label.set_rotation('vertical')
            label.set_color('white')
    plt.show()


def get_star_data(star_data_loc, data_types):
    """Load star data from a CSV file, retaining all entries."""
    return pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)


# MAIN LOGIC

if __name__ == '__main__':
    star_df = get_star_data(STAR_DATA_LOC, DATA_TYPES)
    ra_range = (0, 24)  # Total Right Ascension range
    dec_range = (-90, 90) # Total Declination range
    precision = 15 # 10 is one cell per degree of declination and degree of RA (180 * 360 cells). 1 is chunky cells 10 deg tall and 10 deg wide. (18 * 24 cells)
    ra_divisions = int(4 * 9 * precision)  # Number of divisions along RA
    dec_divisions = int(2 * 9 * precision)  # Number of divisions along DEC
    sky_grid = SkyGrid(ra_range, dec_range, ra_divisions, dec_divisions)
    sky_grid.assign_stars(star_df, verbose=True)
    # plot_sky_regions_density(sky_grid)
    # plot_sky_regions_binary(sky_grid, 1200, 1600)
    smoothed = sky_grid.smooth_densities()
    plot_sky_regions_smoothed_density(sky_grid, smoothed, f'{precision}_precision_smoothed.png')