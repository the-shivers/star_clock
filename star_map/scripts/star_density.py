"""
The original purpose of this file was to map out the milky way so it could be rendered on a map.
The methodology was going to be looking at star density, segmenting the celestial sphere
into "cells" based on declination and right ascension, then if the star density in that cell
met some cutoff, it would be considered part of the milky way. This went okay, but ultimately 
I found some better milky way outlines that I will be using instead. 
"""

import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import pandas as pd
import numpy as np

# STAR_DATA_LOC = 'star_map/data/athyg_24_reduced_m10.csv'
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
        self.local_star_count = 0
        self.local_area = self.calculate_area()
        self.local_density = None  # To be calculated later with actual star count
        self.cumulative_star_count = 0  # To store the count including its neighbors
        self.cumulative_area = 0  # To store the area including its neighbors
        self.cumulative_density = None  # To be calculated later with cumulative star count
        self.ra_midpoint = (ra_min + ra_max) / 2
        self.dec_midpoint = (dec_min + dec_max) / 2

    def calculate_area(self):
        return (self.ra_max - self.ra_min) * (self.dec_max - self.dec_min) * np.cos(np.radians((self.dec_min + self.dec_max) / 2))

    def update_density(self):
        if self.local_area > 0:
            self.local_density = self.local_star_count / self.local_area

    def update_smoothed_values(self, cumulative_star_count, cumulative_area):
        self.cumulative_star_count = cumulative_star_count
        self.cumulative_area = cumulative_area
        if self.cumulative_area > 0:
            self.cumulative_density = self.cumulative_star_count / self.cumulative_area

    def __repr__(self):
        return (
            f"SkyRegion(RA: {self.ra_min}-{self.ra_max}h, DEC: {self.dec_min}-{self.dec_max}°, "
            f"Local Stars: {self.local_star_count}, Local Area: {self.local_area:.2f}, Local Density: {self.local_density:.2f}, "
            f"Cumulative Stars: {self.cumulative_star_count}, Cumulative Area: {self.cumulative_area:.2f}, Cumulative Density: {self.cumulative_density:.2f})"
        )

class SkyGrid:
    def __init__(self, ra_range, dec_range, ra_divisions, dec_divisions):
        self.ra_divisions = np.linspace(ra_range[0], ra_range[1], ra_divisions + 1)
        self.dec_divisions = np.linspace(dec_range[0], dec_range[1], dec_divisions + 1)
        self.regions = self.create_regions(ra_divisions, dec_divisions)

    def create_regions(self, ra_divisions, dec_divisions):
        regions = []
        for i in range(ra_divisions):
            row = []
            for j in range(dec_divisions):
                ra_min = self.ra_divisions[i]
                ra_max = self.ra_divisions[i + 1]
                dec_min = self.dec_divisions[j]
                dec_max = self.dec_divisions[j + 1]
                region = SkyRegion(ra_min, ra_max, dec_min, dec_max)
                row.append(region)
            regions.append(row)
        return regions

    def assign_stars(self, star_data, verbose=False):
        """Assigns each star to the appropriate region based on its RA and DEC."""
        if verbose:
            total_stars = len(star_data)
            ten_percent_chunk = total_stars // 10  # Determine the size of each 10% chunk of the data
            print("Starting star assignment...")
        for index, star in enumerate(star_data.itertuples(index=True, name='Star')):
            if verbose and index % ten_percent_chunk == 0 and index != 0:
                percent_complete = (index // ten_percent_chunk) * 10
                print(f"{percent_complete}% complete.")
            # Use attribute access for namedtuple
            ra_index = np.searchsorted(self.ra_divisions, star.ra, side='right') - 1
            dec_index = np.searchsorted(self.dec_divisions, star.dec, side='right') - 1
            if 0 <= ra_index < len(self.regions) and 0 <= dec_index < len(self.regions[0]):
                region = self.regions[ra_index][dec_index]
                region.local_star_count += 1
        if verbose:
            print("Star assignment completed.")

    def calculate_densities(self):
        for row in self.regions:
            for region in row:
                region.update_density()
    
    def smooth_densities(self, verbose=False):
        if verbose:
            print("Starting density smoothing...")
            total_regions = len(self.regions)
            ten_percent_chunk = total_regions // 10  # Determine the size of each 10% chunk of the data
        for i, row in enumerate(self.regions):
            if verbose and i % ten_percent_chunk == 0 and i != 0:
                percent_complete = (i // ten_percent_chunk) * 10
                print(f"{percent_complete}% complete.")
            for j, region in enumerate(row):
                cumulative_star_count = 0
                cumulative_area = 0
                # Examine self and neighboring cells
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < len(self.regions) and 0 <= nj < len(self.regions[0]):
                            neighbor = self.regions[ni][nj]
                            cumulative_star_count += neighbor.local_star_count
                            cumulative_area += neighbor.local_area
                # Update the region with cumulative values
                region.update_smoothed_values(cumulative_star_count, cumulative_area)
        print('Density smoothing complete!')

    def __repr__(self):
        return f"SkyGrid(RA Divisions: {len(self.ra_divisions)-1}, DEC Divisions: {len(self.dec_divisions)-1})"


# FUNCTIONS
    
# def plot_sky_regions_density(sky_grid, density_type='local'): # polar
#     fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))
#     axs[0].set_title('Northern Celestial Hemisphere')
#     axs[1].set_title('Southern Celestial Hemisphere')

#     if density_type == 'local':
#         densities = [region.local_density for row in sky_grid.regions for region in row if region.local_density is not None]
#     elif density_type == 'cumulative':
#         densities = [region.cumulative_density for row in sky_grid.regions for region in row if region.cumulative_density is not None]
#     else:
#         raise ValueError("density_type must be 'local' or 'cumulative'")
    
#     cmap = plt.cm.viridis
#     norm = plt.Normalize(min(densities), max(densities))

#     for ax in axs:
#         ax.set_ylim(0, 1)  # Set the limits of the plot to [0, 1]
#         ax.spines['polar'].set_visible(False)  # Remove the polar axes border
    
#     for row in sky_grid.regions:
#         for region in row:
#             theta_start = np.radians((region.ra_min / 24 * 360) + 7.5)  # Shift by half an hour to left align
#             theta_end = np.radians((region.ra_max / 24 * 360) + 7.5)
#             width = theta_end - theta_start
#             r_outer = 1 - (abs(region.dec_max) / 90)
#             r_inner = 1 - (abs(region.dec_min) / 90)
#             ax_index = 0 if region.dec_midpoint >= 0 else 1
#             ax = axs[ax_index]
#             density = region.local_density if density_type == 'local' else region.cumulative_density
#             color = cmap(norm(density))
#             if ax_index == 0:
#                 theta_start, theta_end = 2 * np.pi - theta_end, 2 * np.pi - theta_start  # Correct direction for northern hemisphere
#             ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor=color)
    
#     northern_theta_ticks = np.concatenate([[0], np.linspace(0, 2 * np.pi, 9)[:-1][::-1][:-1]])
#     southern_theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]  # Not reversed
#     theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
#     r_ticks = np.array([1-(2/9), 1-(4/9), 1-(6/9), 1-(8/9)])  # Corrected for 20°, 40°, 60°, 80° declinations
#     r_labels = ["20°", "40°", "60°", "80°"]
#     for i, ax in enumerate(axs):
#         ax.set_theta_zero_location('N')  # Set the 0 degree to the top of the plot
#         if i == 0:
#             ax.set_xticks(northern_theta_ticks)
#             ax.set_xticklabels(theta_labels)
#         else:
#             ax.set_xticks(southern_theta_ticks)
#             ax.set_xticklabels(theta_labels)
#         ax.set_yticks(r_ticks)
#         ax.set_yticklabels(r_labels, fontsize=8)  # Set font size smaller for degree labels
#         for label in ax.get_yticklabels():
#             label.set_rotation('vertical')  # Rotate the degree labels to vertical
#             label.set_color('white')
#     plt.show()

def get_star_data(star_data_loc, data_types):
    """Load star data from a CSV file, retaining all entries."""
    return pd.read_csv(star_data_loc, usecols=data_types.keys(), dtype=data_types)

def plot_contours(sky_grid, density_type='cumulative', levels=4, sigma=1, title='Star Density Contours', extra_points = []):
    # Determine the grid size
    ra_size = len(sky_grid.ra_divisions) - 1
    dec_size = len(sky_grid.dec_divisions) - 1
    
    # Initialize an empty array with NaNs which won't be plotted
    density_array = np.full((dec_size, ra_size), np.nan)
    
    # Map each region's density to the appropriate place in the array
    ra_midpoints = []
    dec_midpoints = []
    for i, row in enumerate(sky_grid.regions):
        ra_midpoints.append(row[0].ra_midpoint)
        for j, region in enumerate(row):
            dec_midpoints.append(region.dec_midpoint)
            # Convert RA and DEC midpoints to array indices
            ra_idx = int((region.ra_midpoint - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
            dec_idx = int((region.dec_midpoint - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
            # Update the array with the region
            # Update the array with the region's density value
            if density_type == 'cumulative':
                density_array[dec_idx, ra_idx] = region.cumulative_density
            elif density_type == 'local':
                density_array[dec_idx, ra_idx] = region.local_density

    # Further smoothed grid
    further_smoothed_grid = gaussian_filter(density_array, sigma=sigma)
    
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_title(title)

    # Generate contour lines
    contour = ax.contour(further_smoothed_grid, levels=levels, cmap='viridis')
    
    plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f', colors='black', manual=False, inline_spacing=5, use_clabeltext=True)
    # Set the tick marks for RA and DEC
    num_ra_ticks = 12
    num_dec_ticks = 12
    ra_indices = np.linspace(0, density_array.shape[1] - 1, num_ra_ticks, dtype=int)
    dec_indices = np.linspace(0, density_array.shape[0] - 1, num_dec_ticks, dtype=int)
    ra_labels = [f'{ra_midpoints[idx]:.1f}h' for idx in ra_indices]
    dec_labels = [f'{dec_midpoints[idx]:.1f}°' for idx in dec_indices]
    plt.xticks(ra_indices, ra_labels)
    plt.yticks(dec_indices, dec_labels)
    plt.xlabel('Right Ascension (RA)')
    plt.ylabel('Declination (DEC)')

    # Plot extra points
    for ra, dec in extra_points:
        ra_idx = int((ra - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
        dec_idx = int((dec - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
        plt.scatter(ra_idx, dec_idx, color='red', s=100, edgecolor='white', zorder=5)  # zorder for layering on top

    plt.tight_layout()
    plt.show()


    
def plot_density_heatmap(sky_grid, density_type='cumulative', extra_points=[]):  # Cartesian, not polar
    # Determine the grid size
    ra_size = len(sky_grid.ra_divisions) - 1
    dec_size = len(sky_grid.dec_divisions) - 1
    
    # Initialize an empty array with NaNs which won't be plotted
    density_array = np.full((dec_size, ra_size), np.nan)
    
    # Map each region's density to the appropriate place in the array
    ra_midpoints = []
    dec_midpoints = []
    for i, row in enumerate(sky_grid.regions):
        ra_midpoints.append(row[0].ra_midpoint)
        for j, region in enumerate(row):
            dec_midpoints.append(region.dec_midpoint)
            ra_idx = int((region.ra_midpoint - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
            dec_idx = int((region.dec_midpoint - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
            if density_type == 'cumulative':
                density_array[dec_idx, ra_idx] = region.cumulative_density
            elif density_type == 'local':
                density_array[dec_idx, ra_idx] = region.local_density
    
    # Plot the density heatmap
    plt.figure(figsize=(10, 6))
    plt.imshow(density_array, cmap='viridis', origin='lower', aspect='auto', interpolation='nearest')
    plt.colorbar(label='Density')

    # Plot extra points
    for ra, dec in extra_points:
        ra_idx = int((ra - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
        dec_idx = int((dec - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
        plt.scatter(ra_idx, dec_idx, color='red', s=100, edgecolor='white', zorder=5)  # zorder for layering on top

    # Set the tick marks for RA and DEC
    num_ra_ticks = 12
    num_dec_ticks = 12
    ra_indices = np.linspace(0, density_array.shape[1] - 1, num_ra_ticks, dtype=int)
    dec_indices = np.linspace(0, density_array.shape[0] - 1, num_dec_ticks, dtype=int)
    ra_labels = [f'{ra_midpoints[idx]:.1f}h' for idx in ra_indices]
    dec_labels = [f'{dec_midpoints[idx]:.1f}°' for idx in dec_indices]
    plt.xticks(ra_indices, ra_labels)
    plt.yticks(dec_indices, dec_labels)
    plt.xlabel('Right Ascension (RA)')
    plt.ylabel('Declination (DEC)')
    plt.title('Sky Region Density Heatmap')
    plt.show()


def plot_cartesian_density_with_contours(sky_grid, density_type='cumulative', levels=4, sigma=1, title='Star Density Heatmap with Contours', extra_points=[]): # Cartesian
    ra_size = len(sky_grid.ra_divisions) - 1
    dec_size = len(sky_grid.dec_divisions) - 1
    density_array = np.full((dec_size, ra_size), np.nan)
    ra_midpoints = []
    dec_midpoints = []
    for i, row in enumerate(sky_grid.regions):
        ra_midpoints.append(row[0].ra_midpoint)
        for j, region in enumerate(row):
            dec_midpoints.append(region.dec_midpoint)
            ra_idx = int((region.ra_midpoint - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
            dec_idx = int((region.dec_midpoint - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
            if density_type == 'cumulative':
                density_array[dec_idx, ra_idx] = region.cumulative_density
            elif density_type == 'local':
                density_array[dec_idx, ra_idx] = region.local_density
    smoothed_density = gaussian_filter(density_array, sigma=sigma)
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_title(title)
    img = ax.imshow(smoothed_density, cmap='viridis', origin='lower', aspect='auto', interpolation='nearest')
    fig.colorbar(img, ax=ax, label='Star Density')
    contours = ax.contour(smoothed_density, levels=levels, colors='white', alpha=0.5, linewidths=0.5)
    # plt.clabel(contours, inline=True, fontsize=8, fmt='%1.0f', colors='black', manual=False, inline_spacing=5, use_clabeltext=True)
    for (ra, dec), name in extra_points:
        ra_idx = int((ra - sky_grid.ra_divisions[0]) / (sky_grid.ra_divisions[-1] - sky_grid.ra_divisions[0]) * ra_size)
        dec_idx = int((dec - sky_grid.dec_divisions[0]) / (sky_grid.dec_divisions[-1] - sky_grid.dec_divisions[0]) * dec_size)
        ax.scatter(ra_idx, dec_idx, color='red', s=15, edgecolor='white', zorder=5)
        ax.text(ra_idx, dec_idx, name, color='white', ha='left', va='center', fontsize=4)
    num_ra_ticks = 12
    num_dec_ticks = 12
    ra_indices = np.linspace(0, smoothed_density.shape[1] - 1, num_ra_ticks, dtype=int)
    dec_indices = np.linspace(0, smoothed_density.shape[0] - 1, num_dec_ticks, dtype=int)
    ra_labels = [f'{ra_midpoints[idx]:.1f}h' for idx in ra_indices]
    dec_labels = [f'{dec_midpoints[idx]:.1f}°' for idx in dec_indices]
    ax.set_xticks(ra_indices)
    ax.set_yticks(dec_indices)
    ax.set_xticklabels(ra_labels)
    ax.set_yticklabels(dec_labels)
    ax.set_xlabel('Right Ascension (RA)')
    ax.set_ylabel('Declination (DEC)')
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'star_map/pics/cartesian_contours_{levels}_{sigma}_{density_type}.png', format='png', dpi=300)


def plot_polar_density_with_contours(sky_grid, density_type='local', levels=4, sigma=3, extra_points=[]): # Polar plot!
    fig, axs = plt.subplots(1, 2, subplot_kw=dict(polar=True), figsize=(12, 6))
    axs[0].set_title('Northern Celestial Hemisphere')
    axs[1].set_title('Southern Celestial Hemisphere')
    # Determine the type of density to plot
    if density_type == 'local':
        densities = [region.local_density for row in sky_grid.regions for region in row if region.local_density is not None]
    elif density_type == 'cumulative':
        densities = [region.cumulative_density for row in sky_grid.regions for region in row if region.cumulative_density is not None]
    else:
        raise ValueError("density_type must be 'local' or 'cumulative'")
    cmap = plt.cm.viridis
    norm = plt.Normalize(min(densities), max(densities))
    theta_resolution = len(sky_grid.regions)
    r_resolution = len(sky_grid.regions[0])
    theta = np.linspace(0, 2 * np.pi, theta_resolution)
    r = np.linspace(0, 1, r_resolution)
    theta_grid, r_grid = np.meshgrid(theta, r)
    # Initialize two separate grids for northern and southern hemispheres
    density_grid_north = np.zeros_like(theta_grid)
    density_grid_south = np.zeros_like(theta_grid)
    for ax in axs:
        ax.set_ylim(0, 1)
        ax.spines['polar'].set_visible(False)
    for row in sky_grid.regions:
        for region in row:
            theta_start = np.radians(region.ra_min / 24 * 360)
            theta_end = np.radians(region.ra_max / 24 * 360)
            width = theta_end - theta_start
            ax_index = 0 if region.dec_midpoint >= 0 else 1
            ax = axs[ax_index]
            density = region.local_density if density_type == 'local' else region.cumulative_density
            if region.dec_midpoint >= 0:
                grid_to_fill = density_grid_north
                r_inner = 1 - (abs(region.dec_max) / 90)
                r_outer = 1 - (abs(region.dec_min) / 90)
                theta_start, theta_end = 2 * np.pi - theta_end, 2 * np.pi - theta_start
            else:
                grid_to_fill = density_grid_south
                r_outer = 1 - (abs(region.dec_max) / 90)
                r_inner = 1 - (abs(region.dec_min) / 90)
            color = cmap(norm(density))
            mask = (theta_grid >= theta_start) & (theta_grid <= theta_end) & (r_grid >= r_inner) & (r_grid <= r_outer)
            grid_to_fill[mask] = density
            ax.bar(theta_start, height=(r_outer - r_inner), width=width, bottom=r_inner, color=color, edgecolor=color)
    for ax, grid in zip(axs, [density_grid_north, density_grid_south]):
        # Apply Gaussian filter if needed
        filtered_grid = gaussian_filter(grid, sigma=sigma)
        contour = ax.contour(theta_grid, r_grid, filtered_grid, levels=levels, colors='white', linewidths=0.5, alpha=0.5,)
        # plt.clabel(contour, inline=True, fontsize=8)
    for (ra, dec), name in extra_points:
        theta = np.radians(ra / 24 * 360)
        if dec >= 0:  # Adjust theta for northern hemisphere stars
            theta = 2 * np.pi - theta
        r = 1 - (abs(dec) / 90)
        ax_index = 0 if dec >= 0 else 1
        ax = axs[ax_index]
        ax.scatter(theta, r, color='red', s=15, edgecolor='white', zorder=5)
        ax.text(theta, r, name, color='white', ha='left', va='center', fontsize=4)
    # Set tick marks for RA and DEC correctly for both hemispheres
    northern_theta_ticks = np.concatenate([[0], np.linspace(0, 2 * np.pi, 9)[:-1][::-1][:-1]])
    southern_theta_ticks = np.linspace(0, 2 * np.pi, 9)[:-1]
    theta_labels = ["{}h".format(h) for h in range(0, 24, 3)]
    r_ticks = np.array([7/9, 5/9,  3/9, 1/9])
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
        ax.set_yticklabels(r_labels)
        for label in ax.get_yticklabels():
            label.set_color('white')
            label.set_rotation('vertical')
    # plt.show()
    plt.savefig(f'star_map/pics/polar_contours_{levels}_{sigma}_{density_type}.png', format='png', dpi=300)


if __name__ == '__main__':
    star_df = get_star_data(STAR_DATA_LOC, DATA_TYPES)
    ra_range = (0, 24)  # Total Right Ascension range
    dec_range = (-90, 90) # Total Declination range
    precision = 10 # 10 is one cell per degree of declination and degree of RA (180 * 360 cells). 1 is chunky cells 10 deg tall and 10 deg wide. (18 * 24 cells)
    ra_divisions = int(4 * 9 * precision)  # Number of divisions along RA
    dec_divisions = int(2 * 9 * precision)  # Number of divisions along DEC
    sky_grid = SkyGrid(ra_range, dec_range, ra_divisions, dec_divisions)
    sky_grid.assign_stars(star_df, verbose=True)
    sky_grid.calculate_densities()
    sky_grid.smooth_densities(verbose=True)
    # plot_sky_regions_density(sky_grid, density_type='cumulative')
    # plot_sky_regions_binary(sky_grid, 1200, 1600)

    stars_to_plot = [
        ((16.49012987, -26.43194608), '  Antares'),

        ((20.69053151, 45.28033423), '  Deneb'),
        ((19.84630057, 8.86738491), '  Altair'),
        ((18.61560722, 38.78299311), '  Vega'),

        ((5.27813768, 45.99902927), '  Capella'),
        ((7.75537884, 28.02631031), '  Pollux'),
        ((7.65514946, 5.22750767), '  Procyon'),
        ((6.7525694, -16.71314306), '  Sirius'),
        ((5.24229757, -8.20163919), '  Rigel'),
        ((4.5986668, 16.50976164), '  Aldebaran'),

        ((5.91952477, 7.40703634), '  Betelgeuse'),
        ((5.67931245, -1.94257841), '  Alnitak'),
        ((5.60355905, -1.20191725), '  Alnilam'),
        ((5.53344367, -0.29908071), '  Mintaka'),

        # ((0, 40), 'Calibration Point')
    ]

    # Contouring
    # ra_midpoints, dec_midpoints = sky_grid.get_midpoints()
    # plot_contours(sky_grid, density_type='cumulative', levels=[1000, 2000, 3200], sigma=1, title='Star Density Contours', extra_points = [(16.5, -26.4), (19.8, 8.9), (20.7, 45.3), (18.6, 38.8)])
    # plot_density_heatmap(sky_grid, density_type='cumulative', extra_points = [(16.5, -26.4), (19.8, 8.9), (20.7, 45.3), (18.6, 38.8)])
    # plot_density_heatmap_with_contours(sky_grid, density_type='cumulative', levels=[1050, 2400], sigma=1, title='Star Density Heatmap with Contours', extra_points=[(16.5, -26.4), (19.8, 8.9), (20.7, 45.3), (18.6, 38.8)])
    # plot_density_heatmap_with_contours(sky_grid, density_type='cumulative', levels=[1050, 2400], sigma=1, title='Star Density Heatmap with Contours', extra_points=[(16.5, -26.4), (19.8, 8.9), (20.7, 45.3), (18.6, 38.8)])
    # plot_polar_density_with_contours(sky_grid, density_type='local', levels=[1300, 1800, 2800], sigma=1.2, extra_points=stars_to_plot)
    plot_cartesian_density_with_contours(sky_grid, density_type='local', levels=[1200, 1800, 2800], sigma=1, extra_points=stars_to_plot)