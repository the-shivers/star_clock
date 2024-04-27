import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

MAG_LIMIT = 5
COLS = ['hip', 'con', 'proper', 'ra', 'dec', 'mag', 'ci']
STAR_DATA_LOC = 'star_clock/star_map/data/star_databases/athyg_v24.csv'

data = pd.read_csv(STAR_DATA_LOC)
data = data[pd.notnull(data['hip'])]
data['hip'] = data['hip'].astype(int)
data = data[data['mag'] <= MAG_LIMIT]
data = data[COLS]
data['theta'] = np.radians(data['ra'] * 15)
data['r'] = abs(data['dec'] / 90)
north = data[data['dec'] >= 0]

















# def plot_stereographic_projection(data, compress_factor=1):
#     # Cartesian star coords
#     x = data['r'] * np.cos(data['theta'])
#     y = data['r'] * np.sin(data['theta'])
    
#     # Create the plot
#     plt.figure(figsize=(8, 8))
#     ax = plt.gca()
#     ax.scatter(x, y, color='white', s=1)  # 's' controls the size of the points
#     ax.set_facecolor('black')
    
#     # Label points if the 'proper' column has a name
#     for i, txt in enumerate(data['proper']):
#         if pd.notnull(txt):
#             ax.annotate(txt, (x[i], y[i]), textcoords="offset points", xytext=(0,5), ha='center', color='white', fontsize=8)
    
#     # Add gridlines for specific declinations
#     decs = [20, 40, 60, 80]  # Declination values for gridlines
#     for dec in decs:
#         dec_rad = np.radians(dec)
#         r = np.sqrt(2 * np.arctan(np.cos(dec_rad))) * compress_factor
#         circle = plt.Circle((0, 0), r, color='white', fill=False, linestyle='--', linewidth=0.5)
#         ax.add_patch(circle)
    
#     ax.set_xlim(-np.sqrt(2)*compress_factor, np.sqrt(2)*compress_factor)
#     ax.set_ylim(-np.sqrt(2)*compress_factor, np.sqrt(2)*compress_factor)
#     ax.set_aspect('equal', adjustable='box')  # Keep the aspect ratio of the plot as a circle
    
#     plt.xticks([])
#     plt.yticks([])  # Remove axis marks
    
#     plt.title('Stereographic Projection of Stars')
#     plt.show()


# # Assuming 'filtered_data' is already loaded and contains the 'r' and 'theta' columns
# data = pd.read_csv('star_clock/star_map/data/star_databases/athyg_v24.csv')
# filtered_data = get_hemisphere_disk(data, 'north', 2.3, 5)
# plot_stereographic_projection(filtered_data)

