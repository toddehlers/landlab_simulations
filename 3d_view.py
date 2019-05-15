#!/usr/bin/env python

import sys
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from netCDF4 import Dataset

import numpy as np

from PIL import Image, ImageChops

font_color = 'white'
bg_color = 'black'

mpl.rcParams['text.color'] = font_color
mpl.rcParams['axes.labelcolor'] = font_color
mpl.rcParams['xtick.color'] = font_color
mpl.rcParams['ytick.color'] = font_color
mpl.rcParams['savefig.facecolor'] = bg_color

# Clean up file names:
# mmv 'output*__*.nc' 'output_#2.nc'
#
# python3 ../../../3d_view.py 0.8 azucar/*.nc
#
# Create movie:
# ffmpeg -framerate 10 -i output_02%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -s 800x600 output.mp4
# ffmpeg -framerate 10 -pattern_type glob -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -s 800x600 output.mp4


def process(filename):
    print('processing file: {}'.format(filename))
    file_base, file_extension = os.path.splitext(filename)

    nc_file = Dataset(filename, mode='r')

    scale_factor = 0.001

    x = nc_file.variables['x'][:] * scale_factor
    y = nc_file.variables['y'][:] * scale_factor
    z = nc_file.variables['topographic__elevation'][:][0] * scale_factor

    # print('x: {}, y: {}, z: {}'.format(len(x), len(y), len(z)))

    z_min = 0.0

    # For Pan de Azucar: 0.8
    # For Nahuelbuta: 2.0
    # For La Campana: 3.0
    # For Santa Gracia: 1.2
    z_max = 3.0

    num_of_vals = 200

    z_scale = 10.0

    fig = plt.figure(figsize=(10, 7))
    fig.patch.set_facecolor(bg_color)
    fig.set_facecolor(bg_color)

    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(0.0, z_scale)
    ax.view_init(40, 225)
    ax.set_xlabel('X [km]')
    ax.set_ylabel('Y [km]')
    ax.set_zlabel('')
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor(bg_color)
    ax.yaxis.pane.set_edgecolor(bg_color)
    ax.zaxis.pane.set_edgecolor(bg_color)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_facecolor(bg_color)
    ax.set_zticks([])

    surface = ax.plot_surface(x, y, z, cmap=cm.terrain, linewidth=0, antialiased=True, vmin=z_min, vmax=z_max, rcount=num_of_vals, ccount=num_of_vals)
    cbar_axes = fig.add_axes([0.85, 0.12, 0.03, 0.4])
    cbar = fig.colorbar(surface, label='Elevation [km]', cax=cbar_axes)

    image_file = '{}.png'.format(file_base)
    plt.savefig(image_file, dpi=100, bbox_inches='tight')
    plt.close()

    im = Image.open(image_file)
    w, h = im.size
    im = im.crop((50, 200, w, h))
    im = ImageChops.offset(im, 0, -25)
    im.save(image_file)


for filename in sys.argv[1:]:
    process(filename)
