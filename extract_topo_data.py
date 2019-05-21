#!/usr/bin/env python3

import sys
import os

from netCDF4 import Dataset

import numpy as np

def process(filename):
    print('processing file: {}'.format(filename))
    file_base, file_extension = os.path.splitext(filename)

    nc_file = Dataset(filename, mode='r')

    x = nc_file.variables['x'][:]
    y = nc_file.variables['y'][:]
    z = nc_file.variables['topographic__elevation'][:][0]

    output_file = "{}.txt".format(file_base)

    length1 = len(x)
    length2 = len(x[0])

    with open(output_file, 'w') as f:
        for i in range(0, length1):
            for j in range(0, length2):
                f.write('{}, {}, {}\n'.format(x[i][j], y[i][j], z[i][j]))

for filename in sys.argv[1:]:
    process(filename)
