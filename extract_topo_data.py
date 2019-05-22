#!/usr/bin/env python3

import sys
import os
import datetime

from netCDF4 import Dataset

import numpy as np

def process(filename):
    print('processing file: {}'.format(filename))
    file_base, file_extension = os.path.splitext(filename)

    nc_file = Dataset(filename, mode='r')

    scale = 0.0001

    x = nc_file.variables['x'][:] * scale
    y = nc_file.variables['y'][:] * scale
    z = nc_file.variables['topographic__elevation'][:][0] * scale

    output_file = "{}.obj".format(file_base)

    length1 = len(x)
    length2 = len(x[0])

    with open(output_file, 'w') as f:
        f.write('# Creation time: {}\n'.format(datetime.datetime.now()))
        f.write('o terrain\n')

        for i in range(0, length1):
            for j in range(0, length2):
                f.write('v {} {} {}\n'.format(x[i][j], y[i][j], z[i][j]))

        f.write('s off\n')

        for i in range(0, length1 - 1):
            for j in range(0, length1  - 1):
                f1 = (i * length1) + j + 1
                f2 = (i * length1) + j + 2
                f3 = ((i + 1) * length1) + j + 1
                f4 = ((i + 1) * length1) + j + 2

                f.write('f {} {} {}\n'.format(f1, f2, f3))
                f.write('f {} {} {}\n'.format(f2, f4, f3))


for filename in sys.argv[1:]:
    process(filename)
