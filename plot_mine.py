#!/usr/bin/python3
import numpy as np
from mayavi.mlab import *
import sys
import os

global args
args = sys.argv

BLOCK_ID = 0
X_VAL = 1
Y_VAL = 2
Z_VAL = 3

x = []
y = []
z = []

values = []
max_val = 0

plot_blocks = []

count = -1

with open(args[1], 'r') as infile:
    for line in infile:
        count = count + 1
        split_line = line.split();
        # plot_blocks.append(count)
        if int(split_line[Z_VAL]) > 12 and int(split_line[Z_VAL]) < 16:
            plot_blocks.append(count)
        x.append(int(split_line[X_VAL]))
        y.append(int(split_line[Y_VAL]))
        z.append(int(split_line[Z_VAL]))

with open(args[2], 'r') as infile:
    for data in infile:
        if int(data) > max_val:
            max_val = int(data)
        values.append(int(data))

values = [-1 if x==max_val else x for x in values]

prev_max = max(values)

values = [x/prev_max for x in values]

nodes = points3d([x[i] for i in plot_blocks],[y[i] for i in plot_blocks] , [z[i] for i in plot_blocks], mode='cube', scale_factor=1.0)
nodes.glyph.scale_mode = 'scale_by_vector'

nodes.mlab_source.dataset.point_data.scalars = [values[i] for i in plot_blocks]

show()