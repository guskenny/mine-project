#!/usr/bin/python3
import numpy as np
from mayavi import mlab
import sys
import os

global args
args = sys.argv

BLOCK_ID = 0
X_VAL = 0
Y_VAL = 1
Z_VAL = 2

coords = [[],[],[]]

values = []
max_data = 0

plot_blocks = []
no_plot = []

with open(args[1], 'r') as infile:
    for line in infile:
        split_line = line.split();
        coords[X_VAL].append(int(split_line[X_VAL+1]))
        coords[Y_VAL].append(int(split_line[Y_VAL+1]))
        coords[Z_VAL].append(int(split_line[Z_VAL+1]))

with open(args[2], 'r') as infile:
    for data in infile:
        if int(data) > max_data:
            max_data = int(data)
        values.append(int(data))

periods = [x for x in values]

no_plot = [0 for x in range(len(periods))]

for b in range(len(periods)):
    if periods[b] == max_data:
        no_plot[b] = 1

values = [0 if x==max_data else x for x in values]

prev_max = max(values)

max_coords = [max(coords[X_VAL]), max(coords[Y_VAL]), max(coords[Z_VAL])]

scalars = np.arange(len(coords[X_VAL]))

t = np.linspace(-510, 510, prev_max+1)                                              
RGB = np.round(np.clip(np.dstack([t, 510-np.abs(t), -t, 255*np.ones(prev_max+1)]), 0, 255)).astype(np.uint8)[0]

colours = RGB[values]

values = [x/prev_max for x in values]

# base_plt.glyph.color_mode = 'color_by_scalar' # Color by scalar
# base_plt.module_manager.scalar_lut_manager.lut.table = colours

fig = mlab.figure(size=(800,800))

@mlab.animate(delay=1000)
def anim():
    f = mlab.gcf()
    while True:
        for p in range(prev_max): 
            plot_blocks.clear()
            for b in range(len(periods)):
                if no_plot[b] == 0 and periods[b] == p:
                    plot_blocks.append(b)

            # print("period: ",p,", prev_max: ",prev_max,", plot_blocks: ",len(plot_blocks))    

            # mlab.draw(figure=f)

            plt.mlab_source.reset(x=[coords[X_VAL][i] for i in plot_blocks], y=[coords[Y_VAL][i] for i in plot_blocks], z=[coords[Z_VAL][i] for i in plot_blocks],scalars=scalars[plot_blocks])   

            plt.module_manager.scalar_lut_manager.lut.table = RGB[len(plot_blocks)*[p]]

            # plt.mlab_source.dataset.point_data.scalars = [values[i] for i in plot_blocks]

            yield


for b in range(len(coords[Z_VAL])):
    if no_plot[b] == 0 and coords[Z_VAL][b] <= max_coords[Z_VAL] and coords[Z_VAL][b] > max_coords[Z_VAL]-3:
        plot_blocks.append(b)

base_plt=mlab.points3d([coords[X_VAL][i] for i in plot_blocks], [coords[Y_VAL][i] for i in plot_blocks], [coords[Z_VAL][i] for i in plot_blocks], mode='2dsquare', scale_factor=1.0,opacity=0.1)
base_plt.glyph.scale_mode = 'scale_by_vector'
base_plt.mlab_source.dataset.point_data.scalars = [values[i] for i in plot_blocks]

plot_blocks.clear()
for b in range(len(coords[Z_VAL])):
    if no_plot[b] == 0 and periods[b] == 0:
        plot_blocks.append(b)


plt=mlab.points3d([coords[X_VAL][i] for i in plot_blocks], [coords[Y_VAL][i] for i in plot_blocks], [coords[Z_VAL][i] for i in plot_blocks],scalars[plot_blocks], mode='cube', scale_factor=1.0)
plt.glyph.scale_mode = 'scale_by_vector'
plt.glyph.color_mode = 'color_by_scalar' # Color by scalar
plt.module_manager.scalar_lut_manager.lut.table = RGB[len(plot_blocks)*[0]]

anim()
mlab.show()