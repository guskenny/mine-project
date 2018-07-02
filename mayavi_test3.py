#!/usr/bin/python3
# Imports
import numpy as np
from mayavi import mlab

# Primitives
N = 100 # Number of points

ones = np.ones(N)
scalars = np.arange(N) # Key point: set an integer for each point

# # Define color table (including alpha), which must be uint8 and [0,255]
# colors = (np.random.random((N, 4))*255).astype(np.uint8)
# colors[:,-1] = 255 # No transparency


# number of discrete colours
n = 10

t = np.linspace(-510, 510, n)                                              
RGB = np.round(np.clip(np.dstack([t, 510-np.abs(t), -t, 255*np.ones(n)]), 0, 255)).astype(np.uint8)[0]

# print(RGB)

idxes = np.random.randint(n, size=N)

colors = RGB[idxes]
# print(colors)



# Define coordinates and points
x, y, z = (np.random.random(N)*100).astype(np.uint8),(np.random.random(N)*100).astype(np.uint8),(np.random.random(N)*100).astype(np.uint8)  # Assign x, y, z values to match color

@mlab.animate(delay=1000)
def anim():
    f = mlab.gcf()
    while True:
        for i in range(N):
            plot_points = [(i+j)%N for j in range(10)]
            # print("x: ",x[plot_points],", y: ",y[plot_points], ", z: ",z[plot_points], ", s: ", scalars[plot_points])
            pts.mlab_source.reset(x=x[plot_points], y=y[plot_points], z=z[plot_points],scalars = scalars[plot_points])
            pts.module_manager.scalar_lut_manager.lut.table = colors[plot_points]   
            yield

pts = mlab.points3d(x, y, z, scalars, mode='sphere',scale_factor=10.0) # Create points
pts.glyph.color_mode = 'color_by_scalar' # Color by scalar
pts.glyph.scale_mode = 'scale_by_vector'
# Set look-up table and redraw
pts.module_manager.scalar_lut_manager.lut.table = colors


base_pts = mlab.points3d(x, y, z, scalars, mode='sphere',opacity=0.2,scale_factor=10.0) # Create points
base_pts.glyph.color_mode = 'color_by_scalar' # Color by scalar
base_pts.glyph.scale_mode = 'scale_by_vector'
# Set look-up table and redraw
base_pts.module_manager.scalar_lut_manager.lut.table = colors


anim()
mlab.show()